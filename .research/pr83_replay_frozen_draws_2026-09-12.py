"""Replay every available frozen first proposal twice, retaining invalid results."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

def sha(data):
    return hashlib.sha256(data).hexdigest()

def canonical(data):
    return (json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()

def read(path):
    return json.loads(path.read_text())

def write(path, data):
    with path.open("x") as stream:
        json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    path.chmod(0o600)

os.umask(0o077)
root, model_plan_path, private = map(lambda p: Path(p).resolve(), sys.argv[1:])
sys.path.insert(0, str(root))
from scripts.task_campaign import read_plan
from sle.algorithms.common import runtime_source_sha256, task_package_sha256
from sle.registry import find_task
from sle.evaluate import evaluate_candidate
plan = read_plan(model_plan_path)
git = lambda *a: subprocess.check_output(["git", "-c", "core.commitGraph=false", "-C", str(root), *a], text=True).strip()
assert not git("status", "--porcelain")
current_revision = git("rev-parse", "HEAD")
assert root not in private.parents and private != root
assert all(os.environ.get(k) == "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"))
spec = find_task(plan["bindings"]["task_id"], include_uncertified=True)
assert runtime_source_sha256() == plan["bindings"]["runtime_source_sha256"]
task_path = spec.task_dir.relative_to(root).as_posix()
changed = git("diff", "--name-only", plan["source_provenance"]["git_revision"], "HEAD", "--", task_path).splitlines()
assert set(changed) <= {task_path + "/TASK_CARD.yaml", task_path + "/references/known_best.md"}
current_package = task_package_sha256(spec)
private.mkdir(mode=0o700, exist_ok=False)
frozen = []
unavailable = []
for cell in plan["cells"]:
    folder = Path(cell["workdir"])
    events = [json.loads(line) for line in (folder / "trajectory.jsonl").read_text().splitlines() if line.strip()]
    if len(events) != 2 or events[1].get("candidate_sha256") is None:
        unavailable.append(cell["seed"])
        continue
    event = events[1]
    code = [x for x in read(folder / "checkpoint.json")["evaluated_candidates"] if x["step"] == 1]
    assert len(code) == 1
    candidate_bytes = code[0]["program"].encode()
    assert sha(candidate_bytes) == event["candidate_sha256"]
    receipt_id = event["algorithm_metadata"]["evaluation_request_id"]
    receipt = read(folder / "evaluation_ledger/receipts" / (receipt_id + ".json"))
    assert receipt["metrics"] == event["metrics"] and receipt["metrics_sha256"] == sha(canonical(event["metrics"]))
    candidate = private / ("candidate_%s.py" % cell["seed"])
    with candidate.open("xb") as stream:
        stream.write(candidate_bytes)
    candidate.chmod(0o600)
    frozen.append((cell["seed"], candidate, event["candidate_sha256"], receipt["metrics"]))
replay_plan = {"schema_version":1,"source_revision":current_revision,"original_bindings":plan["bindings"],
               "current_task_package_sha256":current_package,"documentation_files_changed":changed,
               "original_model_plan_sha256":plan["plan_sha256"],"timeout_s":plan["timeout_s"],
               "model_calls":0,"planned_evaluations":len(frozen)*2,"unavailable_seed_labels":unavailable,
               "candidates":[{"seed_label":seed,"candidate_sha256":digest,"original_metrics_sha256":sha(canonical(metrics))}
                             for seed,_,digest,metrics in frozen],"driver_sha256":sha(Path(__file__).read_bytes())}
write(private / "plan.json", replay_plan)
rows = []
for seed, candidate, candidate_hash, original_metrics in frozen:
    for repeat in range(2):
        started = time.monotonic()
        metrics = evaluate_candidate(spec, candidate, timeout_s=plan["timeout_s"])
        raw = private / ("metrics_%s_%s.json" % (seed,repeat))
        write(raw,metrics)
        row = {"seed_label":seed,"repeat":repeat,"candidate_sha256":candidate_hash,
               "full_metrics_sha256":sha(canonical(metrics)),"raw_file_sha256":sha(raw.read_bytes()),
               "matches_original_full_metrics":metrics == original_metrics,"seconds":time.monotonic()-started,
               "aggregate_metrics":{k:v for k,v in metrics.items() if type(v) in (float,int,bool)},
               "observed_worlds":len(metrics.get("per_instance",[])),
               "valid_worlds":sum(bool(r.get("valid")) for r in metrics.get("per_instance",[])),
               "candidate_failure_kind":metrics.get("candidate_failure_kind")}
        rows.append(row)
        write(private / ("receipt_%s_%s.json" % (seed,repeat)),row)
        print(seed,repeat,metrics.get("valid"),metrics.get("combined_score"),row["matches_original_full_metrics"],flush=True)
unchanged = not git("status", "--porcelain") and git("rev-parse", "HEAD") == replay_plan["source_revision"]
assert unchanged and task_package_sha256(spec) == current_package
report = {"schema_version":1,"scope":"fixed retained proposals; no retry, candidate selection, or new model call",
          "plan":replay_plan,"source_unchanged":unchanged,"records":rows,
          "every_available_candidate_metrics_match_original":all(x["matches_original_full_metrics"] for x in rows),
          "all_three_available":len(frozen)==3,"unavailable_seed_labels":unavailable}
write(private / "public-replay-review.json",report)
