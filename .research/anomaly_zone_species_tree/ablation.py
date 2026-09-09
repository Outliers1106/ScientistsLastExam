import os as _os, tempfile as _tempfile
TASK = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "benchmarks", "Biology", "AnomalyZoneSpeciesTree"))
OUT = _tempfile.gettempdir()
import sys, json, importlib.util, collections, numpy as np
sys.path.insert(0, TASK + "/verification")
import evaluator as ev, msc
SRC = open(TASK + "/verification/reference_quartet_consensus.py").read()
BASE = open(TASK + "/solution.py").read()
FULL = msc.FULL_MASK
def compatible(a, b): return (a & b) == 0 or (a & b) == a or (a & b) == b or (a | b) == FULL
def variant(name, patches, inject=None):
    src = SRC
    for old, new in patches:
        assert old in src, (name, old)
        src = src.replace(old, new)
    ns = {"_oracle_splits": lambda locus: msc.splits_of(msc.parse_newick(locus["nj_tree"])[1])}
    if inject: ns.update(inject)
    exec(compile(src, "variant_" + name, "exec"), ns)
    return ns["infer_species_tree"]
def greedy_from_counts(counts):
    # greedy consensus needs split frequencies, not quartet counts; handled by patching the loop instead
    raise NotImplementedError
CORR = "splits = _splits(_neighbour_joining(_gamma_distances(seqs, shape)))"
rows = []
variants = [
 ("reference", []),
 ("oracle NJ trees, no gamma correction", [(CORR, "splits = _oracle_splits(locus)")]),
 ("slow loci", [('CHOSEN_CLASS = "medium"', 'CHOSEN_CLASS = "slow"')]),
 ("fast loci", [('CHOSEN_CLASS = "medium"', 'CHOSEN_CLASS = "fast"')]),
 ("300-site loci (all 300 in the catalogue)", [("CHOSEN_SITES = 800", "CHOSEN_SITES = 300")]),
 ("2000-site loci (96 of them)", [("CHOSEN_SITES = 800", "CHOSEN_SITES = 2000")]),
 ("never refusing", [("REFUSAL_STATISTIC = 115.0", "REFUSAL_STATISTIC = 1e9")]),
 ("refusal threshold 80", [("REFUSAL_STATISTIC = 115.0", "REFUSAL_STATISTIC = 80.0")]),
 ("refusal threshold 160", [("REFUSAL_STATISTIC = 115.0", "REFUSAL_STATISTIC = 160.0")]),
 ("constant branch lengths 0.1", [("lengths[mask] = float(np.median(estimates)) if estimates else MIN_LENGTH", "lengths[mask] = 0.1")]),
 ("half the budget", [('budget = int(problem["locus_budget"])', 'budget = int(problem["locus_budget"]) // 2')]),
 ("quarter of the budget", [('budget = int(problem["locus_budget"])', 'budget = int(problem["locus_budget"]) // 4')]),
 ("oracle trees and fast loci", [(CORR, "splits = _oracle_splits(locus)"), ('CHOSEN_CLASS = "medium"', 'CHOSEN_CLASS = "fast"')]),
 ("oracle trees, never refusing", [(CORR, "splits = _oracle_splits(locus)"), ("REFUSAL_STATISTIC = 115.0", "REFUSAL_STATISTIC = 1e9")]),
]
# greedy consensus variant: patch the quartet argmax to a greedy consensus of split frequencies
GREEDY_PATCH = [
 ("        for q, quartet in enumerate(QUARTETS):\n            topology = _quartet_topology(splits, quartet)",
  "        _SPLIT_FREQ.update(splits.keys())\n        for q, quartet in enumerate(QUARTETS):\n            topology = _quartet_topology(splits, quartet)"),
 ("    best = int(np.argmax(score))\n", "    best = _greedy_index(_SPLIT_FREQ, split_sets)\n    _SPLIT_FREQ.clear()\n"),
]
def _greedy_index(freq, split_sets):
    chosen = []
    for m, _ in sorted(freq.items(), key=lambda kv: (-kv[1], kv[0])):
        if all(compatible(m, c) for c in chosen):
            chosen.append(m)
        if len(chosen) == 5: break
    target = set(chosen)
    for i, s in enumerate(split_sets):
        if set(s) == target: return i
    return 0
variants.append(("greedy consensus instead of quartets", GREEDY_PATCH))
inject = {"_SPLIT_FREQ": collections.Counter(), "_greedy_index": _greedy_index}
for name, patches in variants:
    fn = variant(name.replace(" ", "_"), patches, inject)
    m = ev.evaluate(fn)
    rows.append((name, m))
    print("%-42s dev %.3f held %.3f topo %.2f bl %.2f fdr %.2f refusal %.2f cov %.2f loci %.0f" % (
        name, m["combined_score"], m["heldout_mechanism_score"], m["development_topology_rate"], m["development_branch_length_score"],
        m["development_false_discovery_rate"], m["development_correct_refusal_rate"], m["development_discovery_coverage"], m["development_mean_loci_bought"]), flush=True)
ns = {}; exec(compile(BASE, "baseline", "exec"), ns); m = ev.evaluate(ns["infer_species_tree"]); rows.append(("baseline: fast long loci, concatenation, never refusing", m))
print("%-42s dev %.3f held %.3f topo %.2f bl %.2f fdr %.2f refusal %.2f cov %.2f" % ("baseline", m["combined_score"], m["heldout_mechanism_score"], m["development_topology_rate"], m["development_branch_length_score"], m["development_false_discovery_rate"], m["development_correct_refusal_rate"], m["development_discovery_coverage"]))
m = ev.evaluate(lambda p, s: {"verdict": "reticulate", "confidence": 1.0}); rows.append(("declining everything", m))
print("declining everything dev %.3f held %.3f fdr %.2f" % (m["combined_score"], m["heldout_mechanism_score"], m["development_false_discovery_rate"]))
json.dump([(n, {k: v for k, v in m.items() if k != "per_instance"}) for n, m in rows], open(OUT + "/anomaly_zone_ablation.json", "w"), indent=1)
