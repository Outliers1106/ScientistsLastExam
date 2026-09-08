"""Regression checks for BatchEffectDiscovery."""
from __future__ import annotations
import importlib.util
import sys
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
BIO=ROOT/"benchmarks/Biology"

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module

TASKS=(('BatchEffectDiscovery', 'analyze_expression', 'reference_analysis.py'),)

class BiologyExpansionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.loaded={}
        for task,entry,refname in TASKS:
            ev=load(BIO/task/"verification/evaluator.py","ev_"+task)
            base=load(BIO/task/"solution.py","base_"+task)
            ref=load(BIO/task/"verification"/refname,"ref_"+task)
            cls.loaded[task]=(ev,entry,ev.evaluate(getattr(base,entry)),ev.evaluate(getattr(ref,entry)))

    def test_withdrawn_family_scores_are_valid_and_bounded(self):
        for task,(_,_,baseline,reference) in self.loaded.items():
            self.assertEqual(baseline["valid"],1.0,task)
            self.assertAlmostEqual(baseline["combined_score"],0.0,places=12,msg=task)
            self.assertEqual(reference["valid"],1.0,task)
            self.assertGreater(reference["combined_score"],0.0,task)
            self.assertLessEqual(reference["combined_score"],1.0,task)

    def test_evaluators_are_deterministic(self):
        for task,entry,refname in TASKS:
            ev=load(BIO/task/"verification/evaluator.py","det_ev_"+task)
            ref=load(BIO/task/"verification"/refname,"det_ref_"+task)
            one=ev.evaluate(getattr(ref,entry)); two=ev.evaluate(getattr(ref,entry))
            self.assertEqual(one,two,task)

    def test_discovery_metrics_do_not_depend_on_process_hash_seed(self):
        script="""
import importlib.util, json, pathlib, sys
def load(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module
root=pathlib.Path(sys.argv[1]); entry=sys.argv[2]
ev=load(root/'verification/evaluator.py','ev')
base=load(root/'solution.py','base')
ref=load(root/'verification'/sys.argv[3],'ref')
print(json.dumps([ev.evaluate(getattr(base,entry)),ev.evaluate(getattr(ref,entry))],sort_keys=True))
"""
        for task,entry,ref in TASKS:
            results=[]
            for seed in ("1","2","7"):
                env={**os.environ,"PYTHONHASHSEED":seed,"OPENBLAS_NUM_THREADS":"1"}
                results.append(subprocess.check_output(
                    [sys.executable,"-c",script,str(BIO/task),entry,ref],env=env,text=True,timeout=30))
            self.assertEqual(results,[results[0]]*len(results),task)

    def test_malformed_candidates_fail_closed(self):
        bad=(lambda *args:{},lambda *args:None,lambda *args:"bad",lambda *args:[],
             lambda *args:{"unexpected":1},lambda *args:float("nan"),
             lambda *args:(_ for _ in ()).throw(RuntimeError("boom")))
        for task,(ev,_,_,_) in self.loaded.items():
            for candidate in bad:
                metrics=ev.evaluate(candidate)
                self.assertEqual(metrics["valid"],0.0,task)
                self.assertEqual(metrics["combined_score"],0.0,task)

    def test_discovery_tasks_publish_axes_and_denominators(self):
        keys={"development_mechanism_score","development_false_discovery_rate",
              "development_false_discovery_count","development_claim_count",
              "development_correct_refusal_rate","development_discovery_coverage"}
        for task in ('BatchEffectDiscovery',):
            self.assertTrue(keys<=set(self.loaded[task][3]),task)

    def test_blanket_abstention_scores_zero(self):
        ev=self.loaded['BatchEffectDiscovery'][0]
        self.assertEqual(ev.evaluate(lambda p, callback: {"discoveries":[],"abstain":True,"reason_code":"not_identifiable"})["combined_score"],0.)

    def test_batch_scores_effects_and_reason_codes(self):
        ev=self.loaded["BatchEffectDiscovery"][0]
        ref=load(BIO/"BatchEffectDiscovery"/"verification"/"reference_analysis.py","batch_mut_ref")
        def zero_effect(problem,measure):
            out=ref.analyze_expression(problem,measure)
            out["discoveries"]=[{**row,"effect":0.0} for row in out["discoveries"]]
            return out
        def wrong_reason(problem,measure):
            out=ref.analyze_expression(problem,measure)
            out["reason_code"]="supported"
            return out
        reference=self.loaded["BatchEffectDiscovery"][3]
        self.assertLess(ev.evaluate(zero_effect)["combined_score"],reference["combined_score"])
        self.assertLess(ev.evaluate(wrong_reason)["combined_score"],reference["combined_score"])

    def test_batch_confounded_world_has_no_cross_cell(self):
        ev=self.loaded["BatchEffectDiscovery"][0]; lab=ev._Lab("confounded",103)
        self.assertEqual(lab.available,[[0,0],[1,1]])
        with self.assertRaises(ValueError): lab(1,0)

    def test_batch_sampling_layout_does_not_disclose_effect_presence(self):
        ev=self.loaded["BatchEffectDiscovery"][0]
        layouts=[]
        for kind in ("supported","null","confounded"):
            rows=ev._initial(kind,101)
            layouts.append([(r["batch"],r["condition"],r["library_size"]) for r in rows])
        self.assertEqual(layouts[0],layouts[1])
        self.assertEqual(layouts[0],layouts[2])
        self.assertEqual(ev._Lab("supported",101).available,ev._Lab("null",101).available)

    def test_batch_no_discovery_policies_score_zero(self):
        ev=self.loaded["BatchEffectDiscovery"][0]
        def deny(p,measure):
            return {"discoveries":[],"abstain":False,"reason_code":"no_effect"}
        def layout_only(p,measure):
            confounded=len(p["available_cells"])<3
            return {"discoveries":[],"abstain":confounded,
                    "reason_code":"not_identifiable" if confounded else "no_effect"}
        for candidate in (deny,layout_only):
            self.assertEqual(ev.evaluate(candidate)["combined_score"],0.0)

    def test_task_entrypoints_start_without_pythonpath_from_another_directory(self):
        env=dict(os.environ)
        env.pop("PYTHONPATH",None)
        with tempfile.TemporaryDirectory() as tmp:
            for task,_,_ in TASKS:
                result=subprocess.run([sys.executable,"-I",str(BIO/task/"frontier_eval/run_eval.py"),"--help"],
                                      cwd=tmp,env=env,capture_output=True,text=True,timeout=20)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertIn("--candidate",result.stdout)


def test_single_invalid_instance_zeros_all_aggregate_scores():
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "benchmarks/Biology/BatchEffectDiscovery/verification/evaluator.py"
    spec = importlib.util.spec_from_file_location("invalid_instance_oracle", path)
    oracle = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(oracle)
    reference_spec = importlib.util.spec_from_file_location("reference", path.parents[1] / 'verification/reference_analysis.py')
    reference = importlib.util.module_from_spec(reference_spec)
    reference_spec.loader.exec_module(reference)
    calls = 0
    def malformed(*args):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise ValueError("one invalid instance")
        return reference.analyze_expression(*args)
    result = oracle.evaluate(malformed)
    assert calls > 1
    worlds = result.get("per_world", result.get("per_instance"))
    assert not worlds[0]["valid"] and all(row["valid"] for row in worlds[1:])
    assert result["valid"] == result["combined_score"] == 0
    assert all(value == 0 for key, value in result.items()
               if key.startswith("heldout_") and "score" in key)


def test_supported_gene_counts_vary_across_seeds():
    ev = load(BIO/"BatchEffectDiscovery/verification/evaluator.py", "gene_counts_ev")
    counts = {np.count_nonzero(ev._truth("supported", seed)[0]) for seed in range(100, 200)}
    assert counts == set(range(2, 9))
    for kind in ("null", "confounded"):
        assert not np.any(ev._truth(kind, 101)[0])


def test_heldout_no_discovery_policy_scores_zero():
    ev = load(BIO/"BatchEffectDiscovery/verification/evaluator.py", "heldout_floor_ev")
    def policy(p, measure):
        confounded = len(p["available_cells"]) < 3
        return dict(discoveries=[], abstain=confounded,
                    reason_code="not_identifiable" if confounded else "no_effect")
    result = ev.evaluate(policy)
    assert result["combined_score"] == result["heldout_scientific_score"] == 0
    assert result["development_refusal_world_count"] == result["heldout_refusal_world_count"] == 1


def test_fdr_counts_false_claims_in_supported_worlds_too():
    ev = load(BIO/"BatchEffectDiscovery/verification/evaluator.py", "fdr_ev")
    ref = load(BIO/"BatchEffectDiscovery/verification/reference_analysis.py", "fdr_ref")
    result = ev.evaluate(ref.analyze_expression)
    for prefix, indices in (("development", ev.DEV), ("heldout", ev.HELD)):
        rows = [result["per_world"][i] for i in indices]
        false = sum(row["false_positive"] for row in rows)
        claims = sum(row["claimed"] for row in rows)
        assert result[f"{prefix}_false_discovery_count"] == false
        assert result[f"{prefix}_claim_count"] == claims
        assert result[f"{prefix}_false_discovery_rate"] == false/max(1, claims)
    assert any(row["kind"] == "supported" and row["false_positive"] > 0 for row in result["per_world"])


def test_withdrawal_ablation_batch_column_still_has_no_effect():
    ev = load(BIO/"BatchEffectDiscovery/verification/evaluator.py", "ablation_ev")
    ref = load(BIO/"BatchEffectDiscovery/verification/reference_analysis.py", "ablation_ref")
    def without_batch(p, measure):
        if len(p["available_cells"]) < 3:
            return dict(discoveries=[], abstain=True, reason_code="not_identifiable")
        rows = list(p["initial_samples"])
        for i in range(p["sample_budget"]):
            rows.append(measure(*((0, 1), (1, 0))[i % 2]))
        x = np.array([[1, row["condition"]] for row in rows], float)
        y = np.log1p([row["counts"] for row in rows])
        effect = np.linalg.lstsq(x, y, rcond=None)[0][1]
        hits = [dict(gene=g, effect=float(e)) for g, e in zip(p["gene_ids"], effect) if abs(e) > .55]
        return dict(discoveries=hits, abstain=False, reason_code="supported" if hits else "no_effect")
    full = ev.evaluate(ref.analyze_expression)
    ablated = ev.evaluate(without_batch)
    assert abs(full["combined_score"]-ablated["combined_score"]) < 1e-12
    assert abs(full["heldout_scientific_score"]-ablated["heldout_scientific_score"]) < 1e-12
