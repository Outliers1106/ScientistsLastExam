"""Discovery-contract pins for AnomalyZoneSpeciesTree.

The public score is mechanism, normalised so that declining every world is exactly zero and so is
naming one fixed tree everywhere. A reticulate world is the unsupported case: declining it is
correct, declining everything is not; a tree world is determinable and its tree is the answer.

Two of these tests pin what the construction checkpoints found the hard way: that every tree world
really is in the anomaly zone, and that the sequencing budget is not free.
"""
from __future__ import annotations

import collections
import importlib.util
import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TASK = ROOT / "benchmarks/Biology/AnomalyZoneSpeciesTree"

FIXED_TREE = "((A:1,B:1):0.1,(C:1,D:1):0.1,((E:1,F:1):0.1,(G:1,H:1):0.1):0.1);"


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AnomalyZoneSpeciesTreeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evaluator = _load(TASK / "verification/evaluator.py", "anomaly_zone_oracle")
        cls.baseline = _load(TASK / "solution.py", "anomaly_zone_baseline")
        cls.reference = _load(TASK / "verification/reference_quartet_consensus.py",
                              "anomaly_zone_reference")

    def test_blanket_refusal_in_both_forms_and_a_fixed_tree_all_score_zero(self):
        for submission in ({"verdict": "reticulate", "confidence": 1.0}, {"abstain": True}):
            metrics = self.evaluator.evaluate(lambda _p, _s, s=submission: dict(s))
            self.assertEqual(metrics["valid"], 1.0, submission)
            self.assertEqual(metrics["combined_score"], 0.0, submission)
            self.assertEqual(metrics["development_correct_refusal_rate"], 1.0, submission)
        fixed = self.evaluator.evaluate(
            lambda _p, _s: {"verdict": "tree", "newick": FIXED_TREE, "confidence": 0.9})
        self.assertEqual(fixed["valid"], 1.0)
        self.assertEqual(fixed["combined_score"], 0.0)
        self.assertEqual(fixed["development_false_discovery_rate"], 1.0)

    def test_every_tree_world_is_in_the_anomaly_zone(self):
        """The one property the task is built on. If the species tree were the most frequent gene
        tree, majority vote and concatenation would work and the task would be a different task."""
        ev = self.evaluator
        msc = sys.modules[ev.msc.__name__]
        for spec in ev.DEVELOPMENT_WORLDS + ev.HELDOUT_WORLDS:
            if spec["kind"] != "anomaly":
                continue
            world = ev._world(spec)
            rng = np.random.default_rng(12345)
            counts = collections.Counter()
            for _ in range(20000):
                gene_tree, _ = msc.simulate_gene_tree(world["trees"][0], rng)
                counts[frozenset(gene_tree.unrooted_splits())] += 1
            truth = frozenset(world["truth_splits"])
            top, top_count = counts.most_common(1)[0]
            self.assertNotEqual(top, truth, spec["seed"])
            self.assertLess(counts[truth], top_count, spec["seed"])

    def test_the_sequencing_budget_is_not_free(self):
        """A quarter of the budget must be materially worse than all of it."""
        source = (TASK / "verification/reference_quartet_consensus.py").read_text(encoding="utf-8")
        quarter = source.replace('budget = int(problem["locus_budget"])',
                                 'budget = int(problem["locus_budget"]) // 4')
        self.assertNotEqual(quarter, source)
        namespace = {}
        exec(compile(quarter, "quarter_budget_reference", "exec"), namespace)  # noqa: S102
        full = self.evaluator.evaluate(self.reference.infer_species_tree)["combined_score"]
        starved = self.evaluator.evaluate(namespace["infer_species_tree"])["combined_score"]
        self.assertLess(starved, full - 0.3, "a quarter of the budget scores as well as all of it")

    def test_the_concatenation_baseline_is_valid_and_scores_zero(self):
        baseline = self.evaluator.evaluate(self.baseline.infer_species_tree)
        reference = self.evaluator.evaluate(self.reference.infer_species_tree)
        self.assertEqual(baseline["valid"], 1.0)
        self.assertEqual(reference["valid"], 1.0)
        self.assertEqual(baseline["combined_score"], 0.0)
        self.assertGreater(baseline["development_false_discovery_rate"], 0.8)
        self.assertEqual(baseline["development_correct_refusal_rate"], 0.0)
        self.assertGreater(reference["combined_score"], 0.5)
        self.assertLess(reference["combined_score"], 1.0)
        # The reference names every tree world and declines every reticulate one, so what it
        # leaves on the table is the branch lengths: it inverts the quartet frequencies without
        # correcting for gene-tree estimation error. If that ever reaches the ceiling the task has
        # stopped measuring the axis it was built around.
        self.assertEqual(reference["development_discovery_coverage"], 1.0)
        self.assertEqual(reference["development_correct_refusal_rate"], 1.0)
        self.assertEqual(reference["development_false_discovery_rate"], 0.0)
        self.assertLess(reference["development_branch_length_score"], 0.9)
        for row in reference["per_instance"]:
            self.assertEqual(row["abstained"], row["kind"] == "reticulate", row)

    def test_two_evaluations_are_identical(self):
        first = self.evaluator.evaluate(self.reference.infer_species_tree)
        second = self.evaluator.evaluate(self.reference.infer_species_tree)
        self.assertEqual(json.dumps(first, sort_keys=True), json.dumps(second, sort_keys=True))

    def test_malformed_candidates_score_zero_without_raising(self):
        def raises(_p, _s):
            raise RuntimeError("boom")

        def overspend(problem, sequence):
            for row in problem["catalogue"]:
                sequence(int(row["locus"]))
            return {"verdict": "reticulate"}

        shapes = {
            "raises": raises,
            "none": lambda _p, _s: None,
            "empty": lambda _p, _s: {},
            "string": lambda _p, _s: "tree",
            "verdict_bad": lambda _p, _s: {"verdict": "network"},
            "tree_no_newick": lambda _p, _s: {"verdict": "tree"},
            "newick_int": lambda _p, _s: {"verdict": "tree", "newick": 7},
            "newick_garbage": lambda _p, _s: {"verdict": "tree", "newick": "((A,B"},
            "missing_taxon": lambda _p, _s: {"verdict": "tree", "newick": FIXED_TREE.replace("H:1", "A:1")},
            "no_lengths": lambda _p, _s: {"verdict": "tree", "newick": "((A,B),(C,D),((E,F),(G,H)));"},
            "negative_length": lambda _p, _s: {"verdict": "tree", "newick": FIXED_TREE.replace(":0.1,((E", ":-0.1,((E")},
            "nan_length": lambda _p, _s: {"verdict": "tree", "newick": FIXED_TREE.replace("G:1", "G:nan")},
            "unresolved": lambda _p, _s: {"verdict": "tree", "newick": "(A:1,B:1,C:1,D:1,E:1,F:1,G:1,H:1);"},
            "confidence_nan": lambda _p, _s: {"verdict": "tree", "newick": FIXED_TREE, "confidence": float("nan")},
            "confidence_str": lambda _p, _s: {"verdict": "tree", "newick": FIXED_TREE, "confidence": "high"},
            "overspend": overspend,
            "bad_locus": lambda _p, s: s(9999),
            "float_locus": lambda _p, s: s(1.0),
        }
        self.assertGreaterEqual(len(shapes), 12)
        for name, candidate in shapes.items():
            metrics = self.evaluator.evaluate(candidate)
            self.assertEqual(metrics["valid"], 0.0, name)
            self.assertEqual(metrics["combined_score"], 0.0, name)
            self.assertEqual(metrics["feasibility_rate"], 0.0, name)

    def test_hidden_axes_stay_out_of_the_search_view(self):
        from sle.metric_visibility import SEARCH_VISIBLE_KEYS

        for key in SEARCH_VISIBLE_KEYS:
            self.assertNotIn("heldout", key)
            self.assertNotIn("mechanism", key)
            self.assertNotIn("topology", key)


if __name__ == "__main__":
    unittest.main()
