"""Clock-sync isolation and numerically equivalent LP coordinates, without a full witness run."""
import importlib.util
import platform
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
TASK = ROOT / "benchmarks/ComputerScience/ClockSyncInversion"


def load(name):
    spec = importlib.util.spec_from_file_location("clock_runtime_" + name, TASK / "verification" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_lp_translation_restores_both_optimizer_and_objective():
    ref = load("reference_envelope_lp")
    horizon = 14400.0
    center = np.array([0.003, -30e-6])
    scale = np.array([1e-6, 1e-6 / horizon])
    a = np.array([[1, horizon], [1, 0], [0, horizon]])
    b = a @ center + np.array([2e-6, 1.5e-6, 1.5e-6])
    c = np.array([-1, -0.5 * horizon])
    bounds = [(float(x), None) for x in center]
    result = ref._solve_centered_lp(c, a, b, bounds, center, scale)
    # In residual microseconds the optimum is the corner (1.5, .5), not
    # (0, 0). A missing c@center term would silently corrupt reported intervals.
    assert result.status == 0
    np.testing.assert_allclose((result.x - center) / scale, [1.5, 0.5], atol=1e-5)
    assert result.fun == pytest.approx(float(c @ center - 1.75e-6), abs=1e-11)


def test_unbounded_lp_is_not_a_model_refusal():
    ref = load("reference_envelope_lp")
    with pytest.raises(RuntimeError, match="optimum or infeasibility"):
        ref._solve_centered_lp(np.array([-1.0]), np.array([[-1.0]]), np.array([0.0]),
                               [(0, None)], np.array([0.0]), np.array([1e-6]))


def test_infeasible_lp_remains_distinct_from_numerical_failure():
    ref = load("reference_envelope_lp")
    result = ref._solve_centered_lp(np.array([1.0]), np.array([[1.0], [-1.0]]),
                                    np.array([0.0, -1e-6]), [(None, None)],
                                    np.array([0.003]), np.array([1e-6]))
    assert result.status == 2


def test_candidate_reset_covers_every_world_and_split_boundary():
    oracle = load("evaluator")

    class Proxy:
        calls = 0
        resets = 0

        def reset_session(self):
            self.resets += 1

        def __call__(self, problem, exchange, wait):
            assert self.resets == self.calls
            self.calls += 1
            return {"verdict": "no_model"}

    candidate = Proxy()
    metrics = oracle.evaluate(candidate)
    count = len(oracle.DEVELOPMENT_WORLDS) + len(oracle.HELDOUT_WORLDS)
    assert candidate.calls == count and candidate.resets == count - 1
    assert all(row["valid"] for row in metrics["per_instance"])


@pytest.mark.skipif(platform.system() != "Linux", reason="requires Linux candidate sandbox")
def test_real_sandbox_worlds_cannot_reuse_globals_or_tmpfs(tmp_path):
    from sle.evaluate import evaluate_candidate
    from sle.spec import load_task_spec
    candidate = tmp_path / "candidate.py"
    candidate.write_text('''from pathlib import Path
seen = False
def identify(problem, exchange, wait):
    global seen
    marker = Path('/tmp/clock_previous_world')
    if seen or marker.exists():
        return []
    seen = True
    marker.write_text('visited')
    return {'verdict': 'no_model'}
''')
    metrics = evaluate_candidate(load_task_spec(TASK), candidate, timeout_s=300)
    assert metrics["valid"] == 1 and metrics["combined_score"] == 0
    assert len(metrics["per_instance"]) == 18
    assert all(row["valid"] for row in metrics["per_instance"])
    assert metrics["development_valid_count"] == 12
    assert metrics["heldout_valid_count"] == 6
