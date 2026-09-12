import importlib.util
from pathlib import Path
import sys
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
TASK = ROOT / 'benchmarks/ComputerScience/JuntaVariableDiscovery'


def oracle():
    spec = importlib.util.spec_from_file_location('junta_admission_oracle', TASK / 'verification/evaluator.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_false_discovery_uses_valid_claim_denominator():
    ev = oracle()
    rows = []
    for kind, valid, abstained, false in [
        ('supported', True, False, False), ('supported', True, False, True),
        ('unsupported', True, True, False), ('unsupported', False, True, False),
    ]:
        rows.append(dict(kind=kind, valid=valid, abstained=abstained, false_discovery=false,
                         mechanism_score=0.0, sharpness=0.0, correct_refusal=valid and abstained,
                         confidence_calibration_score=0.0, queries_used=0))
    summary = ev._split_summary(rows)
    assert summary['false_discovery_rate'] == 0.5
    assert summary['false_discovery_count'] == 1
    assert summary['false_discovery_denominator'] == 2
    assert summary['all_world_false_claim_fraction'] == 0.25
    empty = ev._split_summary([dict(r, abstained=True, false_discovery=False) for r in rows])
    assert empty['false_discovery_denominator'] == 0
    assert empty['false_discovery_rate'] == 0.0


@pytest.mark.parametrize('failed_split', ['development', 'heldout'])
def test_one_failed_world_cannot_produce_a_valid_task(failed_split):
    ev = oracle()
    original = ev._evaluate_world

    def one_failure(candidate, spec, split, index):
        row = original(candidate, spec, split, index)
        if split == failed_split and index == 0:
            row.update(valid=False, mechanism_score=0.0, correct_refusal=False)
        return row

    with patch.object(ev, '_evaluate_world', side_effect=one_failure):
        result = ev.evaluate(lambda p, q: {'verdict': 'no_junta'})
    assert result['valid'] == 0.0
    assert result['combined_score'] == 0.0


@pytest.mark.skipif(sys.platform != 'linux', reason='real CandidateProxy requires Linux')
def test_real_proxy_resets_world_state_and_keeps_state_across_queries(tmp_path):
    from sle.secure_eval import CandidateProxy
    ev = oracle()
    source = '''import os
calls = 0
def identify(problem, query):
    global calls
    calls += 1
    marker = '/tmp/junta_world_marker'
    stale = calls != 1 or hasattr(os, '_junta_marker') or os.path.exists(marker)
    os._junta_marker = 1
    with open(marker, 'w') as handle:
        handle.write('same world')
    for _ in range(2):
        answers = query([[0] * problem['n']])
        if len(answers) != 1 or calls != 1 or os._junta_marker != 1:
            return None
        with open(marker) as handle:
            if handle.read() != 'same world':
                return None
    return None if stale else {'verdict': 'no_junta'}
'''
    candidate = tmp_path / 'candidate.py'
    candidate.write_text(source)
    with CandidateProxy(candidate, 'identify', timeout_s=60) as proxy:
        with patch.object(proxy, 'reset_session', wraps=proxy.reset_session) as reset:
            result = ev.evaluate(proxy)
            assert reset.call_count == 18
    assert result['valid'] == 1.0
    assert all(row['valid'] and row['queries_used'] == 2 for row in result['per_instance'])
