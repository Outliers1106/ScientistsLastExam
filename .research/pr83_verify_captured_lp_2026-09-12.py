"""Independent fixed captured-LP check; private data never enters the checkout."""
import bisect
import hashlib
import importlib.util
import json
import os
from fractions import Fraction
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import linprog

os.umask(0o077)
root, raw, out = map(Path, sys.argv[1:])
expected = "988052ec0c41008a99f081a60a31fbde63d3640b5ad12d2b63a9dc1ad4386855"
assert hashlib.sha256(raw.read_bytes()).hexdigest() == expected
spec = importlib.util.spec_from_file_location("current_reference", root / "benchmarks/ComputerScience/ClockSyncInversion/verification/reference_envelope_lp.py")
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)
data = json.loads(raw.read_text())
A, b, c, center, scale = (np.array(data[key]) for key in ("A", "b", "c", "center", "scale"))
bounds = data["bounds"]
clock_count = int(np.sum(scale < 1e-6))
kept = ref._clock_row_indices(A, b, clock_count)
assert len(kept) < len(A)
groups = {}
for index in kept:
    row = A[index]
    theta, rate, prop = row[:clock_count], row[clock_count:2*clock_count], row[2*clock_count:]
    active = np.flatnonzero(theta)
    if not len(active):
        continue
    t = float(rate[active[0]] / theta[active[0]])
    if not np.array_equal(rate, theta * t):
        continue
    key = tuple(theta) + tuple(prop)
    groups.setdefault(key, []).append((Fraction(t), Fraction(float(b[index]))))
for points in groups.values():
    points.sort()
proven = 0
for index in sorted(set(range(len(A))) - set(kept)):
    row = A[index]
    theta, rate, prop = row[:clock_count], row[clock_count:2*clock_count], row[2*clock_count:]
    active = np.flatnonzero(theta)
    t = Fraction(float(rate[active[0]] / theta[active[0]]))
    points = groups[tuple(theta) + tuple(prop)]
    times = [p[0] for p in points]
    right = bisect.bisect_left(times, t)
    limit = Fraction(float(b[index]))
    if right < len(points) and times[right] == t:
        assert points[right][1] <= limit
    else:
        assert 0 < right < len(points)
        t0, b0 = points[right-1]
        t1, b1 = points[right]
        weight = (t-t0)/(t1-t0)
        assert 0 <= weight <= 1
        assert (1-weight)*b0 + weight*b1 <= limit
    proven += 1

objective = c * scale
objective_scale = max(float(np.max(np.abs(objective))), 1e-12)
scaled_A = A * scale / 1e-6
scaled_b = (b - A @ center) / 1e-6
shifted_bounds = [(None if lo is None else (lo-center[k])/scale[k],
                   None if hi is None else (hi-center[k])/scale[k])
                  for k, (lo,hi) in enumerate(bounds)]
results = []
for name, rows in (("original_all_rows", list(range(len(A)))), ("reduced_rows", kept)):
    answer = linprog(objective/objective_scale, A_ub=scaled_A[rows], b_ub=scaled_b[rows],
                     bounds=shifted_bounds, method="highs",
                     options={"dual_feasibility_tolerance":1e-9,"primal_feasibility_tolerance":1e-9})
    assert answer.status in (0,2)
    results.append({"method":name,"status":int(answer.status),
                    "objective": None if answer.fun is None else float(answer.fun*objective_scale+c@center)})
own = ref._solve_centered_lp(c, A[kept], b[kept], bounds, center, scale)
results.append({"method":"current_revised_simplex_reduced","status":int(own.status),
                "objective": None if own.status else float(own.fun)})
assert len({r["status"] for r in results}) == 1
if own.status == 0:
    assert max(abs(r["objective"] - own.fun) for r in results) <= 1e-9
report = {"schema_version":1,"scope":"one retained failing development LP; no new task evaluation or score",
          "captured_lp_sha256":expected,"historical_solver_status":data["status"],
          "original_rows":len(A),"retained_rows":len(kept),"columns":len(c),
          "exact_rational_implications_proved":proven,"all_removed_rows_proved":proven == len(A)-len(kept),
          "independent_lp_results":results,"passed":True,
          "driver_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
with out.open("x") as handle:
    json.dump(report,handle,indent=2,sort_keys=True,allow_nan=False)
    handle.write("\n")
print(json.dumps({k:v for k,v in report.items() if k != "independent_lp_results"},sort_keys=True))
