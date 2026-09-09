import os as _os, tempfile as _tempfile
TASK = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "benchmarks", "Biology", "AnomalyZoneSpeciesTree"))
OUT = _tempfile.gettempdir()
import sys, collections, numpy as np
sys.path.insert(0, TASK + "/verification")
import evaluator as ev, msc
N = 100000
for short in [(0.12, 0.2), (0.1, 0.16)]:
    ev.SHORT_BRANCH = short
    print("SHORT_BRANCH", short)
    for spec in [s for s in ev.DEVELOPMENT_WORLDS + ev.HELDOUT_WORLDS if s["kind"] == "anomaly"]:
        w = ev._world(spec); truth = frozenset(w["truth_splits"])
        rng = np.random.default_rng(5); ctr = collections.Counter()
        for i in range(N):
            g, _ = msc.simulate_gene_tree(w["trees"][0], rng); ctr[frozenset(g.unrooted_splits())] += 1
        freqs = sorted(ctr.values(), reverse=True)
        top = freqs[0] / N; mine = ctr[truth] / N
        se = np.sqrt(top * (1 - top) / N)
        print("  %d species %.4f top %.4f rank %d  gap/se %.1f" % (spec["seed"], mine, top, freqs.index(ctr[truth]) + 1, (top - mine) / se))
