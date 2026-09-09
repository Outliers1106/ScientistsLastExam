"""Shortcut probe: low-dimensional strategies that use the oracle's gene trees as given."""
import os as _os, tempfile as _tempfile
TASK = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "benchmarks", "Biology", "AnomalyZoneSpeciesTree"))
OUT = _tempfile.gettempdir()
import sys, json, itertools, collections, time, numpy as np
sys.path.insert(0, TASK + "/verification")
import evaluator as ev, msc
SETS, TABLE = msc.all_topologies()
FULL = msc.FULL_MASK
CELLS = [(rc, s) for rc in ev.RATE_CLASSES for s in ev.LOCUS_SITES]
WORLDS = [("development", i, s) for i, s in enumerate(ev.DEVELOPMENT_WORLDS)] + [("heldout", i, s) for i, s in enumerate(ev.HELDOUT_WORLDS)]

t0 = time.time()
worlds = {}
cache = {}
for split, index, spec in WORLDS:
    w = ev._world(spec); worlds[(split, index)] = w
    for rc, sites in CELLS:
        n = w["budget"] // ev.LOCUS_COST[sites]
        idx = [r["locus"] for r in w["catalogue"] if r["rate_class"] == rc and r["sites"] == sites][:n]
        splits, dists = [], []
        for i in idx:
            loc = ev._locus(w, i)
            _, adj = msc.parse_newick(loc["nj_tree"]); s = msc.splits_of(adj); splits.append(s)
            aln = msc.from_strings([loc["alignment"][t] for t in msc.TAXA]); dists.append(msc.jc_distances(aln))
        counts = np.zeros((70, 3))
        for s in splits:
            for q, t in enumerate(msc.quartet_signature(s)):
                if t >= 0: counts[q, t] += 1
        cache[(split, index, rc, sites)] = {"splits": splits, "mean_d": np.mean(dists, axis=0), "counts": counts}
print("cache built %.0fs" % (time.time() - t0), flush=True)

def compatible(a, b):
    return (a & b) == 0 or (a & b) == a or (a & b) == b or (a | b) == FULL
def greedy(splits_list):
    freq = collections.Counter(m for s in splits_list for m in s)
    chosen = {}
    for m, _ in sorted(freq.items(), key=lambda kv: (-kv[1], kv[0])):
        if all(compatible(m, c) for c in chosen):
            chosen[m] = 1.0
        if len(chosen) == 5: break
    return chosen
def mode(splits_list):
    c = collections.Counter(frozenset(s) for s in splits_list)
    best = max(c.items(), key=lambda kv: (kv[1], sorted(kv[0])))[0]
    return {m: 1.0 for m in best}
def avg_nj(mean_d):
    return msc.splits_of(msc.neighbour_joining(mean_d))
def quartet_cons(counts):
    score = counts[np.arange(70)[None, :], TABLE].sum(axis=1)
    return {m: 1.0 for m in SETS[int(np.argmax(score))]}
def chi35(counts):
    o = np.sort(counts, axis=1)[:, ::-1]
    x = (o[:, 1] - o[:, 2]) ** 2 / np.maximum(1.0, o[:, 1] + o[:, 2])
    return max(x[[q for q, quart in enumerate(msc.QUARTETS) if t in quart]].sum() for t in range(8))
def discordance(splits_list, topo):
    ref = set(topo)
    return float(np.mean([len(ref - set(s)) / 5.0 for s in splits_list]))
def quartet_lengths(topo, counts):
    tot = counts.sum(axis=1); out = {}
    for m in topo:
        est = []
        for q, t in msc.branch_length_quartets(topo, m):
            if tot[q] > 0:
                p = counts[q, t] / tot[q]; est.append(max(0.01, -np.log(max(1e-9, 1.5 * (1 - p)))))
        out[m] = float(np.median(est)) if est else 0.01
    return out

TOPOS = ["greedy", "mode", "avg_nj", "quartets"]
REFUSALS = [("never", None)] + [("discord", th) for th in (0.5, 0.6, 0.7)] + [("chi", th) for th in (80.0, 115.0, 160.0)]
LENGTHS = [0.05, 0.1, 0.2, 0.4, "quartet"]

topo_cache = {}
def topology(key, method):
    k = key + (method,)
    if k not in topo_cache:
        c = cache[key]
        topo_cache[k] = {"greedy": lambda: greedy(c["splits"]), "mode": lambda: mode(c["splits"]),
                         "avg_nj": lambda: avg_nj(c["mean_d"]), "quartets": lambda: quartet_cons(c["counts"])}[method]()
    return topo_cache[k]

results = []
for (rc, sites), method, (rname, th), length in itertools.product(CELLS, TOPOS, REFUSALS, LENGTHS):
    records = {"development": [], "heldout": []}
    for split, index, spec in WORLDS:
        key = (split, index, rc, sites); c = cache[key]; w = worlds[(split, index)]
        topo = topology(key, method)
        abstain = False
        if rname == "discord": abstain = discordance(c["splits"], topo) > th
        elif rname == "chi": abstain = chi35(c["counts"]) > th
        if abstain or len(topo) != 5:
            m = ev._metrics(w, None, 0.9, True)
        else:
            lengths = quartet_lengths(topo, c["counts"]) if length == "quartet" else {k: length for k in topo}
            m = ev._metrics(w, lengths, 0.9, False)
        m.update({"kind": w["kind"], "abstained": abstain, "valid": True, "confidence_calibration_score": 0.0, "loci_bought": 0, "budget_spent": 0})
        records[split].append(m)
    dev = ev._split_summary(records["development"]); held = ev._split_summary(records["heldout"])
    results.append({"class": rc, "sites": sites, "topology": method, "refusal": rname, "threshold": th, "length": length,
                    "dev": dev["normalized_mechanism"], "heldout": held["normalized_mechanism"], "topology_rate": dev["topology_rate"],
                    "fdr": dev["false_discovery_rate"], "refusal_rate": dev["correct_refusal_rate"]})
results.sort(key=lambda r: -r["dev"])
print("strategies", len(results))
for r in results[:25]:
    print("dev %.3f held %.3f topo %.2f fdr %.2f ref %.2f | %s %d %s %s %s %s" % (r["dev"], r["heldout"], r["topology_rate"], r["fdr"], r["refusal_rate"], r["class"], r["sites"], r["topology"], r["refusal"], r["threshold"], r["length"]))
print("best without quartet consensus:", max((r for r in results if r["topology"] != "quartets"), key=lambda r: r["dev"]))
print("best without quartet consensus or quartet lengths:", max((r for r in results if r["topology"] != "quartets" and r["length"] != "quartet"), key=lambda r: r["dev"]))
print("best fast-loci strategy:", max((r for r in results if r["class"] == "fast"), key=lambda r: r["dev"]))
print("best never-refuse strategy:", max((r for r in results if r["refusal"] == "never"), key=lambda r: r["dev"]))
json.dump(results, open(OUT + "/anomaly_zone_probe_results.json", "w"), indent=1)
print("total %.0fs" % (time.time() - t0))
