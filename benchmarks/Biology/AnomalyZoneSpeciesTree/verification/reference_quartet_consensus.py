"""Truth-blind reference for AnomalyZoneSpeciesTree: re-estimate every gene tree under the
published model, assemble the species tree from quartet majorities, read branch lengths off the
quartet frequencies, and decline when the minority quartets are too unequal for any tree.

Reads only the public problem and the budgeted sequencing campaign. The ideas it is built on:

    loci        buy medium-rate loci of 800 sites. Slow loci leave a branch of a tenth of a
                coalescent unit unresolved in most gene trees; fast loci resolve it and saturate
                the long distances. The medium class is the compromise, and 800 sites at that
                rate buys two hundred and forty loci, which is what the quartet frequencies need.
    gene trees  the oracle's neighbour-joining tree is discarded. It uses the plain Jukes-Cantor
                distance, and the sites evolve with gamma-distributed rates of published shape,
                so the distances are recomputed with the gamma correction and neighbour joining
                is run again. Without this the two fast-evolving species attract each other,
                and - worse - the attraction is a systematic imbalance between minority quartet
                topologies, which is exactly the signature of reticulation.
    topology    count the three topologies of every quartet across the gene trees and pick the
                one of the 10395 unrooted eight-taxon trees that agrees with the most quartet
                observations. No quartet is anomalous under the coalescent, so this is
                consistent where the most frequent gene tree, the greedy consensus and
                concatenation are not.
    refusal     under any species tree the two minority topologies of a quartet have equal
                probability. For each species, sum (n1 - n2)^2 / (n1 + n2) over the 35 quartets
                that contain it; a hybrid species lifts this far above what sampling alone does.
    lengths     an internal branch of t coalescent units gives its quartets a majority topology
                with probability 1 - (2/3) exp(-t); invert that on the median frequency over the
                quartets whose path is exactly that branch.

Deliberately not at the ceiling. Three things are left: the gene trees are distance trees, and
every quartet is counted with the same weight whether the gene tree resolved it well or barely;
the branch-length inversion ignores gene-tree estimation error, which flattens the frequencies
towards a third and biases every length downward; and the refusal is a fixed threshold on one
statistic rather than a test of the fitted tree against the whole quartet spectrum.
"""
from __future__ import annotations

import itertools
import math

import numpy as np

CHOSEN_CLASS = "medium"
CHOSEN_SITES = 800
REFUSAL_STATISTIC = 115.0
MIN_LENGTH = 0.01
PENDANT_LENGTH = 1.0
SATURATED = 10.0


# ----------------------------------------------------------------------------------------------
# distances and neighbour joining
# ----------------------------------------------------------------------------------------------

def _encode(alignment, taxa):
    table = {"A": 0, "C": 1, "G": 2, "T": 3}
    return np.array([[table.get(ch, 0) for ch in alignment[t]] for t in taxa], dtype=np.int8)


def _gamma_distances(seqs, shape):
    n = seqs.shape[0]
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            p = float(np.mean(seqs[i] != seqs[j]))
            if p >= 0.75:
                value = SATURATED
            else:
                value = 0.75 * shape * ((1.0 - 4.0 * p / 3.0) ** (-1.0 / shape) - 1.0)
            d[i, j] = d[j, i] = value
    return d


def _neighbour_joining(distance):
    d = np.array(distance, dtype=float)
    n = d.shape[0]
    active = list(range(n))
    adjacency = {i: [] for i in range(n)}
    nxt = n
    while len(active) > 2:
        m = len(active)
        sub = d[np.ix_(active, active)]
        totals = sub.sum(axis=1)
        q = (m - 2) * sub - totals[:, None] - totals[None, :]
        np.fill_diagonal(q, np.inf)
        i, j = divmod(int(np.argmin(q)), m)
        if i > j:
            i, j = j, i
        a, b = active[i], active[j]
        dij = sub[i, j]
        la = max(0.5 * dij + (totals[i] - totals[j]) / (2.0 * (m - 2)), 0.0)
        lb = max(dij - la, 0.0)
        u = nxt
        nxt += 1
        adjacency[u] = [(a, la), (b, lb)]
        adjacency[a].append((u, la))
        adjacency[b].append((u, lb))
        d = np.pad(d, ((0, 1), (0, 1)))
        for k in active:
            if k not in (a, b):
                d[u, k] = d[k, u] = 0.5 * (d[a, k] + d[b, k] - dij)
        active = [k for k in active if k not in (a, b)] + [u]
    a, b = active
    adjacency[a].append((b, max(float(d[a, b]), 0.0)))
    adjacency[b].append((a, max(float(d[a, b]), 0.0)))
    return adjacency


# ----------------------------------------------------------------------------------------------
# splits and quartets
# ----------------------------------------------------------------------------------------------

N = 8
FULL = (1 << N) - 1
QUARTETS = tuple(itertools.combinations(range(N), 4))


def _side(adjacency, start, blocked):
    mask = 0
    stack = [(start, blocked)]
    while stack:
        node, parent = stack.pop()
        if node < N:
            mask |= 1 << node
        for nbr, _ in adjacency[node]:
            if nbr != parent:
                stack.append((nbr, node))
    return mask


def _canonical(mask):
    return (FULL ^ mask) if mask & 1 else mask


def _splits(adjacency):
    found = {}
    for node in adjacency:
        for nbr, length in adjacency[node]:
            if node < nbr:
                mask = _canonical(_side(adjacency, nbr, node))
                if 2 <= bin(mask).count("1") <= N - 2:
                    found[mask] = length
    return found


def _quartet_topology(splits, quartet):
    a, b, c, d = quartet
    for mask in splits:
        inside = [(mask >> t) & 1 for t in (a, b, c, d)]
        if sum(inside) != 2:
            continue
        if inside[0] == inside[1]:
            return 0
        if inside[0] == inside[2]:
            return 1
        return 2
    return -1


def _all_trees():
    """Every unrooted binary tree on eight taxa by stepwise addition: edges, splits, quartets."""
    trees = [[(0, 8), (1, 8), (2, 8)]]
    nxt = 9
    for taxon in range(3, N):
        grown = []
        for edges in trees:
            for k, (u, v) in enumerate(edges):
                grown.append(edges[:k] + edges[k + 1:] + [(u, nxt), (nxt, v), (nxt, taxon)])
        trees = grown
        nxt += 1
    table = np.zeros((len(trees), len(QUARTETS)), dtype=np.int8)
    split_sets = []
    for row, edges in enumerate(trees):
        adjacency = {}
        for u, v in edges:
            adjacency.setdefault(u, []).append((v, 1.0))
            adjacency.setdefault(v, []).append((u, 1.0))
        splits = _splits(adjacency)
        split_sets.append(splits)
        table[row] = [_quartet_topology(splits, q) for q in QUARTETS]
    return trees, split_sets, table


# ----------------------------------------------------------------------------------------------
# branch lengths from quartet frequencies
# ----------------------------------------------------------------------------------------------

def _parts(splits, side):
    """Maximal clades strictly inside `side` (the subtrees hanging off the branch's endpoint)."""
    candidates = set()
    for mask in list(splits) + [1 << t for t in range(N)]:
        for oriented in (mask, FULL ^ mask):
            if oriented and oriented != side and not (oriented & ~side & FULL):
                candidates.add(oriented)
    parts, covered = [], 0
    for mask in sorted(candidates, key=lambda m: -bin(m).count("1")):
        if not mask & covered:
            parts.append(mask)
            covered |= mask
    return parts


def _members(mask):
    return [t for t in range(N) if (mask >> t) & 1]


def _branch_quartets(splits, mask):
    result = []
    for pa, pb in itertools.combinations(_parts(splits, mask), 2):
        for pc, pd in itertools.combinations(_parts(splits, FULL ^ mask), 2):
            for a in _members(pa):
                for b in _members(pb):
                    for c in _members(pc):
                        for d in _members(pd):
                            quartet = tuple(sorted((a, b, c, d)))
                            result.append((QUARTETS.index(quartet), _quartet_topology(splits, quartet)))
    return result


def _render(edges, lengths, labels):
    adjacency = {}
    for (u, v), length in zip(edges, lengths):
        adjacency.setdefault(u, []).append((v, length))
        adjacency.setdefault(v, []).append((u, length))
    root = max(adjacency)

    def render(node, parent):
        kids = [(nbr, length) for nbr, length in adjacency[node] if nbr != parent]
        if not kids:
            return labels[node]
        return "(" + ",".join("%s:%.6f" % (render(nbr, node), length) for nbr, length in kids) + ")"

    return render(root, None) + ";"


# ----------------------------------------------------------------------------------------------
# the candidate
# ----------------------------------------------------------------------------------------------

_CACHE = {}


def infer_species_tree(problem, sequence):
    taxa = list(problem["taxa"])
    shape = float(problem["gamma_shape"])
    budget = int(problem["locus_budget"])
    catalogue = problem["catalogue"]

    chosen = [row for row in catalogue
              if row["rate_class"] == CHOSEN_CLASS and int(row["sites"]) == CHOSEN_SITES]
    spent = 0
    counts = np.zeros((len(QUARTETS), 3))
    loci = 0
    for row in chosen:
        cost = int(row["cost"])
        if spent + cost > budget:
            break
        try:
            locus = sequence(int(row["locus"]))
        except Exception:
            break
        spent += cost
        seqs = _encode(locus["alignment"], taxa)
        splits = _splits(_neighbour_joining(_gamma_distances(seqs, shape)))
        for q, quartet in enumerate(QUARTETS):
            topology = _quartet_topology(splits, quartet)
            if topology >= 0:
                counts[q, topology] += 1
        loci += 1
    if loci == 0:
        return {"verdict": "reticulate", "confidence": 0.1}

    # Refusal: minority imbalance summed over the quartets containing each species.
    ordered = np.sort(counts, axis=1)[:, ::-1]
    imbalance = (ordered[:, 1] - ordered[:, 2]) ** 2 / np.maximum(1.0, ordered[:, 1] + ordered[:, 2])
    worst = 0.0
    for taxon in range(N):
        idx = [q for q, quartet in enumerate(QUARTETS) if taxon in quartet]
        worst = max(worst, float(imbalance[idx].sum()))
    if worst > REFUSAL_STATISTIC:
        return {"verdict": "reticulate", "confidence": 0.7}

    if "trees" not in _CACHE:
        _CACHE["trees"] = _all_trees()
    edges_list, split_sets, table = _CACHE["trees"]
    score = counts[np.arange(len(QUARTETS))[None, :], table].sum(axis=1)
    best = int(np.argmax(score))
    edges, splits = edges_list[best], split_sets[best]

    lengths = {}
    totals = counts.sum(axis=1)
    for mask in splits:
        estimates = []
        for q, topology in _branch_quartets(splits, mask):
            if totals[q] <= 0:
                continue
            p = counts[q, topology] / totals[q]
            estimates.append(max(MIN_LENGTH, -math.log(max(1e-9, 1.5 * (1.0 - p)))))
        lengths[mask] = float(np.median(estimates)) if estimates else MIN_LENGTH

    adjacency = {}
    for u, v in edges:
        adjacency.setdefault(u, []).append(v)
        adjacency.setdefault(v, []).append(u)
    rendered = []
    for u, v in edges:
        adj = {n: [(m, 1.0) for m in adjacency[n]] for n in adjacency}
        mask = _canonical(_side(adj, v, u))
        rendered.append(lengths.get(mask, PENDANT_LENGTH))
    return {"verdict": "tree", "newick": _render(edges, rendered, taxa), "confidence": 0.75}
