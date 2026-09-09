# AnomalyZoneSpeciesTree — reference results

Every number here is produced by running code in this directory or the construction scripts in
`.research/anomaly_zone_species_tree/` at the repository root. Nothing is copied from a table.

## Reproducing

```
python3 frontier_eval/run_eval.py --candidate verification/reference_quartet_consensus.py \
    --metrics-out /tmp/metrics.json
```

## Reference — `verification/reference_quartet_consensus.py`

Truth-blind: it reads only the public problem and the budgeted sequencing campaign.

| metric | development | held out |
|---|---|---|
| mechanism score (normalized) | **0.877** | 0.929 |
| topology rate | 1.00 | 1.00 |
| branch-length score | 0.75 | 0.86 |
| false discovery rate | 0.00 | 0.00 |
| correct refusal rate | 1.00 | 1.00 |
| discovery coverage | 1.00 | 1.00 |
| confidence calibration | 0.96 | — |
| loci bought | 240 of 2700 | 240 of 2700 |
| budget spent | 480 of 480 | 480 of 480 |

Its design: buy medium-rate loci of 800 sites in catalogue order until the budget is spent; discard
the sequencing centre's gene trees and re-estimate each one by neighbour joining on the
Jukes-Cantor distance corrected for the published gamma shape; count the three topologies of each
of the 70 quartets across the gene trees; decline when, for some species, the sum over the 35
quartets containing it of the squared difference between the two minority counts divided by their
sum exceeds 115; otherwise take the one of the 10395 unrooted eight-taxon trees that agrees with
the most quartet observations, and read each internal branch length from the median majority
frequency of the quartets whose path is exactly that branch, inverting `p = 1 - (2/3) exp(-t)`.

**The reference is deliberately not at the ceiling.** Its branch-length score is 0.75 on the
development split. Three things are left on the table: every quartet is counted with the same
weight whether the gene tree resolved it well or barely; the inversion ignores gene-tree estimation
error, which flattens the frequencies towards a third and biases every length downward; and the
refusal is a fixed threshold on one statistic rather than a test of the fitted tree against the
whole quartet spectrum. One evaluation of the reference takes about seven seconds on a laptop.

## Model draws — not performed

No frontier model has seen this task. It was built on a machine without a model endpoint, so the
admission bar (no first proposal reaches the reference) is not measured. `lineage.calibration_runs`
on the task card is empty and `calibration_evidence_status` says so. The first draw should follow
the exemplar's protocol: three seeds, three proposals each, `greedy_rewrite`, budget 3, against
the reference's 0.877.

## Baseline — `solution.py`

Buys the longest fast loci until the budget is spent (96 loci of 2000 sites), concatenates them,
computes the plain Jukes-Cantor distance on the whole alignment, builds one neighbour-joining tree
and reports its internal lengths in substitutions per site as if they were coalescent units.
Never declines, confidence 0.9.

| metric | development | held out |
|---|---|---|
| combined score | **0.0000** | 0.0000 |
| raw mechanism | 0.043 | — |
| topology rate | 0.12 | 0.25 |
| branch-length score | 0.01 | 0.01 |
| false discovery rate | 0.92 | 0.83 |
| correct refusal rate | 0.00 | 0.00 |

Confidently wrong rather than empty: eleven of twelve development worlds are false discoveries.
Concatenation is inconsistent in the anomaly zone, the uncorrected distance on the fastest loci
joins the two fast-evolving species, and the hybrid worlds get a tree like everything else. The
one topology it does recover is scored on branch lengths a hundredfold too small.

## Difficulty ladder

One choice of the reference changed at a time. Held out is the normalized held-out mechanism.

| strategy | score | topology | branch lengths | false discovery | refusal | coverage | loci | held out |
|---|---|---|---|---|---|---|---|---|
| gamma-corrected gene trees, 240 medium 800-site loci, quartet consensus, quartet lengths, imbalance refusal (reference) | **0.877** | 1.00 | 0.75 | 0.00 | 1.00 | 1.00 | 240 | 0.929 |
| same, the sequencing centre's uncorrected gene trees | 0.316 | 0.38 | 0.26 | 0.42 | 1.00 | 0.75 | 240 | 0.467 |
| same, slow loci | 0.796 | 0.88 | 0.72 | 0.08 | 1.00 | 1.00 | 240 | 0.683 |
| same, fast loci | 0.542 | 0.75 | 0.58 | 0.25 | 0.75 | 1.00 | 240 | 0.904 |
| same, 300-site loci (all 300 in the catalogue) | 0.553 | 0.75 | 0.61 | 0.25 | 0.75 | 1.00 | 300 | 0.874 |
| same, 2000-site loci (96 of them) | 0.290 | 0.62 | 0.46 | 0.42 | 0.50 | 1.00 | 96 | 0.617 |
| same, never declining | 0.377 | 1.00 | 0.75 | 0.33 | 0.00 | 1.00 | 240 | 0.429 |
| same, refusal threshold 80 | 0.877 | 1.00 | 0.75 | 0.00 | 1.00 | 1.00 | 240 | 0.696 |
| same, refusal threshold 160 | 0.877 | 1.00 | 0.75 | 0.00 | 1.00 | 1.00 | 240 | 0.929 |
| same, constant branch lengths of 0.1 | 0.744 | 1.00 | 0.49 | 0.00 | 1.00 | 1.00 | 240 | 0.705 |
| same, greedy consensus of the gene trees instead of quartets | 0.452 | 0.50 | 0.40 | 0.33 | 1.00 | 1.00 | 240 | 0.694 |
| same, half the budget | 0.528 | 0.88 | 0.68 | 0.25 | 0.50 | 1.00 | 120 | 0.913 |
| same, a quarter of the budget | 0.000 | 0.50 | 0.38 | 0.67 | 0.00 | 1.00 | 60 | 0.390 |
| uncorrected gene trees and fast loci | 0.100 | 0.25 | 0.20 | 0.58 | 0.75 | 0.38 | 240 | 0.223 |
| uncorrected gene trees, never declining | 0.000 | 0.38 | 0.26 | 0.75 | 0.00 | 1.00 | 240 | 0.404 |
| concatenation of the longest fast loci, plain Jukes-Cantor, never declining (baseline) | 0.000 | 0.12 | 0.01 | 0.92 | 0.00 | 1.00 | 96 | 0.000 |
| declining everything | 0.000 | — | — | 0.67 | 1.00 | 0.00 | 0 | 0.000 |

The ladder is not a difficulty measurement; that would be a frontier draw. What it shows is where
the score lives: re-estimating the gene trees is worth 0.56, declining 0.50, quartets over greedy
consensus 0.43, the locus length 0.32 or 0.59, the rate class 0.08 or 0.34, half the budget 0.35
and the branch lengths 0.13. The refusal threshold has slack on the development split, where the
statistic is at most 69.5 on a tree world and at least 171.6 on a reticulate one; held out a tree
world reaches 82.6, which is why 80 loses there and 160 does not.

Capping the reference at a fixed number of loci, in catalogue order, regardless of budget:

| loci | development | held out |
|---|---|---|
| 25 | 0.000 | 0.000 |
| 50 | 0.000 | 0.130 |
| 100 | 0.527 | 0.682 |
| 240 (reference) | 0.877 | 0.929 |

## Shortcut probe

The question this repository learned to ask after a submitted task turned out to be solvable by a
two-parameter grid search: **how far does a low-dimensional strategy get without the science?**

1260 strategies that take the sequencing centre's gene trees as given were evaluated: one of the
nine catalogue cells (rate class by length, as many loci as the budget and the catalogue allow),
topology by greedy consensus, the most frequent gene tree, neighbour joining on the mean
Jukes-Cantor distance, or quartet consensus; refusal never, on a mean gene-tree discordance
threshold of 0.5, 0.6 or 0.7, or on the minority-imbalance statistic at 80, 115 or 160; branch
lengths constant at 0.05, 0.1, 0.2 or 0.4, or from the quartet frequencies.

| family | strategies | best development | its held out |
|---|---|---|---|
| all | 1260 | **0.838** | 0.604 |
| quartet consensus on the free gene trees | 315 | 0.785 | 0.679 |
| neighbour joining on the mean distance | 315 | 0.838 | 0.604 |
| greedy consensus | 315 | 0.348 | 0.458 |
| most frequent gene tree | 315 | 0.348 | 0.235 |
| no quartet consensus and no quartet lengths | 756 | 0.744 | 0.509 |
| never declining | 180 | 0.338 | 0.104 |
| declining on gene-tree discordance | 540 | 0.234 | 0.000 |
| fast loci | 420 | 0.331 | 0.414 |
| medium loci | 420 | 0.438 | 0.678 |
| 2000-site loci | 420 | 0.555 | 0.657 |
| reference | — | 0.877 | 0.929 |
| ceiling | — | 1.000 | 1.000 |

The top strategy is neighbour joining on the mean Jukes-Cantor distance over the 300 slow 300-site
loci, declining on the imbalance statistic at 115, with quartet-frequency lengths. Averaging
distances over loci is a consistent species-tree method under the coalescent, so this is not a
science-free shortcut; it stays below the reference on both splits because slow short loci resolve
the short branches badly, which flattens the quartet frequencies its branch lengths are read from.
Every strategy that uses the free gene trees on medium or fast loci, where the uncorrected distance
matters, stays under 0.44. The probe's own weakness is that it never re-estimates a gene tree.

## Four construction errors, all found by the checkpoints before any model saw the task

- **There was no long-branch attraction.** The first two-fast-species worlds put the fast species
  on ordinary branches, and the uncorrected distance moved zero to half a per cent of gene trees.
  The analytic four-point condition showed why: attraction needs the fast species' sisters to be
  short and the branch between them shorter. The worlds are now two cherries joined by a branch of
  0.1 to 0.2 coalescent units, with the fast species one from each cherry at five to eight times
  the rate, and the attraction is ten to thirty per cent of gene trees on fast loci.
- **The budget could not buy the answer.** At 80 units no quartet strategy recovered any
  anomaly-zone tree and the refusal statistic did not separate hybrid from tree worlds. A
  sorted-difference statistic also turned out to have a null mean independent of the sample
  size, so it could never be thresholded; it was replaced by the chi-square form. The budget is now
  480 with a catalogue of 2700, and the branches were lengthened until the reference reads them.
- **One anomaly world was not anomalous.** At short branches of 0.12 to 0.2 coalescent units one
  held-out species tree was the most frequent gene tree by simulation. The range is now 0.1 to
  0.16, and every world is checked with one hundred thousand gene trees: the species tree ranks
  fourth, seventh, twenty-seventh, ninth, third and third, each at least three standard errors
  behind the top topology.
- **The reference tied the probe.** At a budget of 400 the reference lost one development world
  by 0.2 quartet votes out of 243, and the probe's best strategy tied it there at 0.779 against
  0.774 while trailing held out at 0.682 against 0.917. Neither reweighting the quartets, mixing
  rate classes nor bootstrap support moved that world off the edge; each flipped a different
  world. The budget was raised to 480, at which the reference recovers every tree world and leads
  the whole probe on both splits.

## Robustness

- Twenty malformed candidate shapes — raising, `None`, an empty mapping, a string, a bad verdict,
  a tree verdict without `newick`, an integer `newick`, garbage Newick, a missing taxon, a
  duplicated taxon, no branch lengths, a negative length, a NaN length, an unresolved tree, a NaN
  confidence, a string confidence, overspending the budget, a locus of 9999, a float locus and a
  Newick string of five thousand characters — all score zero with `valid = 0`, and none raises out
  of the evaluator.
- Two consecutive evaluations of the reference are key-identical, and so are their JSON payloads.
- Declining every world scores exactly 0.0 by construction of the normalisation, in the verdict
  form and in the `abstain: true` form; naming one fixed tree everywhere scores 0.0 with a false
  discovery rate of 1.0.
- Every source file compiles under Python 3.8 syntax, which is what the evaluation host runs.
