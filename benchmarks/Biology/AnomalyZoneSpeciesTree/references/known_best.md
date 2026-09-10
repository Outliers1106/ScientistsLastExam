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
| mechanism score (normalized) | **0.873** | 0.929 |
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
whole quartet spectrum. One evaluation of the reference takes about ten seconds on the
construction laptop and under thirty on the evaluation host, which is what `eval_time_seconds`
says.

## Model draws — not performed

No frontier model has seen this task. It was built on a machine without a model endpoint, so the
admission bar (no first proposal reaches the reference) is not measured. `lineage.calibration_runs`
on the task card is empty and `calibration_evidence_status` says so. The first draw should follow
the exemplar's protocol: three seeds, three proposals each, `greedy_rewrite`, budget 3, against
the reference's 0.873.

## Baseline — `solution.py`

Buys the longest fast loci until the budget is spent (96 loci of 2000 sites), concatenates them,
computes the plain Jukes-Cantor distance on the whole alignment, builds one neighbour-joining tree
and reports its internal lengths in substitutions per site as if they were coalescent units.
Never declines, confidence 0.9.

| metric | development | held out |
|---|---|---|
| combined score | **0.0000** | 0.0000 |
| raw mechanism | 0.000 | — |
| topology rate | 0.00 | 0.25 |
| branch-length score | 0.00 | 0.01 |
| false discovery rate | 1.00 | 0.83 |
| correct refusal rate | 0.00 | 0.00 |

Confidently wrong rather than empty: every development world is a false discovery. Concatenation
is inconsistent in the anomaly zone, the uncorrected distance on the fastest loci joins the two
fast-evolving species, and the hybrid worlds get a tree like everything else. The one held-out
topology it does recover is scored on branch lengths a hundredfold too small.

## Difficulty ladder

One choice of the reference changed at a time. Held out is the normalized held-out mechanism.

| strategy | score | topology | branch lengths | false discovery | refusal | coverage | loci | held out |
|---|---|---|---|---|---|---|---|---|
| gamma-corrected distance gene trees from 240 medium loci, quartet consensus, quartet-frequency lengths, minority-imbalance refusal | **0.873** | 1.00 | 0.75 | 0.00 | 1.00 | 1.00 | 240 | 0.929 |
| same, but the sequencing centre's uncorrected gene trees | 0.316 | 0.38 | 0.26 | 0.42 | 1.00 | 0.38 | 240 | 0.467 |
| same, slow loci instead of medium | 0.891 | 1.00 | 0.78 | 0.00 | 1.00 | 1.00 | 240 | 0.662 |
| same, fast loci instead of medium | 0.639 | 0.88 | 0.65 | 0.17 | 0.75 | 1.00 | 240 | 0.854 |
| same, 300-site loci (all 300 of them) instead of 800 | 0.405 | 0.62 | 0.44 | 0.33 | 0.75 | 1.00 | 300 | 0.839 |
| same, 2000-site loci (96 of them) instead of 800 | 0.293 | 0.62 | 0.46 | 0.42 | 0.50 | 1.00 | 96 | 0.623 |
| same, never declining | 0.373 | 1.00 | 0.75 | 0.33 | 0.00 | 1.00 | 240 | 0.429 |
| same, refusal threshold 80 | 0.758 | 0.88 | 0.64 | 0.08 | 1.00 | 0.88 | 240 | 0.704 |
| same, refusal threshold 160 | 0.873 | 1.00 | 0.75 | 0.00 | 1.00 | 1.00 | 240 | 0.929 |
| same, constant branch lengths of 0.1 | 0.743 | 1.00 | 0.49 | 0.00 | 1.00 | 1.00 | 240 | 0.701 |
| same, greedy consensus of the gene trees instead of quartets | 0.448 | 0.50 | 0.40 | 0.33 | 1.00 | 1.00 | 240 | 0.694 |
| same, half the budget | 0.528 | 0.88 | 0.68 | 0.25 | 0.50 | 1.00 | 120 | 0.916 |
| same, a quarter of the budget | 0.000 | 0.50 | 0.36 | 0.67 | 0.00 | 1.00 | 60 | 0.626 |
| uncorrected gene trees and fast loci | 0.100 | 0.25 | 0.20 | 0.58 | 0.75 | 0.38 | 240 | 0.223 |
| uncorrected gene trees, never declining | 0.000 | 0.38 | 0.26 | 0.75 | 0.00 | 1.00 | 240 | 0.000 |
| concatenation of the longest fast loci, plain Jukes-Cantor, never declining (baseline) | 0.000 | 0.00 | 0.00 | 1.00 | 0.00 | 1.00 | 96 | 0.000 |
| declining everything | 0.000 | — | — | 0.67 | 1.00 | 0.00 | 0 | 0.000 |

The ladder is not a difficulty measurement; that would be a frontier draw. What it shows is where
the score lives: re-estimating the gene trees is worth 0.56, declining 0.50, quartets over greedy
consensus 0.42, the locus length 0.47 or 0.58, fast loci 0.23, half the budget 0.35 and the
branch lengths 0.13. Slow loci gain 0.02 on the development split and lose 0.27 held out, so the
medium class is the robust choice rather than the best one on the visible worlds. The refusal
threshold has slack: on the reference's loci the statistic is at most 91.7 on a tree world and at
least 171.6 on a reticulate one, which is why 160 changes nothing and 80 declines a tree world and loses 0.12.

Why the free gene trees cost 0.56 is worth stating precisely, because the first version of this
document got it wrong. Per locus the free trees are not uniformly worse: on the anomaly-zone
worlds they recover the true gene tree's internal splits as often as the corrected trees do (0.84
against 0.83 on medium loci). On the long-branch worlds they are worse per locus (0.66 against
0.79) and, what matters more, their errors have a direction: they join the two fast species, so
over 240 loci the errors add up instead of cancelling. That is what turns a per-locus estimation
error into a wrong species tree and a spurious reticulation signal.

Capping the reference at a fixed number of loci, in catalogue order, regardless of budget:

| loci | development | held out |
|---|---|---|
| 25 | 0.000 | 0.000 |
| 50 | 0.000 | 0.354 |
| 100 | 0.515 | 0.904 |
| 240 (reference) | 0.873 | 0.929 |

## Long-branch attraction, measured

The quantity is the share of gene trees that join the two fast species as sisters, in the true
gene trees of the loci (that share is incomplete lineage sorting, not attraction), in the
sequencing centre's plain Jukes-Cantor trees and in gamma-corrected trees. The attraction is the
excess of an estimated tree over the truth. 480 loci of 300 sites or 240 of 800, per world and
rate class, from `lba.py`.

| world | class | sites | true | free | excess | corrected | excess |
|---|---|---|---|---|---|---|---|
| 71300105 | slow | 300 | 13.7 | 26.0 | +12.3 | 19.3 | +5.7 |
| 71300105 | slow | 800 | 8.8 | 29.6 | +20.8 | 17.9 | +9.2 |
| 71300105 | medium | 800 | 14.6 | 41.2 | +26.7 | 17.5 | +2.9 |
| 71300105 | fast | 800 | 10.4 | 55.0 | +44.6 | 19.6 | +9.2 |
| 71300106 | slow | 300 | 18.0 | 40.7 | +22.7 | 24.0 | +6.0 |
| 71300106 | slow | 800 | 17.1 | 47.1 | +30.0 | 23.8 | +6.7 |
| 71300106 | medium | 800 | 19.2 | 70.0 | +50.8 | 22.5 | +3.3 |
| 71300106 | fast | 800 | 16.2 | 76.7 | +60.4 | 23.8 | +7.5 |
| 71300107 | slow | 300 | 13.7 | 22.3 | +8.7 | 19.0 | +5.3 |
| 71300107 | slow | 800 | 10.8 | 17.1 | +6.2 | 10.8 | +0.0 |
| 71300107 | medium | 800 | 13.3 | 32.5 | +19.2 | 17.9 | +4.6 |
| 71300107 | fast | 800 | 14.2 | 33.3 | +19.2 | 17.9 | +3.8 |
| 71300108 | slow | 300 | 11.0 | 17.7 | +6.7 | 14.0 | +3.0 |
| 71300108 | slow | 800 | 20.8 | 26.7 | +5.8 | 19.6 | -1.3 |
| 71300108 | medium | 800 | 13.8 | 30.0 | +16.2 | 19.6 | +5.8 |
| 71300108 | fast | 800 | 16.2 | 40.4 | +24.2 | 17.9 | +1.7 |
| 82410203 | slow | 300 | 16.3 | 35.7 | +19.3 | 21.7 | +5.3 |
| 82410203 | slow | 800 | 17.9 | 39.2 | +21.2 | 17.5 | -0.4 |
| 82410203 | medium | 800 | 12.5 | 59.2 | +46.7 | 20.8 | +8.3 |
| 82410203 | fast | 800 | 15.4 | 68.8 | +53.3 | 20.4 | +5.0 |
| 82410204 | slow | 300 | 14.0 | 43.0 | +29.0 | 22.7 | +8.7 |
| 82410204 | slow | 800 | 16.7 | 51.7 | +35.0 | 22.1 | +5.4 |
| 82410204 | medium | 800 | 10.0 | 71.7 | +61.7 | 23.8 | +13.8 |
| 82410204 | fast | 800 | 12.1 | 84.2 | +72.1 | 27.1 | +15.0 |

Percentages of gene trees; the standard error of each share is 2 to 3 points. The excess of the
free trees is positive in every world and every class, and grows with the rate. The gamma
correction removes most of it and not all: the corrected distance is unbiased but noisy on the
long branches, and neighbour joining on noisy distances keeps a small attraction of its own. That
residual is one of the things a better gene-tree estimator would remove.

The same bias, seen from the refusal statistic and from the mean-distance tree:

| class | sites | imbalance on free trees, max over tree worlds | min over reticulate worlds | on corrected trees, max tree | min reticulate |
|---|---|---|---|---|---|
| slow | 300 | 268.6 | 119.4 | 111.5 | 102.3 |
| slow | 800 | 204.8 | 125.7 | 97.2 | 124.8 |
| medium | 800 | 411.9 | 175.8 | 91.7 | 171.6 |
| fast | 800 | 609.3 | 101.3 | 83.0 | 85.9 |

On the free trees no threshold separates a tree world from a hybrid world in any class. On the
corrected trees the reference's class and length separate them with the threshold at 115. On
slow 300-site loci even the corrected trees do not, which is one more reason the reference does
not buy them. Neighbour joining on the mean uncorrected distance recovers the species tree in one
of the six long-branch worlds on slow 300-site loci and in none on the other three cells; on the
mean corrected distance it recovers all six in every cell.

## Anomaly-zone statistics

For every tree world, from `anomaly_stats.py`: the probability that a gene tree has the species
tree's unrooted topology, the probability of the most frequent topology, and the difference in
units of its multinomial standard error. Four hundred thousand gene trees for the anomaly-zone
worlds, one hundred thousand for the long-branch worlds, generator seed 5.

| world | kind | p(species tree) | p(most frequent) | difference / se | rank of the species tree |
|---|---|---|---|---|---|
| 71300101 | anomaly | 0.0086 | 0.0107 | 9.5 | 4 |
| 71300102 | anomaly | 0.0027 | 0.0034 | 5.4 | 9 |
| 71300103 | anomaly | 0.0018 | 0.0023 | 5.2 | 22 |
| 71300104 | anomaly | 0.0035 | 0.0044 | 6.7 | 9 |
| 82410201 | anomaly | 0.0085 | 0.0105 | 9.3 | 4 |
| 82410202 | anomaly | 0.0096 | 0.0113 | 7.3 | 3 |
| 71300105 | long branch | 0.0690 | 0.0690 | 0.0 | 1 |
| 71300106 | long branch | 0.0480 | 0.0480 | 0.0 | 1 |
| 71300107 | long branch | 0.0929 | 0.0929 | 0.0 | 1 |
| 71300108 | long branch | 0.0350 | 0.0350 | 0.0 | 1 |
| 82410203 | long branch | 0.0486 | 0.0486 | 0.0 | 1 |
| 82410204 | long branch | 0.0571 | 0.0571 | 0.0 | 1 |

The six anomaly-zone worlds are anomalous by at least five standard errors. The rank is reported
for orientation only: with thousands of distinct topologies at probabilities of a few thousandths,
it moves between runs, and the difference and its standard error are the measurement. The six
long-branch worlds are built with at most one short internal branch and are outside the zone:
there the species tree is the most frequent gene tree, majority vote is consistent, and what
breaks it is the bias in the free trees, not the coalescent.

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
| all | 1260 | **0.512** | 0.197 |
| quartet consensus on the free gene trees | 315 | 0.437 | 0.221 |
| neighbour joining on the mean distance | 315 | 0.512 | 0.197 |
| greedy consensus | 315 | 0.107 | 0.000 |
| most frequent gene tree | 315 | 0.221 | 0.000 |
| no quartet consensus and no quartet lengths | 756 | 0.500 | 0.190 |
| never declining | 180 | 0.059 | 0.000 |
| declining on gene-tree discordance | 540 | 0.221 | 0.000 |
| slow loci | 420 | 0.512 | 0.197 |
| medium loci | 420 | 0.316 | 0.467 |
| fast loci | 420 | 0.331 | 0.414 |
| 2000-site loci | 420 | 0.316 | 0.000 |
| reference | — | 0.873 | 0.929 |
| ceiling | — | 1.000 | 1.000 |

The top strategy is neighbour joining on the mean Jukes-Cantor distance over the 480 slow 300-site
loci, declining on the imbalance statistic at 115, with quartet-frequency lengths. Averaging
distances over loci is a consistent species-tree method under the coalescent when the distances
are unbiased, so this is not a science-free shortcut; it fails because the uncorrected distances
carry the long-branch bias on every rate class, and it names the wrong tree in five of the six
long-branch worlds. Its held-out score of 0.197 is the same failure on the sealed worlds. Every
strategy that uses the free gene trees on medium or fast loci stays under 0.34. The probe's own
weakness is that it never re-estimates a gene tree, and re-estimating them with the published
site-rate model is the first thing the task asks for.

In the first version of the worlds, with the fast species at five to eight times the rate, the
same strategy reached 0.838 against a reference of 0.877, because the bias was confined to fast
loci and slow loci gave it the topology and the refusal for free. The maintainer's review of the
first pull request caught that, and the fast species are now at ten to fifteen times the rate on
species branches of 0.4 to 0.7 coalescent units.

## Construction errors, found by the checkpoints and by review

Four before any model saw the task, and three by the review of the first pull request.

- **There was no long-branch attraction.** The first two-fast-species worlds put the fast species
  on ordinary branches, and the uncorrected distance moved zero to half a per cent of gene trees.
  The analytic four-point condition showed why: attraction needs the fast species' sisters to be
  short and the branch between them shorter. The worlds became two cherries joined by a branch of
  0.1 to 0.2 coalescent units, with the fast species one from each cherry.
- **The budget could not buy the answer.** At 80 units no quartet strategy recovered any
  anomaly-zone tree and the refusal statistic did not separate hybrid from tree worlds. A
  sorted-difference statistic also turned out to have a null mean independent of the sample
  size, so it could never be thresholded; it was replaced by the chi-square form. The budget is now
  480 with a catalogue of 2700, and the branches were lengthened until the reference reads them.
- **One anomaly world was not anomalous.** At short branches of 0.12 to 0.2 coalescent units one
  held-out species tree was the most frequent gene tree by simulation. The range is now 0.1 to
  0.16.
- **The reference tied the probe.** At a budget of 400 the reference lost one development world
  by 0.2 quartet votes out of 243, and the probe's best strategy tied it there at 0.779 against
  0.774 while trailing held out at 0.682 against 0.917. Neither reweighting the quartets, mixing
  rate classes nor bootstrap support moved that world off the edge; each flipped a different
  world. The budget was raised to 480.
- **The attraction was measured as the wrong quantity, and held in two worlds of six.** The
  first version of this document reported the raw rate at which free trees join the two fast
  species, ten to thirty per cent, most of which is incomplete lineage sorting that the true gene
  trees carry too. Measured as the excess over the true gene trees of the same loci, the
  attraction was positive in only two of the six worlds on fast loci and absent on slow loci,
  which is exactly where the probe's best strategy was buying. The fast species are now at ten to
  fifteen times the rate on species branches of 0.4 to 0.7 coalescent units, and the excess is
  positive in every world and every class, as the table above shows.
- **Every tree world was said to be in the anomaly zone.** Only the six anomaly-zone worlds are;
  the six long-branch worlds have at most one short internal branch and their species tree is the
  most frequent gene tree. The construction script had filtered to the anomaly worlds while the
  task text, the card and a test named after every tree world claimed otherwise. All of them now
  say exactly what was measured, and the anomaly statistics are reported as probabilities with a
  standard error rather than as a rank.
- **A wrong tree at low confidence escaped the false discovery count.** A tree world named wrong
  counted as a false discovery only above a confidence of 0.5, while a tree on a hybrid world or a
  refusal on a tree world counted unconditionally, so a constant confidence of 0.49 lowered the
  reported rate without changing the score. A wrong tree is now a false discovery whatever its
  confidence, and confidence feeds only the calibration axis.

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
  discovery rate of 1.0 at confidence 0.9 and at confidence 0.2 alike.
- Every source file compiles under Python 3.8 syntax, which is what the evaluation host runs.
