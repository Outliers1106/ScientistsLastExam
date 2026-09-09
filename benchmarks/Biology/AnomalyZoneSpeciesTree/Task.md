# AnomalyZoneSpeciesTree — which species tree, when the commonest gene tree is the wrong one?

## 关系与区别 / How this differs from the nearest tasks in this repository

- **`SystemsBiology/GeneNetworkIntervention`** is the other Biology task in the structure cell. It
  recovers a signed regulatory graph from perturbation experiments and then designs an
  intervention. Here the structure is a **phylogeny**: an unrooted binary tree whose leaves are
  species and whose internal branches carry lengths in coalescent units, and the data are
  sequences, not expression levels.
- **`Algorithm/GraphFromDistances`** recovers a graph from exact shortest-path distances. Here the
  distances are estimates from finite sequences under a model with unpublished rate multipliers,
  the loci disagree with each other by construction, and the tree that fits the mean distance is
  not necessarily the species tree.
- **`Physics/HiddenCouplingNetwork`** infers a coupling graph from dynamics it can drive. There is
  no dynamics to drive here: every locus is an independent draw, and the only choice is which
  loci to buy.
- **`PopulationGenetics/DemographicSFS`** is the nearest in biology, a coalescent inversion, but
  it recovers population-size history for one population from a site-frequency spectrum. This
  task recovers the branching order of eight species, where the coalescent is the reason the
  gene trees disagree with the answer.

No task in the Frontier-Eng catalogue concerns phylogenetics, the multispecies coalescent or
species-tree estimation.

## The question

Eight species, one sampled genome each. A catalogue of loci is on offer, each with a length in
sites, a rate class and a price. Each locus you buy returns an alignment of eight sequences and,
as a convenience, the neighbour-joining gene tree the sequencing centre computed for it. From what
you buy, name the unrooted species tree and the length of each of its five internal branches in
coalescent units, or say that no species tree generated these loci because one species is a
hybrid.

## Four ways to be wrong

- **The commonest gene tree is not the species tree.** Every locus is a draw from the
  multispecies coalescent on the species tree. When two adjacent internal branches are short
  enough, the species tree lies in the *anomaly zone*: a gene tree with a different topology is
  the most probable one, and the majority vote over gene trees converges to the wrong answer as
  the number of loci grows. Every tree world here is in that zone, with the species tree ranked
  third to twenty-seventh among gene-tree topologies by probability. Concatenating the loci and
  building one tree is inconsistent for the same reason. What is not anomalous is any **quartet**:
  for any four species, the species-tree quartet is the most probable gene-tree quartet under any
  branch lengths, and its probability is `1 - (2/3) exp(-t)` where `t` is the length of the
  internal path between the two pairs.
- **Two species evolve fast, and the free gene trees join them.** Two non-adjacent species carry a
  rate multiplier that is not published. Sites evolve with gamma-distributed rates of the
  published shape, and the sequencing centre's tree uses the plain Jukes-Cantor distance, which
  ignores that. The omission compresses long distances more than short ones, so on the faster
  loci the two long branches are drawn together: long-branch attraction, in its distance-method
  form. The bias grows with the substitution rate of the locus, so the loci that resolve short
  branches best are also the ones that mislead most.
- **A hybrid species has no species tree.** In some worlds one species descends from two parents,
  and each locus follows one parent or the other with a fixed inheritance probability. Under any
  species tree the two minority topologies of a quartet are equally probable; a hybrid makes them
  unequal on the quartets that contain it. That imbalance is the only evidence, and it is visible
  only in the data: nothing in the public problem says which worlds are reticulate.
- **Short branches cost loci.** The shortest internal branches are a tenth of a coalescent unit.
  The quartet signal for such a branch is a few per cent of the gene trees, so it takes hundreds
  of loci to read it, and a locus of three hundred slow sites often does not resolve it at all.
  The catalogue prices length, not rate class, and the budget buys a fraction of what is on offer.

The reticulate world is the declining case. A tree world is determinable, and declining it is a
false discovery: it claims a mechanism, hybridisation, that is not there.

## Where the budget actually goes

The budget is 480 units. A locus of 300 sites costs 1, of 800 sites 2, of 2000 sites 5, so the
whole budget is 480 short loci, 240 medium ones or 96 long ones, in any mixture. The catalogue
holds 300 loci of each rate class and length, 2700 in all. Halving the reference's spend costs it
0.35 of the score and a quarter of the budget scores zero on the development split, so the
question is not whether to spend but on what: rate class decides whether the short branches are
resolved and whether the long ones are distorted, and length decides how many independent draws
of the coalescent you see.

## What you implement

```python
def infer_species_tree(problem, sequence):
    ...
    return {"verdict": "tree", "newick": "((A:1,B:1):0.12,...);", "confidence": 0.8}
```

### `problem` — every key you are given

| key | meaning |
|---|---|
| `taxa` | the eight species names, `A` to `H` |
| `locus_budget` | how many cost units you may spend on this world (480) |
| `catalogue` | one entry per locus on offer: `{"locus", "sites", "rate_class", "cost"}`; `locus` is the integer you pass to `sequence`, `sites` is 300, 800 or 2000, `rate_class` is `slow`, `medium` or `fast`, `cost` is 1, 2 or 5 |
| `gamma_shape` | shape of the gamma distribution of site rates (0.5) |
| `coalescent_model` | prose: one lineage per species, independent loci under the multispecies coalescent, branch lengths in coalescent units |
| `sequence_model` | prose: Jukes-Cantor with gamma-distributed site rates, three rate classes a factor of three and ten apart, and unpublished per-species rate multipliers |
| `gene_tree_estimate` | prose: what `nj_tree` is and what it does not correct for |
| `reticulation_model` | prose: one hybrid species with a fixed inheritance probability, in some worlds |
| `abstain_when` | prose: when the verdict is `reticulate` |
| `answer_format` | prose: what a tree verdict must contain |

### `sequence(locus)`

Charges the catalogue cost of that locus and returns `{"locus", "sites", "rate_class",
"alignment", "nj_tree"}`. `alignment` maps each taxon to a string of `sites` nucleotides;
`nj_tree` is a Newick string with branch lengths in substitutions per site, unrooted, from
neighbour joining on the plain Jukes-Cantor distance. `locus` must be an integer index into the
catalogue. Buying a locus twice returns the same data and is charged once. Spending past the
budget raises and the world scores zero, so count your units.

### What you return

| key | meaning |
|---|---|
| `verdict` | `"tree"` or `"reticulate"` |
| `newick` | required when the verdict is `tree`: an unrooted binary tree on exactly the eight taxa, each once, with a finite non-negative length on every branch; the five internal lengths are read in coalescent units, the pendant ones are ignored |
| `confidence` | in `[0, 1]`; clipped |
| `abstain` | `True` is accepted as a synonym for the reticulate verdict |

Anything malformed scores that world zero. It is never an infrastructure failure.
`sle.contract_lint` is importable and free to call for shape checks; it costs no locus.

## How you are scored

Twelve development worlds: four in the anomaly zone, four with two fast-evolving species, four
with a hybrid. Each contributes a mechanism score:

- **tree world**: the right unrooted topology scores `0.5 + 0.5 * branch_score`, where the
  branch score is the mean over the five internal branches of `exp(-|ln(estimate / truth)|)`,
  so a length off by a factor of two scores 0.5 on that branch. The wrong topology scores 0 and
  is a false discovery when `confidence >= 0.5`. Declining scores 0 and is a false discovery.
- **reticulate world**: declining scores 1 and is a correct refusal; any tree scores 0 and is a
  false discovery.

`combined_score` is the mean over the development worlds, renormalised so that **declining every
world scores exactly 0.0**. Naming one fixed tree everywhere also scores 0.0.

Reported separately, never averaged into one number:

`development_topology_rate` · `development_branch_length_score` ·
`development_false_discovery_rate` · `development_correct_refusal_rate` ·
`development_discovery_coverage` · `development_confidence_calibration` ·
`development_mean_loci_bought` · `development_mean_budget_spent` · `development_raw_mechanism`

A sealed held-out set of six further worlds, two of each kind, is scored too, under the same
keys with the `heldout_` prefix, and is not visible to a searcher. `per_instance` carries one row
per world.

## What each competence is worth

Ablating the reference, one choice changed at a time:

| strategy | score | topology | branch lengths | false discovery | refusal | coverage | held out |
|---|---|---|---|---|---|---|---|
| gamma-corrected distance gene trees from 240 medium loci, quartet consensus, quartet-frequency lengths, minority-imbalance refusal | **0.877** | 1.00 | 0.75 | 0.00 | 1.00 | 1.00 | 0.929 |
| same, but the sequencing centre's uncorrected gene trees | 0.316 | 0.38 | 0.26 | 0.42 | 1.00 | 0.75 | 0.467 |
| same, slow loci instead of medium | 0.796 | 0.88 | 0.72 | 0.08 | 1.00 | 1.00 | 0.683 |
| same, fast loci instead of medium | 0.542 | 0.75 | 0.58 | 0.25 | 0.75 | 1.00 | 0.904 |
| same, 300-site loci (all 300 of them) instead of 800 | 0.553 | 0.75 | 0.61 | 0.25 | 0.75 | 1.00 | 0.874 |
| same, 2000-site loci (96 of them) instead of 800 | 0.290 | 0.62 | 0.46 | 0.42 | 0.50 | 1.00 | 0.617 |
| same, never declining | 0.377 | 1.00 | 0.75 | 0.33 | 0.00 | 1.00 | 0.429 |
| same, constant branch lengths of 0.1 | 0.744 | 1.00 | 0.49 | 0.00 | 1.00 | 1.00 | 0.705 |
| same, greedy consensus of the gene trees instead of quartets | 0.452 | 0.50 | 0.40 | 0.33 | 1.00 | 1.00 | 0.694 |
| same, half the budget | 0.528 | 0.88 | 0.68 | 0.25 | 0.50 | 1.00 | 0.913 |
| same, a quarter of the budget | 0.000 | 0.50 | 0.38 | 0.67 | 0.00 | 1.00 | 0.390 |
| uncorrected gene trees and fast loci | 0.100 | 0.25 | 0.20 | 0.58 | 0.75 | 0.38 | 0.223 |
| uncorrected gene trees, never declining | 0.000 | 0.38 | 0.26 | 0.75 | 0.00 | 1.00 | 0.404 |
| concatenation of the longest fast loci, plain Jukes-Cantor, never declining (baseline) | 0.000 | 0.12 | 0.01 | 0.92 | 0.00 | 1.00 | 0.000 |
| declining everything | 0.000 | — | — | 0.67 | 1.00 | 0.00 | 0.000 |

Every row costs something real: trusting the free gene trees costs 0.56, never declining 0.50,
greedy consensus 0.43, the wrong locus length 0.32 or 0.59, the wrong rate class 0.08 or 0.34,
half the budget 0.35, and constant branch lengths 0.13. The reference's branch-length score is
0.75, because it inverts the quartet frequencies without correcting for gene-tree estimation
error, and that is where most of what remains lives.

**Low-dimensional shortcuts do not solve this task.** A sweep of 1260 strategies that take the
sequencing centre's gene trees as given, choose one catalogue cell, build the topology by greedy
consensus, the most frequent gene tree, neighbour joining on the mean distance or quartet
consensus, decline on a discordance or imbalance threshold or never, and set the branch lengths to
a constant or from the quartet frequencies, reaches **0.838** on the development split and 0.604
held out, against the reference's 0.877 and 0.929 and a ceiling of 1.0. Without quartet consensus
and quartet lengths the best is 0.744 and 0.509; on fast loci 0.331; never declining 0.338.

## Rules

- Only edit `solution.py`; keep `infer_species_tree(problem, sequence)`.
- `sle.contract_lint` is importable and free to call for shape checks. It costs no locus.
- Do not read `verification/` or `frontier_eval/`.
