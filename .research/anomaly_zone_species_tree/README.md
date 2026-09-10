# AnomalyZoneSpeciesTree construction scripts

The scripts that produced every number in `benchmarks/Biology/AnomalyZoneSpeciesTree/Task.md`,
`TASK_CARD.yaml` and `references/known_best.md`. Run them from the repository root with the
project interpreter; they import the task's oracle directly and never touch the sandbox.

| script | what it produces | runtime |
|---|---|---|
| `run_both.py` | reference and baseline metrics with one row per world | about 15 s |
| `checks.py` | determinism, twenty malformed candidate shapes, blanket abstention, a fixed tree, the reference capped at 25, 50 and 100 loci, Python 3.8 syntax | about 1 min |
| `ablation.py` | the difficulty ladder (one reference choice changed at a time) | about 10 min |
| `probe.py` | the 1260-strategy shortcut probe on the free gene trees | about 1 min |
| `anomaly_stats.py` | for every tree world, the gene-tree probability of the species-tree topology and of the most frequent topology, their difference and its standard error over `AZ_N` simulated gene trees (default one hundred thousand; `AZ_ONLY=<seed>` restricts to one world) | about 20 min at the default, 80 min at four hundred thousand |
| `lba.py` | long-branch attraction measured as the excess rate at which estimated gene trees join the two fast species over the true gene trees of the same loci, for free and gamma-corrected trees, per world and rate class; the minority-imbalance separation on free versus corrected trees; mean-distance neighbour joining; per-locus split recall | about 25 s |

JSON side outputs go to the system temporary directory.
