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
| `anomaly_rank.py` | the rank of each species tree among gene-tree topologies over one hundred thousand simulated gene trees, at the previous and the current short-branch range | about 10 min |

JSON side outputs go to the system temporary directory.
