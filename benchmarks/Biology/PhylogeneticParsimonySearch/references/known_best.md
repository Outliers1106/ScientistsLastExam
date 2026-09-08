# PhylogeneticParsimonySearch: scoring evidence — 2026-09-06

## Scientific endpoints

Zero is the caterpillar tree. One is the sum across sites of distinct-state
count minus one, a lower bound on every tree's Fitch cost. This relaxation may
be unattainable jointly across sites; it is not a published record. Scoring is
now explicitly **clipped**, and average linkage is no longer the endpoint.
The unchanged input-only NNI probe in `verification/headroom_probe.py` improves
both development and held-out scores, which remain bounded by one.

## Reproduction and measurements

The legal `solution.py` baseline scores **0.000000**, valid=1. The input-only
comparison reference is `verification/reference_search.py`.
Reproduce using:

```sh
python -m sle eval --allow-uncertified --task Phylogenetics/PhylogeneticParsimonySearch \
  --candidate benchmarks/Biology/PhylogeneticParsimonySearch/verification/reference_search.py --timeout 300
```

Baseline and reference were validated through the Linux candidate sandbox on
implementation commit `3b62c02`; baseline score was exactly zero and both were valid.

| Solver | Development normalized score | heldout_score | Valid |
| --- | ---: | ---: | ---: |
| Original reference | 0.67836236 | 0.65117581 | 1 |
| Public-input headroom probe | 0.70492678 | 0.69624860 | 1 |

The discovery held-out column is raw scientific quality, not the normalized
development scale. Optimization held-out scores use the same normalization as
development. All reference algorithms are unchanged by this calibration.
Pre-calibration score measurements do not describe this revision.

## Limits and provenance

These original procedural worlds are repository-visible; held-out means excluded
from search feedback, not server-secret. No external datasets or code are
redistributed. Model simplifications and nearest-task overlap are in `Task.md`.
Precision/headroom measurements do not establish expert difficulty: strong
classical comparisons, frontier draws, long-horizon search and external domain
review remain pending. The task stays **candidate**.

Scientific sources: doi:10.1111/j.1096-0031.1999.tb00277.x, doi:10.1186/s12859-018-2009-1.

## Admission review — 2026-09-08

Maintainer probes: a label-only block tree scored 0.667610 versus reference 0.678362; multistart SPR reached about 0.732 development and 0.709 held-out. Label leakage and unreachable normalization headroom remain admission blockers.

These findings are from the [maintainer review](https://github.com/Geniusyingmanji/ScientistsLastExam/pull/41) except the explicitly identified independent geometry reproduction. They supersede any earlier suggestion that a low reference score alone establishes useful headroom. Scientific instance/normalization revisions remain pending; passing software tests does not resolve these blockers. No frontier-model draws were performed in this follow-up.
