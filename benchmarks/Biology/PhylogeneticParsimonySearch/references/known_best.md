# PhylogeneticParsimonySearch: withdrawal after repair — 2026-09-08

## Status and measured repairs

Independently permuted taxon-to-sequence assignments and input rows using a
separate deterministic RNG. The underlying sequence multiset is unchanged.
Regression coverage includes arbitrary renaming, row reordering, the label-only
block probe and twenty additional fixed generator seeds.

On the five scored instances after repair:

| Solver | Development | Held-out |
| --- | ---: | ---: |
| Label-only blocks | 0.065333826 | 0.006547619 |
| UPGMA reference | 0.938471597 | 0.934821429 |
| UPGMA + NNI | 0.943260352 | 0.943452381 |

The original label-only result was 0.667610. Although label leakage is repaired,
removing the structured labels also makes the caterpillar baseline substantially
worse. The unchanged UPGMA reference now closes almost the entire baseline gap;
NNI adds less than one percentage point on both splits. This family is therefore
withdrawn. No new denominator, weaker reference, or reachable-top claim was added
to force a 0.5–0.8 score. The existing sitewise bound remains a relaxation and is
not asserted jointly attainable.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_phylogenetic_parsimony.py
python scripts/check_task_contribution.py --task Phylogenetics/PhylogeneticParsimonySearch --timeout 300
```

The dedicated regression suite reproduces the repaired contract. Scientific
shortcut comparisons are reported above and in the corresponding PR discussion;
they are not frontier-model calibration draws. `verification/reference_*.py`
is the input-only comparison. For parsimony, `verification/shortcut_probe.py`
and `verification/headroom_probe.py` provide the negative and search controls.

## Remaining evidence requirements

Software validity and deterministic repeats do not establish task difficulty.
The current instance family is withdrawn; baseline/reference tests check valid,
bounded execution, **not** the former 0.5–0.8 admission band. This test change
records a failed calibration rather than relaxing admission requirements.
No blind model draws, long-horizon calibration or external domain confirmation
were performed. Seeds are repository-visible; held-out means omitted from search
feedback, not secret. A replacement must document new instances, reachable
endpoints, all shortcut/ablation measurements and independent confirmation.

## History

The original measurements and the maintainer's September 8 findings remain in
[PR #41](https://github.com/Geniusyingmanji/ScientistsLastExam/pull/41).
The branch preserves the fixes and their regression tests for a future redesign;
closing the current PR does not assert that this scientific subject is unusable.
