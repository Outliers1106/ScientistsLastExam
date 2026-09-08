# BatchEffectDiscovery: withdrawal after repair — 2026-09-08

## Status and measured repairs

Fixed the four-gene construction constant: supported worlds now sample 2–8 genes
and include smaller effects. FDR now counts claims in every world, with an explicit
claim denominator; held-out scores use the same normalization rule as development,
and both splits report refusal denominators.

These repairs do not rescue the scientific design. After the repair the reference
scores 0.443361109 development / 0 held-out. Removing the batch column gives
0.443361109 / 0 as well (roundoff below 1e-15). Before the repair both methods
scored 0.511368501 development. The balanced follow-up layout still makes batch
adjustment unnecessary, and `available_cells` identifies the confounded world
without measurements. Variable gene counts remove one shortcut but do not establish
admission or headroom. The family is withdrawn pending a new observation design.

## Reproduction

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_batch_effect.py
python scripts/check_task_contribution.py --task Genomics/BatchEffectDiscovery --timeout 300
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
[PR #38](https://github.com/Geniusyingmanji/ScientistsLastExam/pull/38).
The branch preserves the fixes and their regression tests for a future redesign;
closing the current PR does not assert that this scientific subject is unusable.
