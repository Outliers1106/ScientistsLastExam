# PR85 independent review: engineering repaired, scientific admission held

Reviewed author head: `18f9c06453e943b947817ad66c4bdd9f9de6ebcb`.
Measured repair source: `0a6857652eb3c9c23e7cbdb3af3864da005d04c0`.
The six-program, two-repeat plan was frozen before execution; its SHA-256 is
`e59e03a81a1428855bd5a91d3a6dffd61503c8684801c76485d54db3c2be5089`.
The public scalar projection is `pr85_fixed_comparison_2026-09-12.json`.
All 12 real Linux trusted-driver/CandidateProxy evaluations were valid in all
18 worlds, with no infrastructure failure. Each pair's complete private metrics
was identical. No model proposal, adaptive parameter search, or heldout selection
was performed.

| Fixed method | Development | Heldout | Labels/world |
|---|---:|---:|---:|
| Original baseline | 0 | 0 | 20000 |
| Original reference | 0.547619 | 0.604167 | 17976 |
| Correlation claims / pivotal refusal | 0.333333 | 0.333333 | 19968 |
| Pivotal without refusal | 0 | 0.104167 | 17976 |
| Shared backgrounds, exact test, one label | 0.845238 | 1.000000 | 19981 |
| Shared backgrounds, exact test, majority of three | 0.940476 | 0.979167 | 19923 |

Both new fixed programs have zero observed false discoveries and correct refusal
rate 1 on both splits. They preserve the original delta=1e-3 familywise guarantee:
an irrelevant variable has exact binomial disagreement probability after the
declared independent label noise and majority rule, and Bonferroni control does
not require independence between different variables' tests. Reusing background
labels improves efficiency without additional oracle access or world knowledge.
This comparison therefore does not depend on the previously disputed looser
three-sigma threshold. The construction and derivation were recorded before
measurement in `pr85_admission_plan_2026-09-12.md`.

The engineering repair resets the sandbox between worlds, makes any invalid world
invalidate the headline task result, and defines FDR over valid positive claims
with explicit counts. The previous fraction over all worlds is retained under
`all_world_false_claim_fraction`. World generation, budgets and valid-world scores
are unchanged. Three focused regressions failed before the repair and passed after
it. All four targeted Linux regressions, including the real per-world sandbox
reset, pass (4 passed, 0 skipped, 9.56 seconds). The optional original full oracle suite was stopped by the operator after the
scientific hold was established; it remained incomplete and is neither counted
as passed nor treated as a scientific failure. Its partial log and stop receipt
remain in private operator storage.

**Decision: HOLD.** CONTRIBUTING C12/C14 require a capable complete reference and
hardening when a standard probe exceeds it. The original reference is materially
incomplete. The higher fixed scores are evidence of that omission, not a new
admitted reference, a model-difficulty result, or proof of global saturation.
The new programs and engineering fix are supplied for review; do not weaken them
or choose hidden worlds to preserve the old admission threshold. Reassess the
complete reference and task difficulty before a new immutable D16 draw.

New raw full metrics, source copies and receipts remain outside Git in owner-only
storage. This report contains aggregates and hashes only; historical author
construction records retain their original provenance.
