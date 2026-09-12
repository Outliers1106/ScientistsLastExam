# PR85: fixed engineering and reference-completeness review

Input PR source: `18f9c06453e943b947817ad66c4bdd9f9de6ebcb`.
No model proposal has been requested in this review.

Engineering changes reset the candidate session at all 18 world boundaries,
reject partially invalid runs, and report false discoveries over valid claims
with explicit denominators. World generation, budgets, valid-world numerical
scores and existing candidate programs are unchanged. Three focused regressions
failed on the original implementation and pass after repair. A separate real
CandidateProxy regression checks module, imported-module and file state across
worlds while retaining state across two queries in each world.

The fixed scientific comparison has six candidate programs, each evaluated twice
on Linux with the original 300-second task timeout:

1. Original baseline.
2. Original independent-background normal-threshold reference.
3. Declared correlation-claims/pivotal-refusal probe.
4. Declared pivotal-without-refusal probe.
5. Shared-background pivotal queries with an exact binomial test, one label per point.
6. The same shared-background test with three labels per point and majority voting.

The last two programs are truth-blind, use only the public problem/query interface,
and spend no more than the published label budget. Both retain delta=1e-3
Bonferroni control. With odd repetition count r and independent label noise eta,
majority error eta_r is P[Binomial(r,eta)>r/2]. For an irrelevant variable the
disagreement count across m backgrounds has null Binomial(m,2 eta_r(1-eta_r)).
Different variables can share noisy background answers: the union bound needs
only the individual null marginals, not independence between variables. Programs
5 and 6 thus test reference completeness without relaxing the declared error rate.

There is no parameter grid, adaptive re-selection, held-out selection or world
retuning in this plan. All six programs and their hashes are frozen before the
first Linux measurement. Private full metrics remain outside Git under 0700/0600
permissions. Public results will contain only split aggregates, counts, source
hashes and explicit failures. Model difficulty remains pending until this
reference-completeness check has a defensible outcome.
