# PR83 independent admission review, 2026-09-12

This review starts from author head `c74e9c6188a156e30c158602d2232e5460749851`.
It preserves every failed replay recorded on 2026-09-11. Those 18 invalid results
describe earlier solvers, and are not failures of the current revised-simplex
implementation. Author builder records remain author provenance.

## Fixed numerical and sandbox evidence

The frozen five-program plan ran twice per program using the real Linux
CandidateProxy, Python 3.8.10, NumPy 1.24.4 and SciPy 1.10.1, with all numerical
thread limits one and the original 300-second timeout. Every one of the ten
evaluations was valid in all 18 worlds. Each pair's complete metrics was identical.
The reference took 21.56 and 22.15 seconds. There were no changed worlds, new grid,
model-generated repair, dependency installation or relaxed sandbox rule.

| Fixed method | Development | Held out |
|---|---:|---:|
| Baseline | 0 | 0 |
| Full reference | 0.5688535714 | 0.50286675 |
| Detection-only shortcut | 0.3847268571 | 0.298126 |
| No order-statistic rows | 0.5071945714 | 0.4449615 |
| No empty-queue rows | 0 | 0 |

The declared shortcut stays below the unchanged ten-percent reference margin.
Both substantive ability ablations lower development and held-out scores. The
author's 29-strategy sweep is not relabelled as a new independent sweep.

Two additional runs through `frontier_eval/run_eval.py` reproduced the reference's
full private metrics exactly. Public wrapper output contained only permitted
search metrics. All 21 task/runtime tests passed on Linux without skips, including
the real global-state and private-tmpfs reset check across every world.

The retained 820-row, 26-variable development LP that previously hit an iteration
limit was independently checked. Exact rational convex-combination certificates
establish that all 686 removed rows follow from the 134 retained rows. Independent
HiGHS solves on all rows and reduced rows, and the current revised simplex on
reduced rows, all report infeasibility. No optimum or scientific score is inferred
from that infeasible LP. This is an engineering check of one retained failure,
separate from the complete valid sandbox evaluations.

The initial contribution gate at `c74e9c618` passed structure and runtime but
retained `expected_score: null` and therefore an incomplete shortcut guard.
Commit `91d82c1c36d47ad0a2d386b9bcbdbdccd3f410d0` separately records the independently
measured declarations. A fresh complete gate at that clean revision passed all
16 checks, including runtime and shortcut guard. Each gate made 11 task sandbox
calls; intentional malformed-candidate trials are not valid scientific proposals.

## Novelty and scope

Both external catalogs were retrieved again: the [paper's 47-task appendix](https://arxiv.org/html/2604.12290v1#A1) and
the [repository's 78 rows / 84 expanded entries](https://github.com/Einsia/Frontier-Engineering/blob/e3fa29c193356af2ce1ec8b3d23ab1a2e2410071/TASK_DETAILS.md) at
`e3fa29c193356af2ce1ec8b3d23ab1a2e2410071`. Their hashes match the retained catalog
identities. Comparing task objective, artifact and oracle found no same problem
class. This is catalog coverage, not an audit of every external implementation.
The per-catalog scope and nearest relationships are recorded in
`pr83_frontier_novelty_2026-09-12.json`.

## Frozen first proposals

The independent plan `c5f70345fefbbd9101f6ba9b1126908199b2edf7ad13a1691c5ad54a0213cc0f`
uses clean source `91d82c1c36d47ad0a2d386b9bcbdbdccd3f410d0`, three seed labels
0/1/2, one greedy-rewrite proposal per cell, selection-blind feedback, and
gpt-5.6-sol over Responses with reasoning high, maximum output 16384 tokens and
300-second request/evaluation timeouts. Seed labels do not control provider RNG.
No proposal will be dropped or adaptively retried. All candidate sources, complete
metrics, durable receipts and provider records remain outside git in mode-0700
directories and mode-0600 files.

The completed plan produced three independent proposals, each valid in all 18
worlds. The source-bound ledger review verifies every retained candidate, request
and durable result receipt. Returned-response usage totals 53565 tokens (10821
input, 42744 output); price and unreturned-provider-request usage are unavailable.

| First proposal label | Development | Held out | Dev FDR | Held FDR | Dev / held refusal | Dev / held discovery coverage |
|---|---:|---:|---:|---:|---|---|
| 0 | 0.4496327143 | 0.51835475 | 0.125 | 0 | 0.8 / 1 | 1 / 1 |
| 1 | 0.7907742857 | 0.57135925 | 0 | 0 | 1 / 1 | 1 / 0.75 |
| 2 | 1 | 1 | 0 | 0 | 1 / 1 | 1 / 1 |

**HOLD: D16 failed.** Two of three first proposals exceed the frozen development
reference; one attains the maximum on both splits. These are complete valid
scientific evaluations, not runtime errors or blanket refusal. The first
proposal's held-out score also exceeds the reference; held-out results remain a
separate axis and were not used to change the frozen development threshold.

## The public timing contract defeats the claimed non-identifiability

The contract gives the experimenter exact true start time zero, exact spacing and
exact advancement by each exchange and wait. The evaluator constructs each send
at `t = state.now + SPACING * arange(n)` and returns `T1 = C_i(t) + noise`.
Either endpoint may initiate an exchange. For each sender, therefore,
`T1 - t = theta_i + s_i*t + noise` is a directly observed linear regression. No
unknown propagation delay enters that equation. With sends at multiple times,
the clock offset and skew are point-identifiable under the disclosed sampling
model, instead of being constrained only by propagation/asymmetry bounds.

This derivation uses the public interface and evaluator, without hidden data or
candidate source disclosure. Inspection of the privately retained perfect first
proposal confirms that it uses known send times for clock regression, alongside
separate delay/model checks. It does not obtain evaluator code through the sandbox.
The source locations are `Task.md`'s exchange contract and `verification/evaluator.py`
lines 219, 242, 250 and 255; the timestamps are generated at the known schedule.

A repair requires redesigning the time-control and observation contract and then
re-deriving the identifiability anchor. Raising the reference or changing worlds,
thresholds or solver settings after these draws would not repair this defect.
All failed and successful original records are retained. The published engineering
repair and passing C/runtime checks do not authorize task admission.

The only task-package edits after this comparison document its result. A separate
fixed replay of all three retained programs twice verifies compatibility of those
documentation edits; it does not add model samples. Its completed receipt summary
is recorded separately before final publication. Domain review and certification
remain pending.
