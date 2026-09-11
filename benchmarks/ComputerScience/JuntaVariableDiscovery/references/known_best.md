# JuntaVariableDiscovery: known best

## Reference (truth-blind): pivotal pairs with a family-wise margin

`verification/reference_influence.py`. It reads only the public problem and what `query` returns,
and it works in three stages.

- **Pivotal pairs.** The budget is split evenly over the variables, at ninety per cent of it so a
  malformed count never overspends. For each variable it draws the same number of random
  backgrounds, queries each background and its single-bit flip, and records how often the two
  answers disagree. Flipping a variable against a fixed background is the only way to read its
  influence through a parity or a gate.
- **The noise floor.** Two answers, each flipped with probability eta, disagree with probability
  base = 2 eta (1 - eta) when the variable is irrelevant, and base + Inf (1 - 2 eta)^2 when its
  influence is Inf. So an irrelevant variable's disagreement rate sits at base up to binomial noise,
  and a relevant one sits above it by an amount that shrinks with the noise.
- **The family-wise test.** A variable is named relevant when its disagreement rate is above base by
  more than z standard errors, z the normal quantile of delta / n, so the chance any irrelevant
  variable is named is at most delta by a union bound over the n candidates. When more than
  max_junta_size variables clear the bar, no junta of the published size explains the answers, and
  it declines.

delta is 1e-3 per world.

| split | score | recall | false discovery | refusal | coverage | queries per world |
|---|---|---|---|---|---|---|
| development | 0.548 | 0.548 | 0.00 | 1.00 | 1.00 | 17976 |
| held out | 0.604 | 0.604 | 0.00 | 1.00 | 1.00 | 17976 |

World by world on the development split, the recall runs from 0.417 to 0.750. The relevant
variables it misses are the ones with the smallest influence, the parity and gated-parity
variables, which sit near the edge of what 20000 queries at a 0.10 noise rate can certify above the
family-wise margin. All five worlds that are not juntas are declined, because more than
max_junta_size variables clear the bar.

What it leaves on the table, by design:

- **Its margin is a union bound.** The chance that any of the n candidate variables is named in
  error is held below delta = 1e-3. The disagreement rate of a truly irrelevant variable sits
  tightly at the noise floor, so a tighter margin loses little safety and certifies more of the
  weak relevant variables.

The recoverable headroom is the same test with the margin at three standard errors instead of the
family-wise 3.97, `.research/junta/summary.py`. It scores 0.764 on the development split and 0.783
held out, with no false discovery over twelve run seeds. The cliff is close: at 2.75 standard errors
it scores 0.756 on the graded seed but makes false discoveries on five of the twelve seeds and its
development score collapses to 0.321 on the worst of them; at 2.5 it makes eight.

## Model draws

None. The task was built on a machine without a model endpoint. The card records
`calibration_evidence_status: missing`; the frontier draw and the global evidence refresh are owed
before certification.

## Baseline: degree-one correlation

`solution.py`. One random point per query over the whole budget, each variable's degree-one Fourier
coefficient with the answer, a family-wise margin from zero, decline when more than max_junta_size
clear it.

| split | score | raw | recall | false discovery | refusal | coverage |
|---|---|---|---|---|---|---|
| development | 0.000 | 0.194 | 0.333 | 0.42 | 0.00 | 1.00 |
| held out | 0.000 | - | 0.333 | - | 0.00 | 1.00 |

It recovers the four variables of the additive block, whose degree-one correlation is nonzero, and
so its recall on the junta worlds is 0.333. It recovers none of the seven variables that act only
through a parity or a gate, whose degree-one correlation is exactly zero. A parity of many
variables also has every degree-one coefficient zero, so correlation clears no variable and returns
an empty junta rather than declining, which is a false discovery on every world that is not a
junta. It makes a false discovery in seven of the eighteen worlds.

## Difficulty ladder

`.research/junta/summary.py`, one reference choice changed at a time, each over twelve run seeds.
The graded seeds are the evaluator's own; the mean re-draws every world's label-noise seed. False
discoveries count named irrelevant variables in the junta worlds and claims in the others.

| strategy | development | held out | false discoveries over twelve seeds |
|---|---|---|---|
| **reference**: family-wise margin, z = 3.97 | **0.545** | 0.573 | 0 / 216 |
| margin z = 3.5 | 0.650 | 0.667 | 0 / 216 |
| margin z = 3.25 | 0.696 | 0.726 | 0 / 216 |
| **headroom**: margin z = 3.0 | 0.764 | 0.783 | 0 / 216 |
| margin z = 2.75 | 0.756 | 0.773 | 5 / 216 (dev collapses to 0.321) |
| margin z = 2.5 | 0.757 | 0.828 | 8 / 216 |
| three-quarters of the budget, family-wise margin | 0.469 | 0.476 | 0 / 216 |
| half the budget | 0.428 | 0.439 | 0 / 216 |
| a quarter of the budget | 0.387 | 0.366 | 0 / 216 |

The margin does all the work above the refusal. The family-wise margin at delta / n is a union
bound and is deliberately loose; tightening it to three standard errors recovers the weak relevant
variables and lifts the score to 0.764 with no false discovery, and below 2.75 standard errors the
false discovery rate climbs fast. Under three-quarters of the budget the same full pipeline scores
0.469, so the budget is a real constraint, not slack; that row is the method under-resourced, not a
shortcut.

## Shortcut probe

`.research/junta/summary.py`, ten strategies that do not run the full-budget pivotal test with a
family-wise refusal, in four families, each over twelve run seeds:

| family | strategies | best development | held out | note |
|---|---|---|---|---|
| blind: decline everything, name all twelve, name a random twelve | 3 | 0.000 | 0.000 | the two claiming strategies make a false discovery in every world |
| degree-one correlation at a family-wise margin or a hand-set 0.02 and 0.05 | 3 | 0.000 | 0.000 | correlation cannot see a parity or a gated variable and cannot decline a big parity |
| the pivotal test without the refusal | 1 | 0.000 | 0.073 | it claims every world that is not a junta |
| the pivotal test with a fixed floor margin of 0.02 or 0.03, or correlation for the claims and a pivotal pass only for the refusal | 3 | **0.333** | 0.278 | it makes false discoveries; the only clean strategies decline everything |

The best strategy reaches 0.333 on the development split, 61 per cent of the reference, and it is
not clean: every strategy that makes no false discovery over twelve run seeds scores zero by
declining everything. Correlation-only methods score zero because they miss the parity and gated
variables and claim every big parity. The pivotal test without the refusal claims every world that
is not a junta. A fixed floor margin either misses the weak variables or names irrelevant ones.
Nothing short of the full-budget pivotal test with a calibrated refusal both recovers the weak
relevant variables and declines the big parities.

## Construction errors caught on the way

- The worlds that are not juntas were first an additive threshold of many equal-weight variables.
  Its low-influence variables the reference could not certify, so it failed to decline and made a
  false discovery. They are now a parity of more than the junta size, where every variable has
  influence one and the reference declines by counting significant variables. kappa*, the excess
  influence any junta of the published size must leave uncovered, is the parity size minus twelve,
  from 2 to 6, in every unsupported world and zero in every supported one.
- The budget was first 40000 and the noise rate 0.05. The reference recovered every relevant
  variable and scored near one. The budget is now 20000 and the noise 0.10, so the weakest relevant
  variables sit at the edge of certifiability and the reference scores 0.545.
- A reference at a fixed confidence level was beaten with no false discovery by the same test at a
  tighter margin. The reference is now the family-wise margin at delta / n, and the tighter margin,
  three standard errors, is documented above as the headroom.

## Robustness

`.research/junta/summary.py` re-runs the reference with every world's label-noise seed shifted by
7919 per shift. Over twelve shifts it averages 0.545 on the development split, from 0.500 to 0.595,
and 0.573 held out, from 0.521 to 0.625, with no false discovery in 216 world-runs and every world
that is not a junta declined. The headroom margin averages 0.764 and 0.783 over the same twelve
shifts with none.

`tests/test_junta_variable_discovery.py` checks that a campaign's flip rate follows the published
model against the true function. It checks the recorded relevant sets and kappa*, zero in every
junta world and at least two in every world that is not. It checks that a named irrelevant variable
costs a world, that declining everything scores 0.000 in both forms and naming all twelve or a
random twelve scores 0.000, and that two evaluations of the reference agree key for key. It checks
malformed candidate shapes, including overspending, a patched budget, malformed relevant lists and
malformed calls. All of them score valid 0, combined 0 and feasibility 0 without raising. A
reference evaluation takes under two seconds in process. `.research/junta/summary.py` recomputes
every number on this page.

Not done: the sandbox half of `scripts/check_task_contribution.py` (no Bubblewrap on the build
machine), a frontier-model draw, the global evidence refresh.
