# ClockSyncInversion: how far off is every clock in this network, and is the model even right?

## 关系与区别 / How this differs from the nearest tasks in this repository

- **`Sensors/IMUBiasCalibration`** sits in the same cell. It also runs a charged calibration
  campaign and then either publishes a calibration or declines because the instrument does
  something the model cannot absorb. There the unknowns are point parameters of one instrument.
  Here they are the offsets of a whole network of clocks, and the data can never pin them to a
  point, only to an interval whose width the physics of timing fixes.
- **`Geophysics/GravityInversion`** is also an inversion that many answers fit equally well.
  There the non-uniqueness is resolved by a source library. Here it is part of the answer: you
  report the interval, and the score compares your width with the width that unlimited data
  could not shrink.
- **`ParticlePhysics/DiscrepantMeasurements`** and **`Turbulence/WallClosureDiscovery`** ask
  when the evidence cannot support an answer. Here the model can fail in several ways, and each
  failure is a few microseconds hidden in jitter. It shows only jointly across links and time.

No other task in this repository concerns clocks, networks or timing measurements.

## The question

Every computer keeps time with a crystal oscillator that runs slightly fast or slow. Two
computers learn how far apart their clocks are by exchanging timestamped messages, as NTP and
PTP do. The difficulty is old and fundamental. A message's travel time is unknown. Only the
round trip can be measured, and no timing measurement can tell a clock offset apart from a
difference between the two directions' delays. What the data do determine is how the clocks
drift. The offsets they fix only within bounds, and those bounds come from what is known about
the delays.

You are given a small network. Node 0 keeps true time. You may run two-way timestamp exchanges
over its links, under a budget and within a time horizon. Then bound the offset of every other
clock at the start and at the end of the horizon, or decline. Declining is the right answer when
the exchanges cannot have come from the model described below.

## What makes it hard

- **Offsets hide behind asymmetry.** An interval narrower than the published delay bounds allow
  is a guess. A guess that misses the true offset costs a world.
- **Most probes queue.** A one-way delay is its propagation delay plus queueing. Only the probes
  that find the queue empty carry the propagation floor. At the worst time of day that is one in
  twenty, and timestamp jitter hides which probes they were.
- **Clocks drift.** Rates differ from true time by up to 30 ppm, which is 0.43 s over the horizon.
  Where in time you spend the budget decides how well the rates, and so the end offsets, are
  pinned.
- **The model can be wrong.** In some worlds some part of the model below is false. The exchanges
  then fit no clocks and delays of the model, and every interval is meaningless. The violations
  are a few microseconds, under the same jitter as everything else.

There are twelve development worlds with 5 to 8 nodes and 6 to 11 links. In seven the model
holds and in five it does not. Nothing in the public problem says which worlds are which.

## What you implement

```python
def identify(problem, exchange, wait):
    ...
    return {"verdict": "offsets",
            "intervals": {"1": [[lo_0, hi_0], [lo_H, hi_H]], "2": [...], ...},
            "confidence": 0.8}
```

### `problem`: every key you are given

| key | meaning |
|---|---|
| `nodes` | the number of nodes N, 5 to 8; node 0 keeps true time |
| `links` | the links as pairs `[a, b]` with a < b; the network is connected |
| `lower_bound_s` | `{"a->b": L_ab}` for both directions of every link: the propagation delay from a to b is at least L_ab seconds |
| `calibrated_links` | pairs `[a, b]` whose two directions' propagation delays differ by at most `asymmetry_bound_s` |
| `asymmetry_bound_s` | 1e-6 |
| `timestamp_sigma_s` | per node, the standard deviation of the Gaussian jitter on every timestamp it writes: 2e-7, or 2e-6 on nodes that stamp in software |
| `empty_queue_probability_min` | 0.05, see the delay model |
| `probe_budget` | 20000, the exchanges you may run in this world |
| `horizon_s` | 14400, the time you may use |
| `probe_spacing_s` | 0.25, the time between two probes of one call |
| `epochs_s` | `[0.0, 14400.0]`, the true times at which you bound the offsets |
| `clock_model`, `delay_model`, `timestamp_model`, `exchange_model`, `answer_format` | prose: the model and the contract below |

### The model

- **Clocks.** At true time t node j reads C_j(t) = t + theta_j + s_j t, with theta_j and s_j
  unknown constants. Its offset is x_j(t) = C_j(t) - t. Node 0 has theta_0 = s_0 = 0.
- **Delays.** A message from a to b takes p_ab + Q. The propagation delay p_ab is a constant with
  p_ab >= L_ab. Queueing Q is non-negative. In every direction and at every moment Q = 0 with
  probability at least `empty_queue_probability_min`, independently from probe to probe. On a
  calibrated link, |p_ab - p_ba| <= `asymmetry_bound_s`. Nothing else about Q is promised.
- **Timestamps.** Every timestamp node j writes is C_j at that moment plus independent Gaussian
  noise with standard deviation `timestamp_sigma_s[j]`.

### `exchange(i, j, n)` and `wait(seconds)`

Time starts at 0. `exchange(i, j, n)` runs n two-way exchanges over the link between i and j,
either end first. The k-th starts at the current time plus k times `probe_spacing_s`. Node i
writes T1 and sends, j writes T2 when it receives and T3 when it replies, and i writes T4 when
the reply arrives. The call returns an n x 4 array of (T1, T2, T3, T4), each on its writer's own
clock, and advances the time by n times `probe_spacing_s`. `wait(seconds)` advances the time.
A call that would run past `probe_budget` exchanges in all, start a probe past `horizon_s` or
wait past it raises, and the world scores zero, so count what you spend. Nodes and n are
integers, n >= 1, and there must be a link between i and j.

### What you return

| key | meaning |
|---|---|
| `verdict` | `"offsets"` or `"no_model"` |
| `intervals` | required for offsets: every node j in 1..N-1, as `"j"` or j, maps to two finite intervals `[lo, hi]` in seconds, bounding x_j at `epochs_s[0]` and at `epochs_s[1]` |
| `confidence` | in `[0, 1]`; clipped |
| `abstain` | `True` is accepted as a synonym for the no_model verdict |

Anything malformed scores that world zero, and so does a malformed call to `exchange` or `wait`.
It is never an infrastructure failure.

## How you are scored

The oracle knows the clocks and the delays. With unlimited noiseless exchanges the floor of every
direction, x_b(t) - x_a(t) + p_ab, would be known exactly at every t. That pins every rate but not
the offsets. Any other explanation shifts the offsets by constants and the propagation delays to
match. The shifts are limited only by p_ab >= L_ab and by the asymmetry bounds on calibrated
links. For every clock the oracle solves two linear programs for the range of offsets that
remains. That range is the identifiable width w*_j, the same at both epochs.

- **A world where the model holds.** Intervals that all cover the true offsets score the mean,
  over clocks and epochs, of min(1, w*_j / width). Declining scores 0. If any interval misses its
  true offset, the answer is a false discovery and scores -1.
- **A world where the model does not hold.** Declining scores 1 and is a correct refusal. Any
  intervals are a false discovery and score 0.

A false discovery therefore costs one world wherever it is made. `combined_score` is the mean over
the development worlds, renormalised so that **declining every world scores exactly 0.0**.
Negative values clip to 0. `confidence` feeds only the calibration axis and cannot turn a false
discovery into anything else.

Reported separately, never averaged into one number:

`development_interval_sharpness` · `development_false_discovery_rate` ·
`development_correct_refusal_rate` · `development_discovery_coverage` ·
`development_confidence_calibration` · `development_mean_probes_used` ·
`development_raw_mechanism`

A sealed held-out set of six further worlds, four where the model holds and two where it does not,
is scored too, under the same keys with the `heldout_` prefix, and is not visible to a searcher.
`per_instance` carries one row per world.

## Where the scale sits

The reference scores 0.568 on the development split and 0.501 held out. It declines every world
where the model fails and makes no false discovery. Over sixteen seeds for the queueing and the
jitter, this one and fifteen re-drawn, it averages 0.567 and 0.501 and never makes a false
discovery. It is not the ceiling. Sharper intervals from the same budget exist.

The baseline in `solution.py` scores 0.000. It is NTP with a rate: the minimum delay in each
direction, half their difference, a line through twelve rounds, five microseconds either side. It
misses a true offset in every world where the model holds and claims every world where it does
not. Declining everything scores 0.000.

A sweep of 29 low-effort strategies was scored on the same worlds. They answer blind, take NTP
midpoints with and without a rate at margins from 1 us to 1 ms, bound each round with Cristian's
method, fit a lower envelope to the minimum delays without checking the model, or check the model
first and then answer from a linear program over the minimum delays alone. The best of them
reaches 0.385 on the development split, 68 per cent of the reference, and 0.298 held out.

## Rules

- Only edit `solution.py`; keep `identify(problem, exchange, wait)`.
- NumPy, SciPy and the standard library only. Deterministic CPU code.
- `sle.contract_lint` is importable and free to call for shape checks. It costs no exchange.
- Do not read `verification/` or `frontier_eval/`.
