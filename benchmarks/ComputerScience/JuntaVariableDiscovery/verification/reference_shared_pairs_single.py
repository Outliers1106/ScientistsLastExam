"""Truth-blind shared-background pivotal test with exact Bonferroni control.

Reusing a background's noisy labels across variable flips saves queries. Repeated
labels are combined by majority vote. Different variables' tests need not be
independent: the Bonferroni union bound uses only each binomial null marginal.
Every background/flip batch and repetition is charged by the public query API.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import binom

REPETITIONS = 1
DELTA = 1e-3


def identify(problem, query):
    n = problem['n']
    repetitions = REPETITIONS
    per = problem['budget'] // ((n + 1) * repetitions)
    if per < 1:
        return {'verdict': 'no_junta', 'confidence': 0.0}
    rng = np.random.default_rng(12345)
    points = rng.integers(0, 2, size=(per, n))

    def labels(batch):
        observed = np.asarray(query(np.repeat(batch, repetitions, axis=0)))
        return observed.reshape(per, repetitions).sum(axis=1) > repetitions // 2

    background = labels(points)
    noise = float(binom.sf(repetitions // 2, repetitions, problem['noise_rate']))
    null = 2.0 * noise * (1.0 - noise)
    significant = []
    for index in range(n):
        flipped = points.copy()
        flipped[:, index] ^= 1
        disagreements = int(np.count_nonzero(background != labels(flipped)))
        if binom.sf(disagreements - 1, per, null) <= DELTA / n:
            significant.append(index)
    if len(significant) > problem['max_junta_size']:
        return {'verdict': 'no_junta', 'confidence': 0.9}
    return {'verdict': 'junta', 'relevant': significant, 'confidence': 0.6}
