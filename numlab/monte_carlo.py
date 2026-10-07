"""Monte Carlo estimation and its O(1 / sqrt(N)) convergence.

The law of large numbers makes the sample mean converge to the expectation,
and the central limit theorem gives the error size: about sigma / sqrt(N).
"""

import numpy as np


def estimate_pi(n: int, rng: np.random.Generator) -> float:
    """Throw n random points in the unit square and count those inside the quarter disc."""
    pts = rng.random((n, 2))
    inside = np.sum(pts[:, 0] ** 2 + pts[:, 1] ** 2 <= 1.0)
    return 4.0 * inside / n


def mc_integral(f, a: float, b: float, n: int, rng: np.random.Generator):
    """Monte Carlo estimate of the integral of f on [a, b] with a 95% confidence interval."""
    x = a + (b - a) * rng.random(n)
    y = (b - a) * f(x)
    mean = float(np.mean(y))
    half_width = 1.96 * float(np.std(y, ddof=1)) / np.sqrt(n)
    return mean, (mean - half_width, mean + half_width)
