"""Polynomial interpolation and the Runge phenomenon.

With equally spaced nodes, the interpolation error of 1 / (1 + 25 x^2) on [-1, 1]
grows with the degree. Chebyshev nodes, which cluster near the ends of the interval,
make the error decrease instead.
"""

import numpy as np


def equispaced_nodes(n: int, a: float = -1.0, b: float = 1.0) -> np.ndarray:
    """n + 1 equally spaced nodes on [a, b]."""
    return np.linspace(a, b, n + 1)


def chebyshev_nodes(n: int, a: float = -1.0, b: float = 1.0) -> np.ndarray:
    """n + 1 Chebyshev nodes (roots of T_{n+1}) mapped to [a, b]."""
    k = np.arange(n + 1)
    x = np.cos((2 * k + 1) * np.pi / (2 * (n + 1)))
    return 0.5 * (a + b) + 0.5 * (b - a) * x


def barycentric_weights(nodes: np.ndarray) -> np.ndarray:
    """Weights w_j = 1 / prod_{k != j} (x_j - x_k) of the barycentric Lagrange formula."""
    diff = nodes[:, None] - nodes[None, :]
    np.fill_diagonal(diff, 1.0)
    return 1.0 / np.prod(diff, axis=1)


def lagrange_interpolate(nodes: np.ndarray, values: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Evaluate the interpolating polynomial at x with the (numerically stable) barycentric formula."""
    w = barycentric_weights(nodes)
    x = np.asarray(x, dtype=float)
    d = x[:, None] - nodes[None, :]
    exact = np.isclose(d, 0.0, atol=1e-14)
    d[exact] = 1.0
    t = w / d
    p = (t @ values) / t.sum(axis=1)
    rows, cols = np.nonzero(exact)
    p[rows] = values[cols]
    return p


def runge(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + 25.0 * x ** 2)
