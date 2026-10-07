"""Numerical integration (quadrature) of a function on [a, b].

Composite rules on n sub-intervals of width h = (b - a) / n, plus Gauss-Legendre.
Theoretical error orders for smooth functions:
    midpoint   O(h^2)
    trapezoid  O(h^2)
    Simpson    O(h^4)
    Gauss-Legendre with n points: exact for polynomials of degree <= 2n - 1,
    error decreasing faster than any power of 1/n for analytic functions.
"""

from typing import Callable

import numpy as np

Function = Callable[[np.ndarray], np.ndarray]


def midpoint(f: Function, a: float, b: float, n: int) -> float:
    h = (b - a) / n
    x = a + h * (np.arange(n) + 0.5)
    return float(h * np.sum(f(x)))


def trapezoid(f: Function, a: float, b: float, n: int) -> float:
    x = np.linspace(a, b, n + 1)
    y = f(x)
    h = (b - a) / n
    return float(h * (0.5 * y[0] + np.sum(y[1:-1]) + 0.5 * y[-1]))


def simpson(f: Function, a: float, b: float, n: int) -> float:
    """Composite Simpson rule; n must be even."""
    if n % 2:
        raise ValueError("Simpson's rule needs an even number of sub-intervals")
    x = np.linspace(a, b, n + 1)
    y = f(x)
    h = (b - a) / n
    return float(h / 3 * (y[0] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-1:2]) + y[-1]))


def gauss_legendre(f: Function, a: float, b: float, n: int) -> float:
    """Gauss-Legendre rule with n nodes, mapped from [-1, 1] to [a, b]."""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    x = 0.5 * (b - a) * nodes + 0.5 * (a + b)
    return float(0.5 * (b - a) * np.sum(weights * f(x)))
