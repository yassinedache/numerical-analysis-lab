import math

import numpy as np
import pytest

from numlab import interpolation, monte_carlo


def test_interpolation_passes_through_nodes():
    nodes = interpolation.chebyshev_nodes(8)
    vals = np.sin(nodes)
    assert interpolation.lagrange_interpolate(nodes, vals, nodes) == pytest.approx(vals)


def test_polynomials_are_reproduced_exactly():
    nodes = interpolation.equispaced_nodes(5)
    p = lambda x: 3 * x ** 5 - x ** 2 + 2
    xx = np.linspace(-1, 1, 101)
    assert interpolation.lagrange_interpolate(nodes, p(nodes), xx) == pytest.approx(p(xx), abs=1e-10)


def test_runge_phenomenon():
    xx = np.linspace(-1, 1, 2001)
    fx = interpolation.runge(xx)

    def max_err(nodes):
        return np.max(np.abs(interpolation.lagrange_interpolate(nodes, interpolation.runge(nodes), xx) - fx))

    assert max_err(interpolation.equispaced_nodes(20)) > max_err(interpolation.equispaced_nodes(10))
    assert max_err(interpolation.chebyshev_nodes(20)) < max_err(interpolation.chebyshev_nodes(10))


def test_monte_carlo_pi_and_confidence_interval():
    assert monte_carlo.estimate_pi(200_000, np.random.default_rng(42)) == pytest.approx(math.pi, abs=0.02)
    est, (lo, hi) = monte_carlo.mc_integral(np.sin, 0.0, math.pi, 200_000, np.random.default_rng(7))
    assert lo < 2.0 < hi
