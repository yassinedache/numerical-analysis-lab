import math

import numpy as np
import pytest

from numlab import convergence, quadrature

NS = [4, 8, 16, 32, 64]
HS = [math.pi / n for n in NS]


def errors(rule):
    return [abs(rule(np.sin, 0, math.pi, n) - 2.0) for n in NS]


@pytest.mark.parametrize("rule, order", [
    (quadrature.midpoint, 2),
    (quadrature.trapezoid, 2),
    (quadrature.simpson, 4),
])
def test_composite_rules_have_their_theoretical_order(rule, order):
    assert convergence.observed_order(HS, errors(rule)) == pytest.approx(order, abs=0.15)


def test_simpson_needs_even_n():
    with pytest.raises(ValueError):
        quadrature.simpson(np.sin, 0, 1, 3)


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_gauss_legendre_is_exact_up_to_degree_2n_minus_1(n):
    rng = np.random.default_rng(n)
    p = np.polynomial.Polynomial(rng.normal(size=2 * n))  # polynomial of degree 2n - 1
    exact = p.integ()(2.0) - p.integ()(-1.0)
    assert quadrature.gauss_legendre(p, -1.0, 2.0, n) == pytest.approx(exact, rel=1e-12, abs=1e-12)
