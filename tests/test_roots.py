import math

import pytest

from numlab import convergence, roots

f = lambda x: math.cos(x) - x
df = lambda x: -math.sin(x) - 1
R = 0.7390851332151607


@pytest.mark.parametrize("res", [
    roots.bisection(f, 0, 1),
    roots.fixed_point(math.cos, 1.0),
    roots.newton(f, df, 1.0),
    roots.secant(f, 0, 1),
])
def test_every_method_finds_the_root(res):
    assert res.converged
    assert abs(res.root - R) < 1e-10


def test_bisection_rejects_bad_bracket():
    with pytest.raises(ValueError):
        roots.bisection(f, 1, 2)


def test_newton_is_quadratic():
    err = [abs(x - R) for x in roots.newton(f, df, 1.0).iterates if abs(x - R) > 1e-13]
    assert convergence.iteration_order(err) == pytest.approx(2.0, abs=0.15)


def test_secant_order_is_golden_ratio():
    err = [abs(x - R) for x in roots.secant(f, 0, 1).iterates if abs(x - R) > 1e-13]
    assert convergence.iteration_order(err) == pytest.approx((1 + 5 ** 0.5) / 2, abs=0.15)


def test_bisection_halves_the_error_on_average():
    err = [abs(x - R) for x in roots.bisection(f, 0, 1).iterates]
    assert convergence.linear_rate(err) == pytest.approx(0.5, abs=0.05)
