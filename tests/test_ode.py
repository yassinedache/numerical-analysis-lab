import math

import numpy as np
import pytest

from numlab import convergence, ode

f = lambda t, y: -2 * t * y
EXACT = math.exp(-4.0)
NS = [20, 40, 80, 160, 320]
HS = [2.0 / n for n in NS]


@pytest.mark.parametrize("method, order", [("euler", 1), ("heun", 2), ("rk4", 4)])
def test_global_error_order(method, order):
    err = [abs(ode.solve(f, 0.0, 1.0, 2.0, n, method)[1][-1, 0] - EXACT) for n in NS]
    assert convergence.observed_order(HS, err) == pytest.approx(order, abs=0.2)


def test_systems_are_supported():
    # harmonic oscillator x'' = -x as a first-order system: after one period RK4 returns to the start
    t, Y = ode.solve(lambda t, y: np.array([y[1], -y[0]]), 0.0, [1.0, 0.0], 2 * math.pi, 400, "rk4")
    assert Y[-1] == pytest.approx([1.0, 0.0], abs=1e-6)


def test_unknown_method():
    with pytest.raises(ValueError):
        ode.solve(f, 0, 1, 1, 10, "magic")
