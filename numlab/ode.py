"""One-step methods for the initial value problem y' = f(t, y), y(t0) = y0.

Works for scalar equations and for systems (y as a NumPy array).
Global error orders on a fixed interval:
    explicit Euler   O(h)
    Heun (RK2)       O(h^2)
    Runge-Kutta 4    O(h^4)
"""

from typing import Callable, Tuple

import numpy as np

RHS = Callable[[float, np.ndarray], np.ndarray]


def _euler_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    return y + h * f(t, y)


def _heun_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    k1 = f(t, y)
    k2 = f(t + h, y + h * k1)
    return y + h / 2 * (k1 + k2)


def _rk4_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    k1 = f(t, y)
    k2 = f(t + h / 2, y + h / 2 * k1)
    k3 = f(t + h / 2, y + h / 2 * k2)
    k4 = f(t + h, y + h * k3)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


STEPS = {"euler": _euler_step, "heun": _heun_step, "rk4": _rk4_step}


def solve(f: RHS, t0: float, y0, t_end: float, n_steps: int, method: str = "rk4") -> Tuple[np.ndarray, np.ndarray]:
    """Integrate from t0 to t_end with n_steps equal steps.

    Returns (t, Y) where Y[k] is the approximation of y(t[k]).
    """
    if method not in STEPS:
        raise ValueError(f"unknown method {method!r}, choose from {sorted(STEPS)}")
    step = STEPS[method]
    h = (t_end - t0) / n_steps
    t = t0 + h * np.arange(n_steps + 1)
    y = np.atleast_1d(np.asarray(y0, dtype=float))
    Y = np.empty((n_steps + 1, y.size))
    Y[0] = y
    for k in range(n_steps):
        y = step(f, t[k], y, h)
        Y[k + 1] = y
    return t, Y
