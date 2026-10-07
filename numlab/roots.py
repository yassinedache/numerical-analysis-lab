"""Root-finding methods for f(x) = 0.

Every method returns a RootResult holding the approximate root and the full list
of iterates, so the convergence behaviour can be studied afterwards.
"""

from dataclasses import dataclass, field
from typing import Callable, List

Function = Callable[[float], float]


@dataclass
class RootResult:
    root: float
    iterates: List[float] = field(default_factory=list)
    converged: bool = False

    @property
    def iterations(self) -> int:
        return len(self.iterates)


def bisection(f: Function, a: float, b: float, tol: float = 1e-12, max_iter: int = 200) -> RootResult:
    """Bisection method. Requires f(a) and f(b) of opposite signs.

    Linear convergence: the error bound is halved at every step, |x_n - r| <= (b - a) / 2^(n+1).
    """
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must have opposite signs")
    iterates = []
    for _ in range(max_iter):
        m = (a + b) / 2
        iterates.append(m)
        fm = f(m)
        if fm == 0 or (b - a) / 2 < tol:
            return RootResult(m, iterates, True)
        if fa * fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return RootResult(iterates[-1], iterates, False)


def fixed_point(g: Function, x0: float, tol: float = 1e-12, max_iter: int = 500) -> RootResult:
    """Fixed-point iteration x_{n+1} = g(x_n).

    Converges linearly with rate |g'(r)| when |g'(r)| < 1 near the fixed point r.
    """
    x = x0
    iterates = []
    for _ in range(max_iter):
        x_new = g(x)
        iterates.append(x_new)
        if abs(x_new - x) < tol:
            return RootResult(x_new, iterates, True)
        x = x_new
    return RootResult(x, iterates, False)


def newton(f: Function, df: Function, x0: float, tol: float = 1e-12, max_iter: int = 100) -> RootResult:
    """Newton's method x_{n+1} = x_n - f(x_n) / f'(x_n).

    Quadratic convergence (order 2) near a simple root.
    """
    x = x0
    iterates = []
    for _ in range(max_iter):
        d = df(x)
        if d == 0:
            raise ZeroDivisionError("f'(x) = 0: Newton's method cannot continue")
        x_new = x - f(x) / d
        iterates.append(x_new)
        if abs(x_new - x) < tol:
            return RootResult(x_new, iterates, True)
        x = x_new
    return RootResult(x, iterates, False)


def secant(f: Function, x0: float, x1: float, tol: float = 1e-12, max_iter: int = 100) -> RootResult:
    """Secant method: Newton's method with the derivative replaced by a finite difference.

    Superlinear convergence of order (1 + sqrt(5)) / 2, about 1.618 (the golden ratio).
    """
    iterates = []
    f0, f1 = f(x0), f(x1)
    for _ in range(max_iter):
        if f1 == f0:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        iterates.append(x2)
        if abs(x2 - x1) < tol:
            return RootResult(x2, iterates, True)
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
    return RootResult(iterates[-1] if iterates else x1, iterates, False)
