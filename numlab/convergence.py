"""Tools to measure convergence orders from experimental data."""

from typing import Sequence

import numpy as np


def observed_order(h: Sequence[float], errors: Sequence[float]) -> float:
    """Slope p of the least-squares line through (log h, log error), so that error ~ C h^p."""
    h = np.asarray(h, dtype=float)
    e = np.asarray(errors, dtype=float)
    keep = e > 0
    slope, _ = np.polyfit(np.log(h[keep]), np.log(e[keep]), 1)
    return float(slope)


def linear_rate(errors: Sequence[float]) -> float:
    """Average contraction factor C of a linearly convergent method (e_{n+1} ~ C e_n).

    Fits log(e_n) against n; works even when the errors are not monotone (bisection).
    """
    e = np.asarray([x for x in errors if x > 1e-15], dtype=float)
    if e.size < 3:
        raise ValueError("need at least three non-zero errors")
    slope, _ = np.polyfit(np.arange(e.size), np.log(e), 1)
    return float(np.exp(slope))


def iteration_order(errors: Sequence[float]) -> float:
    """Order q of an iterative method from consecutive errors, using e_{n+1} ~ C e_n^q.

    Uses the last three meaningful errors: q ~ log(e_{n+1} / e_n) / log(e_n / e_{n-1}).
    """
    e = [x for x in errors if x > 1e-15]
    if len(e) < 3:
        raise ValueError("need at least three non-zero errors")
    e0, e1, e2 = e[-3], e[-2], e[-1]
    return float(np.log(e2 / e1) / np.log(e1 / e0))
