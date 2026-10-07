"""numlab: a small numerical analysis lab.

Root finding, numerical integration, ODE solvers, polynomial interpolation and
Monte Carlo methods, each with experiments that measure the theoretical
convergence orders.
"""

from . import convergence, interpolation, monte_carlo, ode, quadrature, roots

__all__ = ["convergence", "interpolation", "monte_carlo", "ode", "quadrature", "roots"]
__version__ = "1.0.0"
