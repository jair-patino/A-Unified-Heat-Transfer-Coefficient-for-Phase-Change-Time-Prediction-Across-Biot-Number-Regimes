"""
Numerical solution of the Stefan problem using the enthalpy method with finite differences.
Reference: Voller and Swaminathan (1997)
"""

import numpy as np

def solve_stefan_1d(geom, Bi, Ste, Theta, N=200, dt=0.01, max_time=1000):
    """
    Solves the one-dimensional Stefan problem.

    Parameters:
    geom : str
        Geometry: 'planar', 'cylinder', 'sphere'
    Bi : float
        Biot number
    Ste : float
        Stefan number
    Theta : float
        Dimensionless temperature difference
    N : int
        Number of nodes
    dt : float
        Time step
    max_time : float
        Maximum simulation time

    Returns:
    time : float
        Total phase change time
    """
    # Implementation of the enthalpy scheme
    # ... (detailed code in the article and supplement)
    pass
