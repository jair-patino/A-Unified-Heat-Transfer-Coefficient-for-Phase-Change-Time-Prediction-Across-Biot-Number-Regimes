"""
Analytical model to predict the phase change time.
"""

import numpy as np

def geometric_factor(geom):
    """
    Returns the geometric factor Φ and the characteristic length Lc.

    geom : str
        Geometry: 'planar', 'cylinder', 'sphere'

    Returns:
    Phi : float
    Lc : float
    """
    factors = {
        'planar': (1.0, 'L/2'),   # L is the thickness
        'cylinder': (0.5, 'R/2'), # R is the radius
        'sphere': (1/3, 'R/3')
    }
    return factors[geom]

def global_U(heff, k, Lc, Phi):
    """
    Global heat transfer coefficient U.

    heff : float
        Effective external coefficient
    k : float
        Thermal conductivity
    Lc : float
        Characteristic length
    Phi : float
        Geometric factor

    Returns:
    U : float
    """
    Bi = heff * Lc / k
    return heff / (1 + Phi * Bi)

def phase_change_time(m, c, L, A, Ti, Tf_dagger, Tinf, U):
    """
    Total phase change time (Equation 9).

    m : float
        Mass
    c : float
        Specific heat
    L : float
        Latent heat
    A : float
        Surface area
    Ti : float
        Initial temperature
    Tf_dagger : float
        Effective phase change temperature
    Tinf : float
        Ambient temperature
    U : float
        Global coefficient

    Returns:
    t_total : float
    """
    t_sensible = (m*c/(U*A)) * np.log(abs(Ti - Tinf)/abs(Tf_dagger - Tinf))
    t_latent = (m*L/(U*A)) / abs(Tinf - Tf_dagger)
    return t_sensible + t_latent
