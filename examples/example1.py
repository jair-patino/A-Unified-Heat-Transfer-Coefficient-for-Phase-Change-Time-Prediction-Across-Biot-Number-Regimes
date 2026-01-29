"""
Example 1: Analytical model for a planar wall.
"""

from src.analytical_model import geometric_factor, global_U, phase_change_time

# Parameters
geom = 'planar'
thickness = 0.01  # m
A = 1.0  # m2
V = A * thickness
m = V * 1000  # kg, assuming density 1000 kg/m3
c = 4186  # J/kgK
L = 334000  # J/kg
k = 0.6  # W/mK
heff = 100  # W/m2K
Ti = 20  # °C
Tf_dagger = 0  # °C
Tinf = -10  # °C

# Geometric factor and characteristic length
Phi, Lc_expr = geometric_factor(geom)
if geom == 'planar':
    Lc = thickness / 2
# For cylinder and sphere, the radius would be needed

# Global coefficient
U = global_U(heff, k, Lc, Phi)

# Phase change time
t_total = phase_change_time(m, c, L, A, Ti, Tf_dagger, Tinf, U)
print(f"Total phase change time: {t_total} s")
