"""
Numerical solver for 1D Stefan problem using enthalpy method with finite differences.
Reference implementation for validation of analytical model.
Author: Jair Patino B.
Date: January 2026
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time

class EnthalpyStefanSolver:
    """
    Solves 1D phase change problems using enthalpy formulation.
    Geometry options: 'planar', 'cylindrical', 'spherical'
    """
    
    def __init__(self, geometry='planar', N=200, max_iter=100, tol=1e-6):
        """
        Initialize solver.
        
        Parameters:
        -----------
        geometry : str
            'planar', 'cylindrical', or 'spherical'
        N : int
            Number of grid points
        max_iter : int
            Maximum Newton-Raphson iterations
        tol : float
            Convergence tolerance (K)
        """
        self.geometry = geometry
        self.N = N
        self.max_iter = max_iter
        self.tol = tol
        
        # Geometry factors for coordinate transformation
        self.geometry_factors = {
            'planar': {'A': lambda r: 1.0, 'V': lambda r: 1.0},
            'cylindrical': {'A': lambda r: 2*np.pi*r, 'V': lambda r: np.pi*r**2},
            'spherical': {'A': lambda r: 4*np.pi*r**2, 'V': lambda r: 4/3*np.pi*r**3}
        }
        
    def setup_grid(self, L, N=None):
        """
        Create computational grid.
        
        Parameters:
        -----------
        L : float
            Characteristic length (m)
        N : int, optional
            Number of grid points (overrides initialization)
        
        Returns:
        --------
        dict with grid parameters
        """
        if N is not None:
            self.N = N
            
        dr = L / (self.N - 1)
        r = np.linspace(0, L, self.N)
        
        # Area and volume factors for control volumes
        A = self.geometry_factors[self.geometry]['A'](r)
        V = self.geometry_factors[self.geometry]['V'](r)
        
        return {
            'r': r,
            'dr': dr,
            'A': A,
            'V': V,
            'N': self.N
        }
    
    def enthalpy_function(self, T, T_pc, L, c_l, c_s, delta_T_mush=0.5):
        """
        Calculate enthalpy and apparent heat capacity.
        
        Parameters:
        -----------
        T : float or array
            Temperature (K or °C)
        T_pc : float
            Phase change temperature (K or °C)
        L : float
            Latent heat (J/kg)
        c_l, c_s : float
            Specific heat of liquid and solid (J/kg·K)
        delta_T_mush : float
            Mushy zone temperature interval (K)
        
        Returns:
        --------
        H : enthalpy (J/kg)
        c_app : apparent heat capacity (J/kg·K)
        f_l : liquid fraction
        """
        T = np.asarray(T)
        H = np.zeros_like(T)
        c_app = np.zeros_like(T)
        f_l = np.zeros_like(T)
        
        # Temperature bounds for mushy zone
        T_solid = T_pc - delta_T_mush/2
        T_liquid = T_pc + delta_T_mush/2
        
        for i, Ti in enumerate(T):
            if Ti <= T_solid:
                # Solid
                H[i] = c_s * (Ti - T_pc)
                c_app[i] = c_s
                f_l[i] = 0.0
            elif Ti >= T_liquid:
                # Liquid
                H[i] = c_l * (Ti - T_pc) + L
                c_app[i] = c_l
                f_l[i] = 1.0
            else:
                # Mushy zone
                f_l[i] = (Ti - T_solid) / delta_T_mush
                H[i] = (c_s * (T_solid - T_pc) + 
                       f_l[i] * L + 
                       (f_l[i]*c_l + (1-f_l[i])*c_s) * (Ti - T_solid))
                c_app[i] = (H[i] - c_s*(T_solid - T_pc)) / (Ti - T_solid)
        
        return H, c_app, f_l
    
    def solve_transient(self, grid, properties, boundary_conditions, 
                       dt=1.0, t_final=3600, save_interval=100):
        """
        Solve transient phase change problem.
        
        Parameters:
        -----------
        grid : dict
            Grid parameters from setup_grid
        properties : dict
            Material properties: rho, k, c_l, c_s, L, T_pc
        boundary_conditions : dict
            h_eff, T_inf, T_initial
        dt : float
            Time step (s)
        t_final : float
            Final time (s)
        save_interval : int
            Save solution every n time steps
        
        Returns:
        --------
        dict with solution history
        """
        # Unpack inputs
        r = grid['r']
        dr = grid['dr']
        N = grid['N']
        A = grid['A']
        
        rho = properties['rho']
        k = properties['k']
        c_l = properties['c_l']
        c_s = properties['c_s']
        L = properties['L']
        T_pc = properties['T_pc']
        
        h_eff = boundary_conditions['h_eff']
        T_inf = boundary_conditions['T_inf']
        T_initial = boundary_conditions['T_initial']
        
        # Initialize
        T = np.ones(N) * T_initial
        time_steps = int(t_final / dt)
        
        # Storage for results
        solution = {
            'time': [],
            'temperature': [],
            'liquid_fraction': [],
            'phase_change_time': None
        }
        
        # Time stepping
        print(f"Starting transient solution: {time_steps} time steps")
        start_time = time.time()
        
        for n in range(time_steps):
            t = n * dt
            
            # Newton-Raphson iteration for implicit scheme
            T_old = T.copy()
            
            for iter in range(self.max_iter):
                # Calculate enthalpy and derivatives
                H, c_app, f_l = self.enthalpy_function(T, T_pc, L, c_l, c_s)
                
                # Build linear system A * T = b
                A_mat = sparse.lil_matrix((N, N))
                b = np.zeros(N)
                
                # Interior points
                for i in range(1, N-1):
                    # Geometric factors
                    if self.geometry == 'planar':
                        r_factor = 1.0
                    elif self.geometry == 'cylindrical':
                        r_factor = r[i]
                    else:  # spherical
                        r_factor = r[i]**2
                    
                    # East and west coefficients
                    k_e = 0.5 * (k + k)  # Constant k assumed
                    k_w = 0.5 * (k + k)
                    
                    a_E = k_e * A[i+1] / dr
                    a_W = k_w * A[i] / dr
                    a_P = rho * c_app[i] * (r_factor * dr) / dt
                    
                    A_mat[i, i+1] = -a_E
                    A_mat[i, i-1] = -a_W
                    A_mat[i, i] = a_P + a_E + a_W
                    b[i] = a_P * T_old[i]
                
                # Boundary conditions
                # Left boundary (center) - symmetry
                A_mat[0, 0] = 1.0
                A_mat[0, 1] = -1.0
                b[0] = 0.0
                
                # Right boundary (surface) - convective
                A_mat[N-1, N-1] = k/dr + h_eff
                A_mat[N-1, N-2] = -k/dr
                b[N-1] = h_eff * T_inf
                
                # Solve linear system
                T_new = spsolve(A_mat.tocsr(), b)
                
                # Check convergence
                error = np.max(np.abs(T_new - T))
                T = T_new
                
                if error < self.tol:
                    break
                    
                if iter == self.max_iter - 1:
                    print(f"Warning: No convergence at time step {n}, error = {error:.2e}")
            
            # Save solution at intervals
            if n % save_interval == 0 or n == time_steps - 1:
                solution['time'].append(t)
                solution['temperature'].append(T.copy())
                solution['liquid_fraction'].append(f_l.copy())
                
                # Check if phase change is complete
                if f_l.mean() < 0.01 and solution['phase_change_time'] is None:
                    solution['phase_change_time'] = t
        
        elapsed = time.time() - start_time
        print(f"Solution completed in {elapsed:.2f} seconds")
        
        return solution
    
    def compute_phase_change_time(self, properties, geometry_params, 
                                 boundary_conditions, dt=0.1):
        """
        Wrapper function to compute total phase change time.
        
        Parameters:
        -----------
        properties : dict
            Material properties
        geometry_params : dict
            'L' for planar thickness, 'R' for cylinder/sphere radius
        boundary_conditions : dict
            h_eff, T_inf, T_initial
        
        Returns:
        --------
        float : total phase change time (s)
        """
        # Set up grid
        if self.geometry == 'planar':
            L = geometry_params['L']
        else:
            L = geometry_params['R']
            
        grid = self.setup_grid(L)
        
        # Solve transient problem
        # Use adaptive time stepping for efficiency
        t_estimate = self.estimate_phase_change_time(properties, geometry_params, 
                                                    boundary_conditions)
        dt = min(dt, t_estimate / 1000)  # Ensure sufficient time resolution
        t_final = t_estimate * 1.5  # Add safety factor
        
        solution = self.solve_transient(grid, properties, boundary_conditions,
                                       dt=dt, t_final=t_final)
        
        return solution['phase_change_time']
    
    def estimate_phase_change_time(self, properties, geometry_params, 
                                  boundary_conditions):
        """
        Quick estimate of phase change time for setting simulation parameters.
        Based on lumped capacitance approximation.
        """
        if self.geometry == 'planar':
            V = geometry_params['L'] * 1.0 * 1.0  # 1 m^2 cross-section
            A = 2.0  # Both sides
        elif self.geometry == 'cylindrical':
            R = geometry_params['R']
            V = np.pi * R**2 * 1.0  # 1 m length
            A = 2 * np.pi * R * 1.0
        else:  # spherical
            R = geometry_params['R']
            V = 4/3 * np.pi * R**3
            A = 4 * np.pi * R**2
        
        h_eff = boundary_conditions['h_eff']
        T_inf = boundary_conditions['T_inf']
        T_initial = boundary_conditions['T_initial']
        
        # Simple lumped capacitance estimate
        m = properties['rho'] * V
        t_sensible = (m * properties['c_l'] / (h_eff * A) * 
                     np.log(abs(T_initial - T_inf) / abs(properties['T_pc'] - T_inf)))
        t_latent = m * properties['L'] / (h_eff * A * abs(T_inf - properties['T_pc']))
        
        return t_sensible + t_latent


def run_validation_case(material='water', geometry='planar', Bi=1.0):
    """
    Example function to run a validation case.
    """
    # Material properties (example for water)
    properties = {
        'water': {
            'rho': 1000,      # kg/m^3
            'k': 0.6,         # W/m·K
            'c_l': 4186,      # J/kg·K
            'c_s': 2100,      # J/kg·K
            'L': 334000,      # J/kg
            'T_pc': 0.0       # °C
        },
        'paraffin': {
            'rho': 900,
            'k': 0.2,
            'c_l': 2200,
            'c_s': 1800,
            'L': 200000,
            'T_pc': 50.0
        }
    }
    
    # Geometry
    L_c = 0.01  # Characteristic length 1 cm
    if geometry == 'planar':
        geometry_params = {'L': L_c * 2}  # Full thickness
    else:
        geometry_params = {'R': L_c * (3 if geometry == 'spherical' else 2)}
    
    # Calculate h_eff from Biot number
    Biot = Bi
    k = properties[material]['k']
    h_eff = Biot * k / L_c
    
    # Boundary conditions
    boundary_conditions = {
        'h_eff': h_eff,
        'T_inf': -20.0,  # °C
        'T_initial': 20.0  # °C
    }
    
    # Create solver
    solver = EnthalpyStefanSolver(geometry=geometry, N=200)
    
    # Solve
    t_total = solver.compute_phase_change_time(
        properties[material], geometry_params, boundary_conditions
    )
    
    print(f"Material: {material}, Geometry: {geometry}, Bi: {Bi}")
    print(f"Total phase change time: {t_total:.1f} s")
    
    return t_total


if __name__ == "__main__":
    # Example usage
    t_water = run_validation_case('water', 'planar', Bi=0.5)
    t_paraffin = run_validation_case('paraffin', 'spherical', Bi=1.0)
