#!/usr/bin/env python3
"""
Enthalpy Method Solver for 1D Stefan Problems

This module implements an enthalpy-based finite difference solver for 
the one-dimensional Stefan problem with convective boundary conditions.

Author: Jair Patiño B.
Date: January 2026
License: MIT
"""

import numpy as np
from typing import Tuple, Dict, Optional
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time
from dataclasses import dataclass


@dataclass
class MaterialProperties:
    """Container for material properties"""
    k: float          # Thermal conductivity [W/(m·K)]
    c: float          # Specific heat [J/(kg·K)]
    L: float          # Latent heat [J/kg]
    rho: float        # Density [kg/m³]
    T_melt: float     # Melting temperature [K]
    name: str = "Unknown"


@dataclass
class Geometry:
    """Container for geometry parameters"""
    shape: str        # 'planar', 'cylindrical', 'spherical'
    L: float          # Characteristic length [m]
    nodes: int = 200  # Number of grid points
    A: float = 1.0    # Surface area [m²]
    V: float = 1.0    # Volume [m³]


class EnthalpySolver:
    """
    1D Enthalpy Method Solver for Stefan Problems
    
    Solves the energy equation in enthalpy form:
        ρ ∂H/∂t = ∇·(k ∇T)
    where H = ∫ ρc dT + fL is the enthalpy.
    
    Reference: Voller & Prakash (1987)
    """
    
    def __init__(self, material: MaterialProperties, geometry: Geometry):
        """
        Initialize the enthalpy solver.
        
        Parameters
        ----------
        material : MaterialProperties
            Material thermophysical properties
        geometry : Geometry
            Geometry and discretization parameters
        """
        self.material = material
        self.geometry = geometry
        
        # Derived properties
        self.alpha = material.k / (material.rho * material.c)  # Thermal diffusivity
        self.Ste = None  # Stefan number (will be set from BCs)
        
        # Grid setup
        self._setup_grid()
        
        # Initialize arrays
        self.T = np.zeros(geometry.nodes)  # Temperature [K]
        self.H = np.zeros(geometry.nodes)  # Enthalpy [J/m³]
        self.f = np.zeros(geometry.nodes)  # Liquid fraction [0-1]
        
        # Solver parameters
        self.max_iter = 100
        self.tolerance = 1e-6
        self.omega = 0.8  # Under-relaxation factor
        
        # Results storage
        self.time_history = []
        self.temp_history = []
        self.interface_position = []
        
    def _setup_grid(self):
        """Setup the computational grid"""
        if self.geometry.shape == 'planar':
            # Linear grid for planar geometry
            self.x = np.linspace(0, self.geometry.L, self.geometry.nodes)
            self.dx = self.x[1] - self.x[0]
            self.vol = self.dx * self.geometry.A  # Volume per node
            
        elif self.geometry.shape == 'cylindrical':
            # Radial grid for cylinder
            self.x = np.linspace(0, self.geometry.L, self.geometry.nodes)
            self.dx = self.x[1] - self.x[0]
            r = self.x
            # Volume elements for cylindrical shells
            self.vol = np.pi * ((r + self.dx/2)**2 - (r - self.dx/2)**2) * 1.0  # per unit length
            
        elif self.geometry.shape == 'spherical':
            # Radial grid for sphere
            self.x = np.linspace(0, self.geometry.L, self.geometry.nodes)
            self.dx = self.x[1] - self.x[0]
            r = self.x
            # Volume elements for spherical shells
            self.vol = (4/3) * np.pi * ((r + self.dx/2)**3 - (r - self.dx/2)**3)
            
        else:
            raise ValueError(f"Unsupported geometry: {self.geometry.shape}")
    
    def _enthalpy_from_temp(self, T: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert temperature to enthalpy and liquid fraction.
        
        Uses an apparent heat capacity method with a small mushy zone.
        
        Parameters
        ----------
        T : np.ndarray
            Temperature array [K]
            
        Returns
        -------
        H : np.ndarray
            Enthalpy [J/m³]
        f : np.ndarray
            Liquid fraction [0-1]
        """
        delta_T = 0.5  # Mushy zone width [K]
        T_m = self.material.T_melt
        
        # Initialize arrays
        H = np.zeros_like(T)
        f = np.zeros_like(T)
        
        # Fully solid
        solid_mask = T <= (T_m - delta_T/2)
        H[solid_mask] = self.material.rho * self.material.c * (T[solid_mask] - T_m)
        
        # Mushy zone
        mushy_mask = np.abs(T - T_m) <= delta_T/2
        if np.any(mushy_mask):
            # Linear interpolation in mushy zone
            f[mushy_mask] = 0.5 + (T[mushy_mask] - T_m) / delta_T
            H[mushy_mask] = (self.material.rho * self.material.c * 
                            (T[mushy_mask] - T_m) + 
                            f[mushy_mask] * self.material.rho * self.material.L)
        
        # Fully liquid
        liquid_mask = T >= (T_m + delta_T/2)
        f[liquid_mask] = 1.0
        H[liquid_mask] = (self.material.rho * self.material.c * 
                         (T[liquid_mask] - T_m) + 
                         self.material.rho * self.material.L)
        
        return H, f
    
    def _temp_from_enthalpy(self, H: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert enthalpy to temperature and liquid fraction.
        
        Inverse of _enthalpy_from_temp.
        
        Parameters
        ----------
        H : np.ndarray
            Enthalpy [J/m³]
            
        Returns
        -------
        T : np.ndarray
            Temperature [K]
        f : np.ndarray
            Liquid fraction [0-1]
        """
        delta_T = 0.5  # Mushy zone width [K]
        T_m = self.material.T_melt
        
        # Initialize arrays
        T = np.zeros_like(H)
        f = np.zeros_like(H)
        
        # Reference enthalpies
        H_solid = -self.material.rho * self.material.c * delta_T/2
        H_liquid = self.material.rho * self.material.c * delta_T/2 + self.material.rho * self.material.L
        
        # Fully solid
        solid_mask = H <= H_solid
        T[solid_mask] = T_m + H[solid_mask] / (self.material.rho * self.material.c)
        f[solid_mask] = 0.0
        
        # Mushy zone
        mushy_mask = (H > H_solid) & (H < H_liquid)
        if np.any(mushy_mask):
            f[mushy_mask] = (H[mushy_mask] - H_solid) / (H_liquid - H_solid)
            T[mushy_mask] = T_m + (f[mushy_mask] - 0.5) * delta_T
        
        # Fully liquid
        liquid_mask = H >= H_liquid
        f[liquid_mask] = 1.0
        T[liquid_mask] = T_m + (H[liquid_mask] - self.material.rho * self.material.L) / \
                         (self.material.rho * self.material.c)
        
        return T, f
    
    def _build_matrix(self, dt: float, h_eff: float) -> sparse.csr_matrix:
        """
        Build the sparse coefficient matrix for the implicit scheme.
        
        Parameters
        ----------
        dt : float
            Time step [s]
        h_eff : float
            Effective heat transfer coefficient [W/(m²·K)]
            
        Returns
        -------
        A : sparse.csr_matrix
            Coefficient matrix (n x n)
        """
        n = self.geometry.nodes
        dx = self.dx
        k = self.material.k
        rho = self.material.rho
        
        # Pre-calculate geometric factors
        if self.geometry.shape == 'planar':
            # Constant cross-section
            factor = k / dx**2
            
        elif self.geometry.shape == 'cylindrical':
            # Cylindrical coordinates (1/r ∂/∂r (r ∂T/∂r))
            r = self.x
            factor = k / (dx**2)
            
        elif self.geometry.shape == 'spherical':
            # Spherical coordinates (1/r² ∂/∂r (r² ∂T/∂r))
            r = self.x
            factor = k / (dx**2)
        
        # Initialize sparse matrix
        data = []
        rows = []
        cols = []
        
        # Internal nodes
        for i in range(1, n-1):
            # Central coefficient
            if self.geometry.shape == 'planar':
                a_p = rho * self.vol / dt + 2 * factor
            elif self.geometry.shape == 'cylindrical':
                a_p = rho * self.vol[i] / dt + factor * (1 + 0.5*dx/r[i]) + factor * (1 - 0.5*dx/r[i])
            elif self.geometry.shape == 'spherical':
                a_p = rho * self.vol[i] / dt + factor * (1 + dx/r[i]) + factor * (1 - dx/r[i])
            
            data.append(a_p)
            rows.append(i)
            cols.append(i)
            
            # West coefficient
            data.append(-factor)
            rows.append(i)
            cols.append(i-1)
            
            # East coefficient
            data.append(-factor)
            rows.append(i)
            cols.append(i+1)
        
        # Boundary conditions
        # Left boundary (center) - symmetry or fixed
        data.append(1.0)  # Dirichlet or symmetry
        rows.append(0)
        cols.append(0)
        
        # Right boundary (surface) - convective
        if self.geometry.shape == 'planar':
            a_surf = rho * self.vol / dt + factor + h_eff * dx / k
        else:
            a_surf = rho * self.vol[-1] / dt + factor + h_eff * dx / k
        
        data.append(a_surf)
        rows.append(n-1)
        cols.append(n-1)
        
        data.append(-factor)
        rows.append(n-1)
        cols.append(n-2)
        
        # Create sparse matrix
        A = sparse.csr_matrix((data, (rows, cols)), shape=(n, n))
        
        return A
    
    def solve(self, T_initial: np.ndarray, T_inf: float, h_eff: float, 
              t_final: float, dt: float = None, adaptive: bool = True) -> Dict:
        """
        Solve the phase change problem.
        
        Parameters
        ----------
        T_initial : np.ndarray
            Initial temperature distribution [K]
        T_inf : float
            Ambient temperature [K]
        h_eff : float
            Effective heat transfer coefficient [W/(m²·K)]
        t_final : float
            Final simulation time [s]
        dt : float, optional
            Time step [s]. If None, computed automatically.
        adaptive : bool
            Whether to use adaptive time stepping.
            
        Returns
        -------
        results : Dict
            Dictionary containing simulation results
        """
        print(f"Starting enthalpy solver for {self.material.name}...")
        print(f"  Geometry: {self.geometry.shape}, Nodes: {self.geometry.nodes}")
        print(f"  T_inf: {T_inf:.1f} K, h_eff: {h_eff:.1f} W/(m²·K)")
        
        # Initialize
        self.T = T_initial.copy()
        self.H, self.f = self._enthalpy_from_temp(self.T)
        
        # Auto time step if not provided
        if dt is None:
            dt = 0.1 * self.dx**2 / self.alpha
            dt = min(dt, 0.1)  # Cap at 0.1s
        
        # Pre-compute matrix (constant for fixed dt)
        A = self._build_matrix(dt, h_eff)
        
        # Time loop
        t = 0.0
        step = 0
        start_time = time.time()
        
        # Initialize results storage
        self.time_history = [t]
        self.temp_history = [self.T.copy()]
        self.interface_position = [self._get_interface_position()]
        
        # Main time loop
        while t < t_final:
            # Adjust dt for final step
            if t + dt > t_final:
                dt = t_final - t
            
            # Source term from previous enthalpy
            b = self.H * self.material.rho / dt
            
            # Add boundary condition at surface
            b[-1] += h_eff * T_inf
            
            # Solve for new enthalpy
            try:
                H_new = spsolve(A, b)
            except:
                # If direct solver fails, use iterative
                H_new, info = sparse.linalg.gmres(A, b, tol=1e-10, maxiter=1000)
                if info != 0:
                    print(f"Warning: Iterative solver did not converge at step {step}")
            
            # Update temperature and liquid fraction
            T_new, f_new = self._temp_from_enthalpy(H_new)
            
            # Under-relaxation for stability
            self.T = self.omega * T_new + (1 - self.omega) * self.T
            self.H = self.omega * H_new + (1 - self.omega) * self.H
            self.f = self.omega * f_new + (1 - self.omega) * self.f
            
            # Update time
            t += dt
            step += 1
            
            # Store results periodically
            if step % 10 == 0 or t >= t_final:
                self.time_history.append(t)
                self.temp_history.append(self.T.copy())
                self.interface_position.append(self._get_interface_position())
            
            # Adaptive time stepping
            if adaptive:
                dT_max = np.max(np.abs(T_new - self.T))
                if dT_max > 10:  # If temperature change too large
                    dt *= 0.8
                    A = self._build_matrix(dt, h_eff)  # Rebuild matrix
                elif dT_max < 1:  # If temperature change too small
                    dt *= 1.2
            
            # Progress indicator
            if step % 100 == 0:
                completion = 100 * t / t_final
                print(f"  Progress: {completion:.1f}%, Time: {t:.1f}s, "
                      f"Interface: {self.interface_position[-1]:.4f}m")
        
        # Compute phase change time
        t_pc = self._compute_phase_change_time()
        
        # Prepare results
        results = {
            'time': np.array(self.time_history),
            'temperature': np.array(self.temp_history),
            'interface_position': np.array(self.interface_position),
            'phase_change_time': t_pc,
            'liquid_fraction': self.f,
            'material': self.material.name,
            'geometry': self.geometry.shape,
            'computation_time': time.time() - start_time,
            'time_steps': step
        }
        
        print(f"Simulation completed in {results['computation_time']:.2f} seconds")
        print(f"Phase change time: {t_pc:.2f} seconds")
        
        return results
    
    def _get_interface_position(self) -> float:
        """Estimate solid-liquid interface position"""
        # Find where liquid fraction crosses 0.5
        if np.all(self.f >= 0.5):
            return 0.0  # All liquid
        elif np.all(self.f <= 0.5):
            return self.geometry.L  # All solid
        
        # Linear interpolation
        idx = np.argmax(self.f < 0.5)
        if idx == 0:
            return 0.0
        
        x1 = self.x[idx-1]
        x2 = self.x[idx]
        f1 = self.f[idx-1]
        f2 = self.f[idx]
        
        if f1 == f2:
            return x1
        
        interface = x1 + (0.5 - f1) * (x2 - x1) / (f2 - f1)
        return interface
    
    def _compute_phase_change_time(self, threshold: float = 0.95) -> float:
        """
        Compute total phase change time.
        
        Parameters
        ----------
        threshold : float
            Liquid fraction threshold for completion (0-1)
            
        Returns
        -------
        t_pc : float
            Phase change time [s]
        """
        # Find when center temperature reaches melting point
        T_center = [T[len(T)//2] for T in self.temp_history]
        times = self.time_history
        
        # Find when temperature crosses melting point
        for i in range(1, len(times)):
            if T_center[i] <= self.material.T_melt < T_center[i-1]:
                # Linear interpolation for more accuracy
                t1 = times[i-1]
                t2 = times[i]
                T1 = T_center[i-1]
                T2 = T_center[i]
                return t1 + (self.material.T_melt - T1) * (t2 - t1) / (T2 - T1)
        
        # If not found, return final time
        return times[-1]
    
    def plot_results(self, results: Dict, save_path: str = None):
        """
        Plot simulation results.
        
        Parameters
        ----------
        results : Dict
            Simulation results from solve() method
        save_path : str, optional
            Path to save the figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Temperature evolution
        ax = axes[0, 0]
        for i in range(0, len(results['temperature']), len(results['temperature'])//5):
            ax.plot(self.x, results['temperature'][i], 
                   label=f"t={results['time'][i]:.1f}s")
        ax.axhline(self.material.T_melt, color='k', linestyle='--', alpha=0.5)
        ax.set_xlabel('Position [m]')
        ax.set_ylabel('Temperature [K]')
        ax.set_title('Temperature Profiles')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Interface position vs time
        ax = axes[0, 1]
        ax.plot(results['time'], results['interface_position'])
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Interface Position [m]')
        ax.set_title('Solid-Liquid Interface Motion')
        ax.grid(True, alpha=0.3)
        
        # Phase change completion
        ax = axes[1, 0]
        liquid_volume = np.array([np.mean(1 - f) for f in results['temperature']])
        ax.plot(results['time'], liquid_volume * 100)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Liquid Fraction [%]')
        ax.set_title('Phase Change Progress')
        ax.grid(True, alpha=0.3)
        
        # Center temperature history
        ax = axes[1, 1]
        center_temp = [T[len(T)//2] for T in results['temperature']]
        ax.plot(results['time'], center_temp)
        ax.axhline(self.material.T_melt, color='k', linestyle='--', alpha=0.5)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Center Temperature [K]')
        ax.set_title('Center Temperature History')
        ax.grid(True, alpha=0.3)
        
        plt.suptitle(f'{self.material.name} - {self.geometry.shape.capitalize()} Geometry', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {save_path}")
        
        plt.show()


# Example usage
if __name__ == "__main__":
    # Define material (water)
    water = MaterialProperties(
        k=0.6,          # Thermal conductivity [W/(m·K)]
        c=4186,         # Specific heat [J/(kg·K)]
        L=334000,       # Latent heat [J/kg]
        rho=1000,       # Density [kg/m³]
        T_melt=273.15,  # Melting temperature [K]
        name="Water"
    )
    
    # Define geometry (planar slab)
    geometry = Geometry(
        shape='planar',
        L=0.01,         # Thickness [m]
        nodes=100
    )
    
    # Create solver
    solver = EnthalpySolver(water, geometry)
    
    # Initial condition (uniform temperature above melting)
    T_initial = np.ones(geometry.nodes) * 278.15  # 5°C above melting
    
    # Solve
    results = solver.solve(
        T_initial=T_initial,
        T_inf=253.15,  # -20°C
        h_eff=100.0,   # Convective coefficient [W/(m²·K)]
        t_final=3600   # 1 hour
    )
    
    # Plot results
    solver.plot_results(results, save_path="water_solidification.png")
