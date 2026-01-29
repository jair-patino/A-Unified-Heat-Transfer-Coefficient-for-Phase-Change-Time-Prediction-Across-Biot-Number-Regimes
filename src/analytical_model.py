"""
Analytical model for phase change time prediction across Biot regimes.
Implements unified heat transfer coefficient approach.
Author: Jair Patino B.
Date: January 2026
"""

import numpy as np
from typing import Union, Tuple

class UnifiedPhaseChangeModel:
    """
    Unified model for phase change time prediction.
    Valid for 0.01 ≤ Bi ≤ 2.
    """
    
    # Geometric factors Φ
    GEOMETRIC_FACTORS = {
        'planar': {'Φ': 1.0, 'Lc_factor': 0.5},      # Lc = thickness/2
        'cylinder': {'Φ': 0.5, 'Lc_factor': 0.5},    # Lc = radius/2
        'sphere': {'Φ': 1/3, 'Lc_factor': 1/3}       # Lc = radius/3
    }
    
    def __init__(self, geometry: str = 'planar'):
        """
        Initialize model for specific geometry.
        
        Parameters:
        -----------
        geometry : str
            'planar', 'cylinder', or 'sphere'
        """
        if geometry not in self.GEOMETRIC_FACTORS:
            raise ValueError(f"Geometry must be one of {list(self.GEOMETRIC_FACTORS.keys())}")
        
        self.geometry = geometry
        self.Φ = self.GEOMETRIC_FACTORS[geometry]['Φ']
        self.Lc_factor = self.GEOMETRIC_FACTORS[geometry]['Lc_factor']
    
    def calculate_characteristic_length(self, dimension: float) -> float:
        """
        Calculate characteristic length Lc = V/A.
        
        Parameters:
        -----------
        dimension : float
            For planar: thickness (m)
            For cylinder/sphere: radius (m)
        
        Returns:
        --------
        Lc : characteristic length (m)
        """
        return dimension * self.Lc_factor
    
    def calculate_global_U(self, h_eff: float, Bi: float) -> float:
        """
        Calculate global heat transfer coefficient U.
        
        Parameters:
        -----------
        h_eff : float
            Effective external heat transfer coefficient (W/m²·K)
        Bi : float
            Biot number
        
        Returns:
        --------
        U : global heat transfer coefficient (W/m²·K)
        """
        return h_eff / (1 + self.Φ * Bi)
    
    def calculate_Biot(self, h_eff: float, Lc: float, k: float) -> float:
        """
        Calculate Biot number.
        
        Parameters:
        -----------
        h_eff : float
            Effective external heat transfer coefficient (W/m²·K)
        Lc : float
            Characteristic length (m)
        k : float
            Thermal conductivity (W/m·K)
        
        Returns:
        --------
        Bi : Biot number
        """
        return h_eff * Lc / k
    
    def calculate_Stefan(self, c: float, delta_T: float, L: float) -> float:
        """
        Calculate Stefan number.
        
        Parameters:
        -----------
        c : float
            Specific heat (J/kg·K)
        delta_T : float
            Temperature difference |T∞ - Tf†| (K or °C)
        L : float
            Latent heat (J/kg)
        
        Returns:
        --------
        Ste : Stefan number
        """
        return c * abs(delta_T) / L
    
    def phase_change_time_total(self, 
                               m: float,           # mass (kg)
                               A: float,           # surface area (m²)
                               U: float,           # global heat transfer coefficient
                               c: float,           # specific heat (J/kg·K)
                               L: float,           # latent heat (J/kg)
                               T_initial: float,   # initial temperature
                               T_inf: float,       # ambient temperature
                               T_f_effective: float # effective phase change temperature
                               ) -> Tuple[float, float, float]:
        """
        Calculate total phase change time with sensible and latent contributions.
        
        Returns:
        --------
        t_total : total phase change time (s)
        t_sensible : sensible cooling/heating time (s)
        t_latent : latent phase change time (s)
        """
        # Ensure temperature differences are positive
        delta_T1 = abs(T_initial - T_inf)
        delta_T2 = abs(T_f_effective - T_inf)
        
        if delta_T1 <= 0 or delta_T2 <= 0:
            raise ValueError("Temperature differences must be positive")
        
        # Sensible time (Eq. 6)
        t_sensible = (m * c / (U * A)) * np.log(delta_T1 / delta_T2)
        
        # Latent time (Eq. 8)
        t_latent = (m * L) / (U * A * abs(T_inf - T_f_effective))
        
        # Total time (Eq. 9)
        t_total = t_sensible + t_latent
        
        return t_total, t_sensible, t_latent
    
    def dimensionless_phase_change_time(self, 
                                       Bi: float,
                                       Ste: float,
                                       Theta: float) -> float:
        """
        Calculate dimensionless phase change time (Fourier number).
        
        Parameters:
        -----------
        Bi : float
            Biot number
        Ste : float
            Stefan number
        Theta : float
            Dimensionless temperature difference |Ti - T∞| / |Tf† - T∞|
        
        Returns:
        --------
        FO_total : dimensionless total time
        """
        return (1 / (self.Φ * Bi * (1 + self.Φ * Bi))) * (np.log(Theta) + 1/Ste)
    
    def effective_phase_change_temperature(self,
                                          T_nucleation: float,
                                          T_equilibrium: float,
                                          alpha: float = 0.69) -> float:
        """
        Calculate effective phase change temperature for supercooling.
        
        Parameters:
        -----------
        T_nucleation : float
            Nucleation temperature (°C)
        T_equilibrium : float
            Equilibrium phase change temperature (°C)
        alpha : float
            Heat partitioning parameter, c_w/(c_w + c_i) ≈ 0.69 for water
        
        Returns:
        --------
        T_f_effective : effective phase change temperature (°C)
        """
        delta_T_sc = T_equilibrium - T_nucleation
        return T_equilibrium - alpha * delta_T_sc
    
    def correction_factor_high_Bi(self, Bi: float) -> float:
        """
        Correction factor for Bi > 2.
        
        Parameters:
        -----------
        Bi : float
            Biot number
        
        Returns:
        --------
        f_c : correction factor
        """
        if Bi <= 2:
            return 1.0
        else:
            return 1 + 0.08 * (Bi - 2)
    
    def calculate_mass_from_geometry(self, 
                                    dimension: float,
                                    rho: float,
                                    length: float = 1.0) -> Tuple[float, float]:
        """
        Calculate mass and surface area for given geometry.
        
        Parameters:
        -----------
        dimension : float
            For planar: thickness (m)
            For cylinder/sphere: radius (m)
        rho : float
            Density (kg/m³)
        length : float
            For planar/cylinder: length (m), for sphere: not used
        
        Returns:
        --------
        m : mass (kg)
        A : surface area (m²)
        """
        if self.geometry == 'planar':
            # Assume unit cross-sectional area
            A = 2.0  # Both sides
            V = dimension * 1.0 * 1.0  # thickness × width × depth (1 m²)
        elif self.geometry == 'cylinder':
            R = dimension
            A = 2 * np.pi * R * length
            V = np.pi * R**2 * length
        else:  # sphere
            R = dimension
            A = 4 * np.pi * R**2
            V = 4/3 * np.pi * R**3
        
        m = rho * V
        return m, A


def example_calculation():
    """
    Example usage of the unified model.
    """
    # Initialize model for planar geometry
    model = UnifiedPhaseChangeModel(geometry='planar')
    
    # Material properties (water)
    rho = 1000       # kg/m³
    k = 0.6          # W/m·K
    c = 4186         # J/kg·K
    L = 334000       # J/kg
    T_pc = 0.0       # °C (equilibrium)
    
    # Geometry
    thickness = 0.02  # 2 cm
    Lc = model.calculate_characteristic_length(thickness)
    
    # Heat transfer conditions
    h_eff = 100      # W/m²·K
    Bi = model.calculate_Biot(h_eff, Lc, k)
    
    # Calculate global U
    U = model.calculate_global_U(h_eff, Bi)
    
    # Calculate mass and area
    m, A = model.calculate_mass_from_geometry(thickness, rho)
    
    # Temperatures
    T_initial = 20.0    # °C
    T_inf = -10.0       # °C
    
    # For supercooling case
    T_nucleation = -8.0  # °C
    T_f_effective = model.effective_phase_change_temperature(
        T_nucleation, T_pc, alpha=0.69
    )
    
    # Calculate phase change time
    t_total, t_sensible, t_latent = model.phase_change_time_total(
        m, A, U, c, L, T_initial, T_inf, T_f_effective
    )
    
    # Print results
    print("=== Unified Phase Change Model Example ===")
    print(f"Geometry: {model.geometry}")
    print(f"Thickness: {thickness*100:.1f} cm")
    print(f"Biot number: {Bi:.3f}")
    print(f"Global U: {U:.2f} W/m²·K")
    print(f"Mass: {m:.3f} kg, Area: {A:.3f} m²")
    print(f"T_f_effective: {T_f_effective:.1f} °C")
    print(f"Sensible time: {t_sensible/60:.1f} min")
    print(f"Latent time: {t_latent/60:.1f} min")
    print(f"Total time: {t_total/60:.1f} min")
    
    # Dimensionless analysis
    Ste = model.calculate_Stefan(c, abs(T_inf - T_f_effective), L)
    Theta = abs(T_initial - T_inf) / abs(T_f_effective - T_inf)
    FO_total = model.dimensionless_phase_change_time(Bi, Ste, Theta)
    
    print(f"\nDimensionless analysis:")
    print(f"Stefan number: {Ste:.3f}")
    print(f"Theta: {Theta:.3f}")
    print(f"Fourier number: {FO_total:.3f}")
    
    return {
        'Bi': Bi, 'U': U, 't_total': t_total,
        't_sensible': t_sensible, 't_latent': t_latent
    }


if __name__ == "__main__":
    # Run example
    results = example_calculation()
    
    # Test different geometries
    print("\n=== Testing different geometries ===")
    
    test_cases = [
        ('planar', 0.02),
        ('cylinder', 0.01),
        ('sphere', 0.01)
    ]
    
    for geometry, dimension in test_cases:
        model = UnifiedPhaseChangeModel(geometry=geometry)
        Lc = model.calculate_characteristic_length(dimension)
        print(f"{geometry.capitalize()}: Lc = {Lc*100:.2f} cm, Φ = {model.Φ:.3f}")
