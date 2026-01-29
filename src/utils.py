"""
Utility functions for phase change model.
Thermophysical properties, conversions, and helper functions.
Author: Jair Patino B.
Date: January 2026
"""

import numpy as np
from typing import Dict, Union, Optional
import warnings

# Constants
R_GAS_CONSTANT = 8.314462618  # J/mol·K
AVOGADRO = 6.02214076e23  # mol⁻¹

def celsius_to_kelvin(T_C: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert Celsius to Kelvin."""
    return T_C + 273.15

def kelvin_to_celsius(T_K: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Convert Kelvin to Celsius."""
    return T_K - 273.15

def calculate_effective_h(h_conv: float, 
                         h_rad: Optional[float] = None,
                         h_cond_ext: Optional[float] = None) -> float:
    """
    Calculate effective external heat transfer coefficient.
    
    Parameters:
    -----------
    h_conv : float
        Convective heat transfer coefficient (W/m²·K)
    h_rad : float, optional
        Radiative heat transfer coefficient (W/m²·K)
    h_cond_ext : float, optional
        External conductive resistance (W/m²·K)
    
    Returns:
    --------
    h_eff : effective heat transfer coefficient (W/m²·K)
    """
    h_eff = h_conv
    
    if h_rad is not None:
        h_eff += h_rad
    
    if h_cond_ext is not None:
        # External conduction resistance in series
        h_eff = 1 / (1/h_eff + 1/h_cond_ext)
    
    return h_eff

def water_properties(T: float = 20.0, phase: str = 'liquid') -> Dict[str, float]:
    """
    Get thermophysical properties of water/ice.
    
    Parameters:
    -----------
    T : float
        Temperature (°C)
    phase : str
        'liquid' or 'solid'
    
    Returns:
    --------
    dict with properties
    """
    T = float(T)
    
    if phase.lower() == 'liquid':
        # Water properties (approximate, from 0-100°C)
        rho = 1000 - 0.2 * (T - 4)  # Maximum density at 4°C
        k = 0.6 + 0.001 * T  # W/m·K
        c = 4186 - 1.5 * T  # J/kg·K
        alpha = k / (rho * c)  # Thermal diffusivity
        beta = 0.000214  # Thermal expansion coefficient (1/K)
        
        return {
            'rho': rho,          # kg/m³
            'k': k,              # W/m·K
            'c': c,              # J/kg·K
            'alpha': alpha,      # m²/s
            'beta': beta,        # 1/K
            'phase': 'liquid',
            'T': T
        }
    
    elif phase.lower() == 'solid':
        # Ice properties (approximate, from -50 to 0°C)
        rho = 917 - 0.1 * (T + 10)  # kg/m³
        k = 2.2 - 0.01 * T  # W/m·K
        c = 2100 + 7.5 * T  # J/kg·K
        alpha = k / (rho * c)
        
        return {
            'rho': rho,
            'k': k,
            'c': c,
            'alpha': alpha,
            'phase': 'solid',
            'T': T
        }
    
    else:
        raise ValueError("Phase must be 'liquid' or 'solid'")

def paraffin_properties(T: float = 25.0, grade: str = 'RT25') -> Dict[str, float]:
    """
    Get properties of paraffin wax (common PCM).
    
    Parameters:
    -----------
    T : float
        Temperature (°C)
    grade : str
        Paraffin grade ('RT25', 'RT42', 'RT50')
    
    Returns:
    --------
    dict with properties
    """
    grades = {
        'RT25': {
            'T_melt': 25.0,      # Melting point (°C)
            'L': 180000,         # Latent heat (J/kg)
            'rho_l': 760,        # Liquid density (kg/m³)
            'rho_s': 880,        # Solid density (kg/m³)
            'c_l': 2200,         # Liquid specific heat (J/kg·K)
            'c_s': 1800,         # Solid specific heat (J/kg·K)
            'k_l': 0.15,         # Liquid thermal conductivity (W/m·K)
            'k_s': 0.2,          # Solid thermal conductivity (W/m·K)
        },
        'RT42': {
            'T_melt': 42.0,
            'L': 200000,
            'rho_l': 780,
            'rho_s': 900,
            'c_l': 2400,
            'c_s': 1900,
            'k_l': 0.17,
            'k_s': 0.22,
        }
    }
    
    if grade not in grades:
        warnings.warn(f"Grade {grade} not found, using RT25")
        grade = 'RT25'
    
    props = grades[grade].copy()
    
    # Add temperature-dependent density
    if T < props['T_melt']:
        props['rho'] = props['rho_s']
        props['k'] = props['k_s']
        props['c'] = props['c_s']
        props['phase'] = 'solid'
    else:
        props['rho'] = props['rho_l']
        props['k'] = props['k_l']
        props['c'] = props['c_l']
        props['phase'] = 'liquid'
    
    props['T'] = T
    props['grade'] = grade
    
    return props

def gallium_properties(T: float = 30.0) -> Dict[str, float]:
    """
    Get properties of gallium (metal PCM).
    """
    T_melt = 29.76  # °C
    
    if T < T_melt:
        # Solid gallium
        return {
            'rho': 5900,          # kg/m³
            'k': 40.6,           # W/m·K
            'c': 370,            # J/kg·K
            'L': 80100,          # J/kg
            'T_melt': T_melt,
            'phase': 'solid',
            'T': T
        }
    else:
        # Liquid gallium
        return {
            'rho': 6095,          # kg/m³
            'k': 30.0,           # W/m·K
            'c': 380,            # J/kg·K
            'L': 80100,          # J/kg
            'T_melt': T_melt,
            'phase': 'liquid',
            'T': T
        }

def calculate_nusselt_number(Re: float, Pr: float, geometry: str = 'sphere') -> float:
    """
    Calculate Nusselt number for convection correlations.
    
    Parameters:
    -----------
    Re : float
        Reynolds number
    Pr : float
        Prandtl number
    geometry : str
        'sphere', 'cylinder', or 'plate'
    
    Returns:
    --------
    Nu : Nusselt number
    """
    if geometry.lower() == 'sphere':
        # Ranz-Marshall correlation for spheres
        if Re < 1e4:
            Nu = 2 + 0.6 * Re**0.5 * Pr**0.33
        else:
            Nu = 2 + (0.4 * Re**0.5 + 0.06 * Re**0.67) * Pr**0.33
    
    elif geometry.lower() == 'cylinder':
        # Churchill-Bernstein correlation
        Re_Pr = Re * Pr
        if Re_Pr < 0.2:
            Nu = 0.3 + (0.62 * Re**0.5 * Pr**0.33) / (1 + (0.4/Pr)**0.67)**0.25
        else:
            Nu = 0.3 + (0.62 * Re**0.5 * Pr**0.33) / (1 + (0.4/Pr)**0.67)**0.25 * (
                1 + (Re/282000)**0.5)**0.8
    
    elif geometry.lower() == 'plate':
        # Flat plate correlation
        if Re < 5e5:  # Laminar
            Nu = 0.664 * Re**0.5 * Pr**0.33
        else:  # Turbulent
            Nu = (0.037 * Re**0.8 - 871) * Pr**0.33
    
    else:
        raise ValueError("Geometry must be 'sphere', 'cylinder', or 'plate'")
    
    return float(Nu)

def calculate_radiative_h(epsilon: float, 
                         T_surface: float, 
                         T_surroundings: float = 293.15) -> float:
    """
    Calculate radiative heat transfer coefficient.
    
    Parameters:
    -----------
    epsilon : float
        Emissivity (0-1)
    T_surface : float
        Surface temperature (K)
    T_surroundings : float
        Surroundings temperature (K)
    
    Returns:
    --------
    h_rad : radiative heat transfer coefficient (W/m²·K)
    """
    sigma = 5.670374419e-8  # Stefan-Boltzmann constant
    
    h_rad = epsilon * sigma * (T_surface**2 + T_surroundings**2) * \
            (T_surface + T_surroundings)
    
    return h_rad

def estimate_convective_h(fluid: str = 'air',
                         flow_type: str = 'natural',
                         delta_T: float = 10.0,
                         L: float = 0.1,
                         velocity: Optional[float] = None) -> float:
    """
    Estimate convective heat transfer coefficient.
    
    Parameters:
    -----------
    fluid : str
        'air', 'water', 'oil'
    flow_type : str
        'natural' or 'forced'
    delta_T : float
        Temperature difference (K)
    L : float
        Characteristic length (m)
    velocity : float, optional
        Flow velocity for forced convection (m/s)
    
    Returns:
    --------
    h_conv : convective heat transfer coefficient (W/m²·K)
    """
    # Typical ranges (W/m²·K)
    typical_values = {
        'air': {
            'natural': 5-25,
            'forced': 10-200
        },
        'water': {
            'natural': 100-1000,
            'forced': 500-10000
        },
        'oil': {
            'natural': 50-500,
            'forced': 200-2000
        }
    }
    
    if fluid not in typical_values:
        raise ValueError(f"Fluid must be one of {list(typical_values.keys())}")
    
    if flow_type not in ['natural', 'forced']:
        raise ValueError("Flow type must be 'natural' or 'forced'")
    
    ranges = typical_values[fluid][flow_type]
    
    if isinstance(ranges, tuple) or isinstance(ranges, list):
        # Return midpoint of range
        h_est = (ranges[0] + ranges[1]) / 2
    else:
        h_est = ranges
    
    # Adjust for temperature difference and length
    if flow_type == 'natural':
        # Natural convection scales with delta_T^0.25 and L^-0.25
        h_est *= (delta_T / 10)**0.25 * (0.1 / L)**0.25
    else:
        # Forced convection scales with velocity^0.8 and L^-0.2
        if velocity is not None:
            h_est *= (velocity / 1.0)**0.8 * (0.1 / L)**0.2
    
    return float(h_est)

def monte_carlo_sensitivity(func, params, n_samples=1000):
    """
    Perform Monte Carlo sensitivity analysis.
    
    Parameters:
    -----------
    func : callable
        Function to analyze
    params : dict
        Dictionary with parameter names and (mean, std) tuples
    n_samples : int
        Number of Monte Carlo samples
    
    Returns:
    --------
    results : dict with sensitivity statistics
    """
    # Generate random samples
    samples = {}
    for name, (mean, std) in params.items():
        samples[name] = np.random.normal(mean, std, n_samples)
    
    # Evaluate function
    outputs = []
    for i in range(n_samples):
        args = {name: samples[name][i] for name in params.keys()}
        outputs.append(func(**args))
    
    outputs = np.array(outputs)
    
    # Calculate statistics
    results = {
        'mean': np.mean(outputs),
        'std': np.std(outputs),
        'min': np.min(outputs),
        'max': np.max(outputs),
        'cv': np.std(outputs) / np.mean(outputs) * 100,  # Coefficient of variation (%)
        'samples': outputs
    }
    
    # Calculate sensitivity coefficients
    sensitivities = {}
    for name in params.keys():
        # Correlation coefficient
        corr = np.corrcoef(samples[name], outputs)[0, 1]
        sensitivities[name] = {
            'correlation': corr,
            'importance': abs(corr)  # Absolute value for importance ranking
        }
    
    results['sensitivities'] = sensitivities
    
    return results

def validate_inputs(**kwargs):
    """
    Validate input parameters for phase change calculations.
    """
    errors = []
    warnings = []
    
    # Check Biot number
    if 'Bi' in kwargs:
        Bi = kwargs['Bi']
        if Bi < 0.01:
            warnings.append(f"Bi={Bi:.3f} is very small, consider using lumped capacitance")
        elif Bi > 2:
            warnings.append(f"Bi={Bi:.3f} > 2, model accuracy decreases")
    
    # Check temperatures
    temp_keys = ['T_initial', 'T_inf', 'T_pc']
    for key in temp_keys:
        if key in kwargs:
            T = kwargs[key]
            if not (-273.15 <= T <= 10000):  # Reasonable range
                warnings.append(f"{key}={T} seems unusual")
    
    # Check positive values
    positive_keys = ['m', 'A', 'c', 'L', 'k', 'rho']
    for key in positive_keys:
        if key in kwargs:
            val = kwargs[key]
            if val <= 0:
                errors.append(f"{key} must be positive, got {val}")
    
    # Check geometry factors
    if 'Phi' in kwargs:
        Phi = kwargs['Phi']
        if not (1/3 <= Phi <= 1):
            warnings.append(f"Φ={Phi} outside typical range [1/3, 1]")
    
    return errors, warnings

def format_scientific(value, decimals=3):
    """
    Format number in scientific notation for display.
    """
    if value == 0:
        return "0"
    
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / 10**exponent
    
    if abs(exponent) <= 2:
        return f"{value:.{decimals}f}"
    else:
        return f"{mantissa:.{decimals}f} × 10^{exponent}"

if __name__ == "__main__":
    # Example usage
    print("=== Utility Functions Examples ===")
    
    # Temperature conversion
    T_C = 20.0
    T_K = celsius_to_kelvin(T_C)
    print(f"{T_C}°C = {T_K:.2f} K")
    
    # Water properties
    water_props = water_properties(25.0, 'liquid')
    print(f"\nWater at 25°C:")
    for key, val in water_props.items():
        if isinstance(val, float):
            print(f"  {key}: {val:.3f}")
    
    # Effective h calculation
    h_eff = calculate_effective_h(h_conv=10, h_rad=5)
    print(f"\nEffective h (conv=10, rad=5): {h_eff:.2f} W/m²·K")
    
    # Validate inputs
    errors, warnings = validate_inputs(Bi=3.0, m=1.0, A=0.1, c=2000, L=100000)
    if warnings:
        print(f"\nWarnings: {warnings}")
    if errors:
        print(f"Errors: {errors}")
    
    # Monte Carlo example
    print("\n=== Monte Carlo Sensitivity Example ===")
    
    def sample_func(x, y):
        return 2*x + 3*y
    
    params = {
        'x': (10, 1),  # mean=10, std=1
        'y': (5, 0.5)   # mean=5, std=0.5
    }
    
    results = monte_carlo_sensitivity(sample_func, params, n_samples=100)
    print(f"Output mean: {results['mean']:.2f}")
    print(f"Output std: {results['std']:.2f}")
    
    print("\nParameter sensitivities:")
    for param, sens in results['sensitivities'].items():
        print(f"  {param}: correlation={sens['correlation']:.3f}")
