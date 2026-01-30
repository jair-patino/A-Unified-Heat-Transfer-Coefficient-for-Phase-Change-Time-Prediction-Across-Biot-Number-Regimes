#!/usr/bin/env python3
"""
Validation Simulations for the Unified Heat Transfer Model

This script performs systematic validation of the analytical model
against numerical solutions for various materials and conditions.

Author: Jair Patiño B.
Date: January 2026
License: MIT
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
import json
from pathlib import Path
from tqdm import tqdm
import seaborn as sns

from enthalpy_solver import EnthalpySolver, MaterialProperties, Geometry


class ValidationStudy:
    """Systematic validation of the analytical model"""
    
    # Pre-defined material properties
    MATERIALS = {
        'water': {
            'k': 0.6,      # Thermal conductivity [W/(m·K)]
            'c': 4186,     # Specific heat [J/(kg·K)]
            'L': 334000,   # Latent heat [J/kg]
            'rho': 1000,   # Density [kg/m³]
            'T_melt': 273.15  # Melting point [K]
        },
        'paraffin': {
            'k': 0.2,
            'c': 2200,
            'L': 216000,
            'rho': 900,
            'T_melt': 318.15  # 45°C
        },
        'gallium': {
            'k': 40.0,
            'c': 370,
            'L': 80100,
            'rho': 6090,
            'T_melt': 302.15  # 29°C
        },
        'hydrated_salt': {
            'k': 0.5,
            'c': 1400,
            'L': 170000,
            'rho': 1450,
            'T_melt': 302.15  # 29°C (CaCl2·6H2O)
        }
    }
    
    # Geometry configurations
    GEOMETRIES = {
        'planar': {'shape': 'planar', 'L': 0.01, 'A': 0.01, 'V': 1e-4},
        'cylindrical': {'shape': 'cylindrical', 'L': 0.01, 'A': 0.00628, 'V': 3.14e-5},
        'spherical': {'shape': 'spherical', 'L': 0.01, 'A': 0.001256, 'V': 4.19e-6}
    }
    
    def __init__(self, output_dir: str = "validation_results"):
        """
        Initialize validation study.
        
        Parameters
        ----------
        output_dir : str
            Directory to save results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize results storage
        self.results = []
        
        # Analytical model coefficients
        self.PHI = {
            'planar': 1.0,
            'cylindrical': 0.5,
            'spherical': 1/3
        }
    
    def analytical_prediction(self, material: str, geometry: str, 
                             T_initial: float, T_inf: float, 
                             h_eff: float) -> float:
        """
        Calculate phase change time using the analytical model.
        
        Parameters
        ----------
        material : str
            Material name
        geometry : str
            Geometry type
        T_initial : float
            Initial temperature [K]
        T_inf : float
            Ambient temperature [K]
        h_eff : float
            Effective heat transfer coefficient [W/(m²·K)]
            
        Returns
        -------
        t_total : float
            Predicted total phase change time [s]
        """
        props = self.MATERIALS[material]
        geom = self.GEOMETRIES[geometry]
        
        # Extract properties
        k = props['k']
        c = props['c']
        L = props['L']
        rho = props['rho']
        T_pc = props['T_melt']
        
        # Characteristic length
        L_c = geom['V'] / geom['A']
        
        # Biot number
        Bi = h_eff * L_c / k
        
        # Geometric factor
        phi = self.PHI[geometry]
        
        # Global heat transfer coefficient
        U = h_eff / (1 + phi * Bi)
        
        # Mass
        m = rho * geom['V']
        
        # Temperature differences (absolute values)
        delta_T1 = abs(T_initial - T_inf)
        delta_T2 = abs(T_pc - T_inf)
        
        # Sensible and latent times
        t_sensible = (m * c / (U * geom['A'])) * np.log(delta_T1 / delta_T2)
        t_latent = (m * L) / (U * geom['A'] * abs(T_inf - T_pc))
        
        return t_sensible + t_latent
    
    def run_numerical_simulation(self, material: str, geometry: str,
                                T_initial: float, T_inf: float,
                                h_eff: float, t_final: float = None) -> Dict:
        """
        Run numerical simulation for given conditions.
        
        Parameters
        ----------
        material : str
            Material name
        geometry : str
            Geometry type
        T_initial : float
            Initial temperature [K]
        T_inf : float
            Ambient temperature [K]
        h_eff : float
            Effective heat transfer coefficient [W/(m²·K)]
        t_final : float, optional
            Maximum simulation time [s]
            
        Returns
        -------
        results : Dict
            Numerical simulation results
        """
        props = self.MATERIALS[material]
        geom_config = self.GEOMETRIES[geometry]
        
        # Create material object
        mat = MaterialProperties(
            k=props['k'],
            c=props['c'],
            L=props['L'],
            rho=props['rho'],
            T_melt=props['T_melt'],
            name=material
        )
        
        # Create geometry object
        geom = Geometry(
            shape=geom_config['shape'],
            L=geom_config['L'],
            nodes=150
        )
        
        # Create solver
        solver = EnthalpySolver(mat, geom)
        
        # Estimate final time if not provided
        if t_final is None:
            # Rough estimate based on analytical model
            t_est = self.analytical_prediction(material, geometry, T_initial, T_inf, h_eff)
            t_final = t_est * 1.5  # Add 50% margin
        
        # Initial temperature distribution
        T_init = np.ones(geom.nodes) * T_initial
        
        # Run simulation
        results = solver.solve(
            T_initial=T_init,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=t_final,
            adaptive=True
        )
        
        return results
    
    def run_single_case(self, case_id: int, material: str, geometry: str,
                       T_initial: float, T_inf: float, h_eff: float) -> Dict:
        """
        Run a single validation case.
        
        Parameters
        ----------
        case_id : int
            Case identifier
        material : str
            Material name
        geometry : str
            Geometry type
        T_initial : float
            Initial temperature [K]
        T_inf : float
            Ambient temperature [K]
        h_eff : float
            Effective heat transfer coefficient [W/(m²·K)]
            
        Returns
        -------
        case_result : Dict
            Results for this case
        """
        print(f"\nRunning Case {case_id}: {material}, {geometry}")
        print(f"  T_initial: {T_initial:.1f}K, T_inf: {T_inf:.1f}K, h_eff: {h_eff:.1f}")
        
        # Analytical prediction
        t_analytical = self.analytical_prediction(material, geometry, T_initial, T_inf, h_eff)
        
        # Numerical simulation
        numerical_results = self.run_numerical_simulation(
            material, geometry, T_initial, T_inf, h_eff
        )
        
        t_numerical = numerical_results['phase_change_time']
        
        # Calculate error
        error = abs(t_analytical - t_numerical) / t_numerical * 100
        
        # Case result
        case_result = {
            'case_id': case_id,
            'material': material,
            'geometry': geometry,
            'T_initial': T_initial,
            'T_inf': T_inf,
            'h_eff': h_eff,
            't_analytical': t_analytical,
            't_numerical': t_numerical,
            'error_percent': error,
            'computation_time': numerical_results['computation_time']
        }
        
        print(f"  Analytical: {t_analytical:.1f}s, Numerical: {t_numerical:.1f}s")
        print(f"  Error: {error:.2f}%, Computation: {numerical_results['computation_time']:.1f}s")
        
        return case_result
    
    def run_bi_sweep(self, material: str = 'water', geometry: str = 'planar'):
        """
        Run Biot number sweep validation.
        
        Parameters
        ----------
        material : str
            Material name
        geometry : str
            Geometry type
        """
        print(f"\n{'='*60}")
        print(f"Biot Number Sweep: {material}, {geometry}")
        print('='*60)
        
        # Fixed conditions
        T_initial = 293.15  # 20°C
        T_inf = 253.15      # -20°C
        
        # Material properties
        k = self.MATERIALS[material]['k']
        L_c = self.GEOMETRIES[geometry]['V'] / self.GEOMETRIES[geometry]['A']
        
        # Biot number range
        Bi_values = np.logspace(-2, np.log10(2), 15)  # 0.01 to 2
        
        for i, Bi in enumerate(tqdm(Bi_values, desc="Bi sweep")):
            # Calculate h_eff from Bi
            h_eff = Bi * k / L_c
            
            # Run case
            case_result = self.run_single_case(
                case_id=i,
                material=material,
                geometry=geometry,
                T_initial=T_initial,
                T_inf=T_inf,
                h_eff=h_eff
            )
            
            self.results.append(case_result)
    
    def run_material_comparison(self, geometry: str = 'planar'):
        """
        Compare different materials at fixed conditions.
        
        Parameters
        ----------
        geometry : str
            Geometry type
        """
        print(f"\n{'='*60}")
        print(f"Material Comparison: {geometry} geometry")
        print('='*60)
        
        # Fixed conditions
        T_inf = 253.15  # -20°C
        h_eff = 50.0    # Moderate convection
        
        # Different materials and initial temperatures
        materials = {
            'water': 293.15,     # 20°C
            'paraffin': 333.15,  # 60°C
            'gallium': 313.15,   # 40°C
            'hydrated_salt': 313.15  # 40°C
        }
        
        for i, (material, T_initial) in enumerate(tqdm(materials.items(), desc="Materials")):
            case_result = self.run_single_case(
                case_id=100 + i,
                material=material,
                geometry=geometry,
                T_initial=T_initial,
                T_inf=T_inf,
                h_eff=h_eff
            )
            
            self.results.append(case_result)
    
    def run_geometry_comparison(self, material: str = 'water'):
        """
        Compare different geometries at fixed conditions.
        
        Parameters
        ----------
        material : str
            Material name
        """
        print(f"\n{'='*60}")
        print(f"Geometry Comparison: {material}")
        print('='*60)
        
        # Fixed conditions
        T_initial = 293.15  # 20°C
        T_inf = 253.15      # -20°C
        h_eff = 50.0        # Moderate convection
        
        geometries = ['planar', 'cylindrical', 'spherical']
        
        for i, geometry in enumerate(tqdm(geometries, desc="Geometries")):
            case_result = self.run_single_case(
                case_id=200 + i,
                material=material,
                geometry=geometry,
                T_initial=T_initial,
                T_inf=T_inf,
                h_eff=h_eff
            )
            
            self.results.append(case_result)
    
    def run_supercooling_study(self):
        """
        Study supercooling effects in water.
        """
        print(f"\n{'='*60}")
        print("Supercooling Study: Water")
        print('='*60)
        
        # Fixed geometry and material
        geometry = 'spherical'  # Droplet shape
        material = 'water'
        
        # Ambient temperature sweep
        T_inf_values = np.array([268.15, 263.15, 258.15, 253.15, 248.15])  # -5 to -25°C
        
        for i, T_inf in enumerate(tqdm(T_inf_values, desc="Supercooling")):
            # Initial temperature slightly above freezing
            T_initial = 275.15  # 2°C
            
            # Moderate convection
            h_eff = 20.0
            
            case_result = self.run_single_case(
                case_id=300 + i,
                material=material,
                geometry=geometry,
                T_initial=T_initial,
                T_inf=T_inf,
                h_eff=h_eff
            )
            
            self.results.append(case_result)
    
    def save_results(self):
        """Save all results to files."""
        # Save as CSV
        df = pd.DataFrame(self.results)
        csv_path = self.output_dir / "validation_results.csv"
        df.to_csv(csv_path, index=False)
        print(f"\nResults saved to {csv_path}")
        
        # Save as JSON
        json_path = self.output_dir / "validation_results.json"
        with open(json_path, 'w') as f:
            json.dump(self.results, f, indent=2, default=float)
        
        # Summary statistics
        summary = {
            'total_cases': len(self.results),
            'mean_error': df['error_percent'].mean(),
            'max_error': df['error_percent'].max(),
            'std_error': df['error_percent'].std(),
            'cases_below_5_percent': len(df[df['error_percent'] < 5]),
            'cases_below_10_percent': len(df[df['error_percent'] < 10]),
            'total_computation_time': df['computation_time'].sum()
        }
        
        summary_path = self.output_dir / "summary_statistics.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"Summary statistics saved to {summary_path}")
        print(f"\nValidation Summary:")
        print(f"  Total cases: {summary['total_cases']}")
        print(f"  Mean error: {summary['mean_error']:.2f}%")
        print(f"  Max error: {summary['max_error']:.2f}%")
        print(f"  Cases with <5% error: {summary['cases_below_5_percent']}")
        print(f"  Total computation time: {summary['total_computation_time']:.1f}s")
        
        return df
    
    def plot_validation_results(self, df: pd.DataFrame = None):
        """
        Plot validation results.
        
        Parameters
        ----------
        df : pd.DataFrame, optional
            Results dataframe. Loads from file if None.
        """
        if df is None:
            csv_path = self.output_dir / "validation_results.csv"
            df = pd.read_csv(csv_path)
        
        # Set style
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette("husl")
        
        # Create figure
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # 1. Error vs Biot number
        ax = axes[0, 0]
        planar_cases = df[df['geometry'] == 'planar']
        ax.scatter(planar_cases['h_eff'] * 0.01 / planar_cases['material'].map(
            lambda x: self.MATERIALS[x]['k']), planar_cases['error_percent'])
        ax.set_xscale('log')
        ax.set_xlabel('Biot Number (Bi)')
        ax.set_ylabel('Error [%]')
        ax.set_title('Error vs Biot Number (Planar Geometry)')
        ax.axhline(y=5, color='r', linestyle='--', alpha=0.5, label='5% threshold')
        ax.axhline(y=10, color='orange', linestyle='--', alpha=0.5, label='10% threshold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. Analytical vs Numerical comparison
        ax = axes[0, 1]
        max_time = df[['t_analytical', 't_numerical']].max().max()
        ax.scatter(df['t_analytical'], df['t_numerical'], c=df['error_percent'], 
                  cmap='viridis', alpha=0.7, edgecolors='k', linewidth=0.5)
        ax.plot([0, max_time], [0, max_time], 'r--', alpha=0.5, label='Perfect agreement')
        ax.set_xlabel('Analytical Prediction [s]')
        ax.set_ylabel('Numerical Solution [s]')
        ax.set_title('Analytical vs Numerical Phase Change Times')
        ax.set_aspect('equal')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add colorbar for error
        sm = plt.cm.ScalarMappable(cmap='viridis', 
                                  norm=plt.Normalize(vmin=df['error_percent'].min(), 
                                                    vmax=df['error_percent'].max()))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Error [%]')
        
        # 3. Material comparison
        ax = axes[1, 0]
        material_groups = df.groupby('material')
        materials = []
        mean_errors = []
        for name, group in material_groups:
            materials.append(name)
            mean_errors.append(group['error_percent'].mean())
        
        bars = ax.bar(materials, mean_errors)
        ax.set_xlabel('Material')
        ax.set_ylabel('Mean Error [%]')
        ax.set_title('Model Accuracy by Material')
        ax.axhline(y=5, color='r', linestyle='--', alpha=0.5)
        
        # Add value labels on bars
        for bar, error in zip(bars, mean_errors):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                   f'{error:.1f}%', ha='center', va='bottom', fontsize=9)
        
        # 4. Geometry comparison
        ax = axes[1, 1]
        geometry_groups = df.groupby('geometry')
        geometries = []
        geom_errors = []
        for name, group in geometry_groups:
            geometries.append(name)
            geom_errors.append(group['error_percent'].mean())
        
        bars = ax.bar(geometries, geom_errors)
        ax.set_xlabel('Geometry')
        ax.set_ylabel('Mean Error [%]')
        ax.set_title('Model Accuracy by Geometry')
        ax.axhline(y=5, color='r', linestyle='--', alpha=0.5)
        
        # Add value labels on bars
        for bar, error in zip(bars, geom_errors):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                   f'{error:.1f}%', ha='center', va='bottom', fontsize=9)
        
        plt.suptitle('Model Validation Results', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # Save figure
        fig_path = self.output_dir / "validation_summary.png"
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        print(f"Validation plot saved to {fig_path}")
        
        plt.show()
        
        # Additional plot: Error distribution
        fig2, ax2 = plt.subplots(figsize=(10, 6))
        ax2.hist(df['error_percent'], bins=20, edgecolor='black', alpha=0.7)
        ax2.set_xlabel('Error [%]')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Error Distribution')
        ax2.axvline(x=5, color='r', linestyle='--', alpha=0.7, label='5% threshold')
        ax2.axvline(x=10, color='orange', linestyle='--', alpha=0.7, label='10% threshold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        dist_path = self.output_dir / "error_distribution.png"
        plt.savefig(dist_path, dpi=300, bbox_inches='tight')
        plt.show()


def main():
    """Main validation routine"""
    print("="*70)
    print("COMPREHENSIVE VALIDATION OF UNIFIED HEAT TRANSFER MODEL")
    print("="*70)
    
    # Create validation study
    study = ValidationStudy(output_dir="validation_results")
    
    # Run different validation studies
    study.run_bi_sweep(material='water', geometry='planar')
    study.run_material_comparison(geometry='planar')
    study.run_geometry_comparison(material='water')
    study.run_supercooling_study()
    
    # Save and plot results
    df = study.save_results()
    study.plot_validation_results(df)
    
    print("\n" + "="*70)
    print("VALIDATION COMPLETED SUCCESSFULLY")
    print("="*70)


if __name__ == "__main__":
    main()
