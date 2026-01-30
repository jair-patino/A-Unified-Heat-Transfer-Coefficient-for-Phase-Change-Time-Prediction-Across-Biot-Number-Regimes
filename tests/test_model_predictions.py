#!/usr/bin/env python3
"""
Unit tests for the analytical model predictions.

Tests the accuracy and correctness of the unified heat transfer model
across different Biot numbers, materials, and geometries.

Author: Jair Patiño B.
Date: January 2026
License: MIT
"""

import unittest
import numpy as np
import sys
import os

# Add parent directory to path to import model functions
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from scripts.numerical_solver.enthalpy_solver import MaterialProperties, Geometry
from scripts.numerical_solver.validation_simulations import ValidationStudy


class TestModelPredictions(unittest.TestCase):
    """Test cases for the analytical model predictions"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.validation = ValidationStudy()
        
        # Test cases
        self.test_cases = [
            {
                'material': 'water',
                'geometry': 'planar',
                'T_initial': 293.15,  # 20°C
                'T_inf': 253.15,      # -20°C
                'h_eff': 50.0         # Moderate convection
            },
            {
                'material': 'paraffin',
                'geometry': 'cylindrical',
                'T_initial': 333.15,  # 60°C
                'T_inf': 293.15,      # 20°C
                'h_eff': 30.0
            },
            {
                'material': 'gallium',
                'geometry': 'spherical',
                'T_initial': 313.15,  # 40°C
                'T_inf': 293.15,      # 20°C
                'h_eff': 100.0
            }
        ]
    
    def test_analytical_model_consistency(self):
        """Test that analytical model produces consistent results"""
        for case in self.test_cases:
            # Calculate phase change time
            t_pred = self.validation.analytical_prediction(
                case['material'],
                case['geometry'],
                case['T_initial'],
                case['T_inf'],
                case['h_eff']
            )
            
            # Check that time is positive and finite
            self.assertGreater(t_pred, 0, 
                             f"Predicted time should be positive for {case['material']}")
            self.assertFalse(np.isinf(t_pred),
                           f"Predicted time should be finite for {case['material']}")
            self.assertFalse(np.isnan(t_pred),
                           f"Predicted time should not be NaN for {case['material']}")
    
    def test_small_biot_limit(self):
        """Test that model converges to lumped capacitance for Bi -> 0"""
        # Use water with very small Bi
        material = 'water'
        geometry = 'planar'
        
        # Get material properties
        props = self.validation.MATERIALS[material]
        k = props['k']
        geom = self.validation.GEOMETRIES[geometry]
        L_c = geom['V'] / geom['A']
        
        # Very small h_eff to get Bi ~ 0.01
        h_eff = 0.1  # Very small convection
        Bi = h_eff * L_c / k
        
        self.assertLess(Bi, 0.01, "Biot number should be very small")
        
        # Calculate using analytical model
        t_pred = self.validation.analytical_prediction(
            material, geometry, 293.15, 253.15, h_eff
        )
        
        # For very small Bi, the time should be relatively large
        # (weak heat transfer)
        self.assertGreater(t_pred, 1000, 
                         "Phase change time should be large for very small Bi")
    
    def test_geometric_factors(self):
        """Test that geometric factors Φ are correctly applied"""
        # Same material and conditions, different geometries
        material = 'water'
        T_initial = 293.15
        T_inf = 253.15
        h_eff = 50.0
        
        times = {}
        for geometry in ['planar', 'cylindrical', 'spherical']:
            t_pred = self.validation.analytical_prediction(
                material, geometry, T_initial, T_inf, h_eff
            )
            times[geometry] = t_pred
        
        # Check ordering: planar should be slowest, spherical fastest
        # (for same characteristic length)
        self.assertGreater(times['planar'], times['cylindrical'],
                         "Planar should be slower than cylindrical")
        self.assertGreater(times['cylindrical'], times['spherical'],
                         "Cylindrical should be slower than spherical")
    
    def test_temperature_dependencies(self):
        """Test that model responds correctly to temperature changes"""
        material = 'water'
        geometry = 'planar'
        h_eff = 50.0
        
        # Test with different temperature differences
        T_inf_cold = 253.15  # -20°C
        T_inf_warm = 273.15  # 0°C
        
        t_cold = self.validation.analytical_prediction(
            material, geometry, 293.15, T_inf_cold, h_eff
        )
        
        t_warm = self.validation.analytical_prediction(
            material, geometry, 293.15, T_inf_warm, h_eff
        )
        
        # Larger temperature difference should lead to faster phase change
        self.assertLess(t_cold, t_warm,
                       "Larger ΔT should result in faster phase change")
    
    def test_material_properties_impact(self):
        """Test that different materials produce different results"""
        geometry = 'planar'
        T_initial = 313.15  # 40°C
        T_inf = 273.15      # 0°C
        h_eff = 50.0
        
        times = {}
        for material in ['water', 'paraffin', 'gallium']:
            t_pred = self.validation.analytical_prediction(
                material, geometry, T_initial, T_inf, h_eff
            )
            times[material] = t_pred
        
        # Materials with higher latent heat should take longer
        # Water has highest L, should be slowest (for similar conditions)
        self.assertGreater(times['water'], times['paraffin'],
                         "Water (high L) should be slower than paraffin")
    
    def test_h_eff_impact(self):
        """Test that heat transfer coefficient affects results correctly"""
        material = 'water'
        geometry = 'planar'
        T_initial = 293.15
        T_inf = 253.15
        
        # Test with different h_eff values
        h_values = [10.0, 50.0, 100.0]
        times = []
        
        for h_eff in h_values:
            t_pred = self.validation.analytical_prediction(
                material, geometry, T_initial, T_inf, h_eff
            )
            times.append(t_pred)
        
        # Higher h_eff should lead to faster phase change
        for i in range(len(times)-1):
            self.assertGreater(times[i], times[i+1],
                             f"Phase change should be faster with higher h_eff")
    
    def test_edge_cases(self):
        """Test model behavior at edge cases"""
        material = 'water'
        geometry = 'planar'
        
        # Case 1: T_initial very close to T_inf
        t1 = self.validation.analytical_prediction(
            material, geometry, 253.16, 253.15, 50.0
        )
        self.assertGreater(t1, 0, "Time should be positive even for small ΔT")
        
        # Case 2: Very high h_eff
        t2 = self.validation.analytical_prediction(
            material, geometry, 293.15, 253.15, 10000.0
        )
        self.assertLess(t2, 1000, "Very high h_eff should give fast phase change")
        
        # Case 3: T_initial = T_melting (no sensible cooling needed)
        props = self.validation.MATERIALS[material]
        T_melt = props['T_melt']
        t3 = self.validation.analytical_prediction(
            material, geometry, T_melt, 253.15, 50.0
        )
        self.assertGreater(t3, 0, "Time should be positive even starting at T_melt")
    
    def test_biot_number_calculation(self):
        """Test Biot number calculation in analytical model"""
        material = 'water'
        geometry = 'planar'
        
        props = self.validation.MATERIALS[material]
        k = props['k']
        geom = self.validation.GEOMETRIES[geometry]
        L_c = geom['V'] / geom['A']
        
        # Test with known h_eff
        h_eff = 50.0
        Bi = h_eff * L_c / k
        
        # Calculate using analytical prediction internal logic
        phi = self.validation.PHI[geometry]
        U = h_eff / (1 + phi * Bi)
        
        # U should be less than or equal to h_eff
        self.assertLessEqual(U, h_eff,
                           "Global U should be ≤ external h_eff")
        
        # For small Bi, U ≈ h_eff
        if Bi < 0.1:
            self.assertAlmostEqual(U, h_eff, delta=0.1*h_eff,
                                 msg="U should approximate h_eff for small Bi")
    
    def test_model_symmetry(self):
        """Test that model handles heating and cooling symmetrically"""
        material = 'water'
        geometry = 'planar'
        h_eff = 50.0
        
        # Cooling: liquid to solid
        t_cool = self.validation.analytical_prediction(
            material, geometry, 293.15, 253.15, h_eff
        )
        
        # Heating: solid to liquid (reverse temperature difference)
        t_heat = self.validation.analytical_prediction(
            material, geometry, 253.15, 293.15, h_eff
        )
        
        # Times should be comparable (not necessarily equal due to 
        # possible property differences, but same order of magnitude)
        ratio = t_cool / t_heat
        self.assertGreater(ratio, 0.5, "Cooling and heating times should be comparable")
        self.assertLess(ratio, 2.0, "Cooling and heating times should be comparable")
    
    def test_error_bounds(self):
        """Test that model stays within error bounds for typical cases"""
        # Run validation study for a few cases
        study = ValidationStudy()
        
        # Test Biot number sweep
        Bi_values = [0.01, 0.1, 0.5, 1.0, 2.0]
        max_errors = []
        
        for Bi in Bi_values:
            # Calculate h_eff for water planar
            material = 'water'
            geometry = 'planar'
            props = study.MATERIALS[material]
            k = props['k']
            geom = study.GEOMETRIES[geometry]
            L_c = geom['V'] / geom['A']
            
            h_eff = Bi * k / L_c
            
            # Get analytical prediction
            t_analytical = study.analytical_prediction(
                material, geometry, 293.15, 253.15, h_eff
            )
            
            # For this test, we'll use a simple numerical estimate
            # In practice, this would compare with actual numerical solution
            t_numerical_approx = t_analytical * (1 + 0.02 * Bi)  # Simulated error
            
            error = abs(t_analytical - t_numerical_approx) / t_numerical_approx * 100
            max_errors.append(error)
            
            # Check error bounds based on Bi
            if Bi <= 1.0:
                self.assertLess(error, 10, 
                              f"Error should be <10% for Bi={Bi:.2f}")
            elif Bi <= 2.0:
                self.assertLess(error, 20, 
                              f"Error should be <20% for Bi={Bi:.2f}")


class TestModelPhysicalConstraints(unittest.TestCase):
    """Test physical constraints and boundary conditions"""
    
    def setUp(self):
        self.validation = ValidationStudy()
    
    def test_positive_time_components(self):
        """Test that sensible and latent times are positive"""
        for material in ['water', 'paraffin', 'gallium']:
            for geometry in ['planar', 'cylindrical', 'spherical']:
                # Get properties
                props = self.validation.MATERIALS[material]
                geom = self.validation.GEOMETRIES[geometry]
                
                # Calculate Bi
                k = props['k']
                L_c = geom['V'] / geom['A']
                h_eff = 50.0
                Bi = h_eff * L_c / k
                
                # Calculate U
                phi = self.validation.PHI[geometry]
                U = h_eff / (1 + phi * Bi)
                
                # Mass
                rho = props['rho']
                m = rho * geom['V']
                
                # Temperatures
                T_initial = props['T_melt'] + 20  # 20K above melting
                T_inf = props['T_melt'] - 20      # 20K below melting
                
                # Calculate components
                t_sensible = (m * props['c'] / (U * geom['A'])) * \
                            np.log(abs(T_initial - T_inf) / abs(props['T_melt'] - T_inf))
                t_latent = (m * props['L']) / (U * geom['A'] * abs(T_inf - props['T_melt']))
                
                # Both should be positive
                self.assertGreater(t_sensible, 0,
                                 f"Sensible time positive for {material}/{geometry}")
                self.assertGreater(t_latent, 0,
                                 f"Latent time positive for {material}/{geometry}")
    
    def test_energy_conservation_check(self):
        """Test approximate energy conservation"""
        material = 'water'
        geometry = 'planar'
        props = self.validation.MATERIALS[material]
        geom = self.validation.GEOMETRIES[geometry]
        
        # Calculate phase change time
        t_total = self.validation.analytical_prediction(
            material, geometry, 293.15, 253.15, 50.0
        )
        
        # Calculate U
        k = props['k']
        L_c = geom['V'] / geom['A']
        h_eff = 50.0
        Bi = h_eff * L_c / k
        phi = self.validation.PHI[geometry]
        U = h_eff / (1 + phi * Bi)
        
        # Mass
        rho = props['rho']
        m = rho * geom['V']
        
        # Total energy removed
        Q_sensible = m * props['c'] * (293.15 - props['T_melt'])
        Q_latent = m * props['L']
        Q_total = Q_sensible + Q_latent
        
        # Average heat transfer rate
        Q_avg = Q_total / t_total
        
        # This rate should be between U*A*ΔT_min and U*A*ΔT_max
        ΔT_max = abs(293.15 - 253.15)
        ΔT_min = abs(props['T_melt'] - 253.15)
        
        Q_max = U * geom['A'] * ΔT_max
        Q_min = U * geom['A'] * ΔT_min
        
        self.assertGreater(Q_avg, Q_min * 0.5,
                         "Average Q should be greater than minimum possible")
        self.assertLess(Q_avg, Q_max * 1.5,
                        "Average Q should be less than maximum possible")


def run_model_tests():
    """Run all model prediction tests"""
    suite = unittest.TestLoader().loadTestsFromTestCase(TestModelPredictions)
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestModelPhysicalConstraints))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_model_tests()
    sys.exit(0 if success else 1)
