#!/usr/bin/env python3
"""
Unit tests for the numerical enthalpy solver.

Tests the correctness, stability, and accuracy of the
enthalpy-based finite difference solver.

Author: Jair Patiño B.
Date: January 2026
License: MIT
"""

import unittest
import numpy as np
import sys
import os
import tempfile
import shutil

# Add parent directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from scripts.numerical_solver.enthalpy_solver import (
    EnthalpySolver, MaterialProperties, Geometry
)
from scripts.numerical_solver.validation_simulations import ValidationStudy


class TestEnthalpySolverBasics(unittest.TestCase):
    """Basic tests for the enthalpy solver"""
    
    def setUp(self):
        """Set up test fixtures"""
        # Create a simple water material
        self.water = MaterialProperties(
            k=0.6,          # Thermal conductivity [W/(m·K)]
            c=4186,         # Specific heat [J/(kg·K)]
            L=334000,       # Latent heat [J/kg]
            rho=1000,       # Density [kg/m³]
            T_melt=273.15,  # Melting temperature [K]
            name="Water"
        )
        
        # Create planar geometry
        self.planar_geom = Geometry(
            shape='planar',
            L=0.01,         # 1 cm thickness
            nodes=50,       # Reduced for faster tests
            A=0.01,         # 100 cm² surface area
            V=1e-4          # Volume
        )
        
        # Create solver
        self.solver = EnthalpySolver(self.water, self.planar_geom)
    
    def test_initialization(self):
        """Test solver initialization"""
        self.assertEqual(self.solver.material.name, "Water")
        self.assertEqual(self.solver.geometry.shape, "planar")
        self.assertEqual(len(self.solver.T), self.planar_geom.nodes)
        self.assertEqual(len(self.solver.H), self.planar_geom.nodes)
        self.assertEqual(len(self.solver.f), self.planar_geom.nodes)
    
    def test_grid_setup(self):
        """Test computational grid setup"""
        # Grid should have correct number of points
        self.assertEqual(len(self.solver.x), self.planar_geom.nodes)
        
        # Grid should be evenly spaced
        dx = self.solver.x[1] - self.solver.x[0]
        for i in range(1, len(self.solver.x)):
            self.assertAlmostEqual(self.solver.x[i] - self.solver.x[i-1], dx,
                                 delta=1e-10,
                                 msg=f"Grid spacing should be uniform at position {i}")
        
        # Grid should span 0 to L
        self.assertAlmostEqual(self.solver.x[0], 0, delta=1e-10)
        self.assertAlmostEqual(self.solver.x[-1], self.planar_geom.L, delta=1e-10)
    
    def test_enthalpy_temperature_conversion(self):
        """Test conversion between temperature and enthalpy"""
        # Test temperatures below melting
        T_solid = np.array([263.15, 268.15])  # -10°C, -5°C
        H_solid, f_solid = self.solver._enthalpy_from_temp(T_solid)
        T_back, f_back = self.solver._temp_from_enthalpy(H_solid)
        
        np.testing.assert_array_almost_equal(T_solid, T_back, decimal=5)
        np.testing.assert_array_almost_equal(f_solid, f_back, decimal=5)
        self.assertTrue(np.all(f_solid == 0))
        
        # Test temperatures above melting
        T_liquid = np.array([278.15, 283.15])  # 5°C, 10°C
        H_liquid, f_liquid = self.solver._enthalpy_from_temp(T_liquid)
        T_back, f_back = self.solver._temp_from_enthalpy(H_liquid)
        
        np.testing.assert_array_almost_equal(T_liquid, T_back, decimal=5)
        np.testing.assert_array_almost_equal(f_liquid, f_back, decimal=5)
        self.assertTrue(np.all(f_liquid == 1))
        
        # Test in mushy zone
        T_mushy = np.array([272.9, 273.15, 273.4])  # Around melting point
        H_mushy, f_mushy = self.solver._enthalpy_from_temp(T_mushy)
        T_back, f_back = self.solver._temp_from_enthalpy(H_mushy)
        
        np.testing.assert_array_almost_equal(T_mushy, T_back, decimal=5)
        np.testing.assert_array_almost_equal(f_mushy, f_back, decimal=5)
        self.assertTrue(np.all((f_mushy >= 0) & (f_mushy <= 1)))
    
    def test_matrix_construction(self):
        """Test sparse matrix construction"""
        dt = 1.0
        h_eff = 50.0
        
        A = self.solver._build_matrix(dt, h_eff)
        
        # Matrix should be square
        n = self.planar_geom.nodes
        self.assertEqual(A.shape, (n, n))
        
        # Matrix should be sparse
        density = A.nnz / (n * n)
        self.assertLess(density, 0.1, "Matrix should be sparse")
        
        # Matrix should be diagonally dominant (for stability)
        diag = A.diagonal()
        row_sums = np.abs(A).sum(axis=1).A1
        for i in range(n):
            self.assertGreaterEqual(diag[i], row_sums[i] - diag[i],
                                  f"Row {i} should be diagonally dominant")


class TestSolverPhysics(unittest.TestCase):
    """Test physical correctness of the solver"""
    
    def setUp(self):
        self.water = MaterialProperties(
            k=0.6, c=4186, L=334000, rho=1000, T_melt=273.15, name="Water"
        )
        self.geometry = Geometry(shape='planar', L=0.01, nodes=30)
        self.solver = EnthalpySolver(self.water, self.geometry)
    
    def test_pure_conduction_no_phase_change(self):
        """Test pure conduction without phase change"""
        # Initial temperature above melting, ambient above melting
        T_initial = np.ones(self.geometry.nodes) * 278.15  # 5°C
        T_inf = 278.15  # Same as initial (no driving force)
        h_eff = 50.0
        
        # Run for short time
        results = self.solver.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=10.0,
            adaptive=False
        )
        
        # Temperature should remain constant (no driving force)
        final_temp = results['temperature'][-1]
        np.testing.assert_array_almost_equal(final_temp, T_initial, decimal=2)
    
    def test_sensible_cooling_only(self):
        """Test sensible cooling without phase change"""
        T_initial = np.ones(self.geometry.nodes) * 278.15  # 5°C
        T_inf = 263.15  # -10°C (below melting but material starts above melting)
        h_eff = 100.0  # Strong convection
        
        results = self.solver.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=50.0,
            adaptive=True
        )
        
        # Temperature should decrease
        final_temp = results['temperature'][-1]
        self.assertTrue(np.all(final_temp < T_initial))
        
        # Center should be warmer than surface
        center_temp = final_temp[len(final_temp)//2]
        surface_temp = final_temp[-1]
        self.assertGreater(center_temp, surface_temp,
                         "Center should be warmer than surface during cooling")
    
    def test_phase_change_progress(self):
        """Test that phase change progresses correctly"""
        T_initial = np.ones(self.geometry.nodes) * 278.15  # 5°C
        T_inf = 253.15  # -20°C
        h_eff = 50.0
        
        results = self.solver.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=200.0,
            adaptive=True
        )
        
        # Liquid fraction should decrease
        final_f = results.get('liquid_fraction', self.solver.f)
        if final_f is not None:
            # Some solidification should occur
            self.assertLess(np.mean(final_f), 0.8,
                          "Liquid fraction should decrease during solidification")
        
        # Interface should move inward
        interface_pos = results['interface_position']
        self.assertTrue(len(interface_pos) > 1)
        self.assertLess(interface_pos[-1], interface_pos[0],
                       "Interface should move from surface to center")
    
    def test_energy_conservation(self):
        """Test approximate energy conservation"""
        T_initial = np.ones(self.geometry.nodes) * 278.15  # 5°C
        T_inf = 253.15  # -20°C
        h_eff = 50.0
        
        results = self.solver.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=100.0,
            adaptive=True
        )
        
        # Calculate energy removed
        initial_energy = np.sum(self.solver._enthalpy_from_temp(T_initial)[0] * 
                               self.solver.vol)
        final_energy = np.sum(self.solver.H * self.solver.vol)
        
        # Energy change
        delta_E = initial_energy - final_energy
        
        # Heat transferred through boundaries (approximate)
        times = results['time']
        temps = results['temperature']
        
        Q_transferred = 0
        for i in range(len(times)-1):
            dt = times[i+1] - times[i]
            T_surf = temps[i][-1]  # Surface temperature
            Q_transferred += h_eff * self.geometry.A * abs(T_surf - T_inf) * dt
        
        # Energy change should be roughly equal to heat transferred
        # Allow 20% tolerance due to numerical approximations
        ratio = delta_E / Q_transferred
        self.assertGreater(ratio, 0.8, "Energy conservation: removed > 80% transferred")
        self.assertLess(ratio, 1.2, "Energy conservation: removed < 120% transferred")


class TestSolverStability(unittest.TestCase):
    """Test numerical stability and robustness"""
    
    def test_large_time_step(self):
        """Test stability with large time steps"""
        water = MaterialProperties(
            k=0.6, c=4186, L=334000, rho=1000, T_melt=273.15, name="Water"
        )
        geometry = Geometry(shape='planar', L=0.01, nodes=20)
        solver = EnthalpySolver(water, geometry)
        
        # Use intentionally large time step
        T_initial = np.ones(geometry.nodes) * 278.15
        T_inf = 253.15
        h_eff = 50.0
        
        # Should not crash
        results = solver.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=10.0,
            dt=10.0,  # Large time step
            adaptive=False
        )
        
        self.assertIsNotNone(results)
        self.assertGreater(len(results['time']), 1)
    
    def test_small_time_step(self):
        """Test accuracy with small time steps"""
        water = MaterialProperties(
            k=0.6, c=4186, L=334000, rho=1000, T_melt=273.15, name="Water"
        )
        geometry = Geometry(shape='planar', L=0.01, nodes=30)
        solver = EnthalpySolver(water, geometry)
        
        T_initial = np.ones(geometry.nodes) * 278.15
        T_inf = 253.15
        h_eff = 50.0
        
        # Run with small time step
        results_small = solver.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=20.0,
            dt=0.1,
            adaptive=False
        )
        
        # Run with adaptive time stepping
        solver2 = EnthalpySolver(water, geometry)
        results_adaptive = solver2.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=20.0,
            dt=None,  # Auto time step
            adaptive=True
        )
        
        # Results should be similar
        t1 = results_small['phase_change_time']
        t2 = results_adaptive['phase_change_time']
        
        # Within 10% difference
        diff = abs(t1 - t2) / max(t1, t2)
        self.assertLess(diff, 0.1, "Small dt and adaptive should give similar results")
    
    def test_extreme_temperature_differences(self):
        """Test solver with extreme temperature conditions"""
        water = MaterialProperties(
            k=0.6, c=4186, L=334000, rho=1000, T_melt=273.15, name="Water"
        )
        geometry = Geometry(shape='planar', L=0.01, nodes=25)
        solver = EnthalpySolver(water, geometry)
        
        # Very large temperature difference
        T_initial = np.ones(geometry.nodes) * 373.15  # 100°C
        T_inf = 173.15  # -100°C
        h_eff = 1000.0  # Very strong convection
        
        results = solver.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=5.0,  # Short time due to rapid cooling
            adaptive=True
        )
        
        self.assertIsNotNone(results)
        self.assertTrue(results['phase_change_time'] > 0)
    
    def test_different_geometries(self):
        """Test solver with different geometries"""
        water = MaterialProperties(
            k=0.6, c=4186, L=334000, rho=1000, T_melt=273.15, name="Water"
        )
        
        geometries = [
            ('planar', Geometry(shape='planar', L=0.01, nodes=20)),
            ('cylindrical', Geometry(shape='cylindrical', L=0.01, nodes=20)),
            ('spherical', Geometry(shape='spherical', L=0.01, nodes=20))
        ]
        
        times = []
        for name, geometry in geometries:
            solver = EnthalpySolver(water, geometry)
            
            T_initial = np.ones(geometry.nodes) * 278.15
            T_inf = 253.15
            h_eff = 50.0
            
            results = solver.solve(
                T_initial=T_initial,
                T_inf=T_inf,
                h_eff=h_eff,
                t_final=100.0,
                adaptive=True
            )
            
            times.append(results['phase_change_time'])
            
            # Should complete without error
            self.assertIsNotNone(results)
        
        # All times should be positive
        self.assertTrue(all(t > 0 for t in times))


class TestValidationStudy(unittest.TestCase):
    """Test the validation study framework"""
    
    def setUp(self):
        self.study = ValidationStudy()
    
    def test_analytical_prediction(self):
        """Test analytical prediction function"""
        # Simple case
        t_pred = self.study.analytical_prediction(
            'water', 'planar', 293.15, 253.15, 50.0
        )
        
        self.assertGreater(t_pred, 0)
        self.assertFalse(np.isnan(t_pred))
        self.assertFalse(np.isinf(t_pred))
    
    def test_material_properties(self):
        """Test material properties database"""
        materials = list(self.study.MATERIALS.keys())
        self.assertIn('water', materials)
        self.assertIn('paraffin', materials)
        self.assertIn('gallium', materials)
        
        # Check that all required properties exist
        for mat_name, props in self.study.MATERIALS.items():
            required = ['k', 'c', 'L', 'rho', 'T_melt']
            for prop in required:
                self.assertIn(prop, props, 
                            f"Material {mat_name} missing property {prop}")
                self.assertIsInstance(props[prop], (int, float),
                                    f"Property {prop} should be numeric")
    
    def test_geometry_configurations(self):
        """Test geometry configurations"""
        geometries = list(self.study.GEOMETRIES.keys())
        self.assertIn('planar', geometries)
        self.assertIn('cylindrical', geometries)
        self.assertIn('spherical', geometries)
        
        for geom_name, config in self.study.GEOMETRIES.items():
            required = ['shape', 'L', 'A', 'V']
            for prop in required:
                self.assertIn(prop, config,
                            f"Geometry {geom_name} missing property {prop}")
    
    def test_phi_factors(self):
        """Test geometric Φ factors"""
        phi = self.study.PHI
        self.assertEqual(phi['planar'], 1.0)
        self.assertEqual(phi['cylindrical'], 0.5)
        self.assertAlmostEqual(phi['spherical'], 1/3, delta=1e-10)
    
    def test_single_case_run(self):
        """Test running a single validation case"""
        with tempfile.TemporaryDirectory() as tmpdir:
            study = ValidationStudy(output_dir=tmpdir)
            
            case_result = study.run_single_case(
                case_id=999,
                material='water',
                geometry='planar',
                T_initial=293.15,
                T_inf=253.15,
                h_eff=50.0
            )
            
            # Check result structure
            required_keys = ['case_id', 'material', 'geometry', 'T_initial',
                           'T_inf', 'h_eff', 't_analytical', 't_numerical',
                           'error_percent', 'computation_time']
            
            for key in required_keys:
                self.assertIn(key, case_result,
                            f"Case result missing key {key}")
            
            # Check values
            self.assertEqual(case_result['case_id'], 999)
            self.assertEqual(case_result['material'], 'water')
            self.assertEqual(case_result['geometry'], 'planar')
            self.assertGreater(case_result['t_analytical'], 0)
            self.assertGreater(case_result['t_numerical'], 0)
            self.assertGreaterEqual(case_result['error_percent'], 0)
            self.assertGreater(case_result['computation_time'], 0)


class TestPerformanceBenchmarks(unittest.TestCase):
    """Performance benchmarks (not strict tests, but useful indicators)"""
    
    def test_solver_speed(self):
        """Benchmark solver speed for typical case"""
        import time
        
        water = MaterialProperties(
            k=0.6, c=4186, L=334000, rho=1000, T_melt=273.15, name="Water"
        )
        geometry = Geometry(shape='planar', L=0.01, nodes=100)
        solver = EnthalpySolver(water, geometry)
        
        T_initial = np.ones(geometry.nodes) * 278.15
        T_inf = 253.15
        h_eff = 50.0
        
        start = time.time()
        results = solver.solve(
            T_initial=T_initial,
            T_inf=T_inf,
            h_eff=h_eff,
            t_final=100.0,
            adaptive=True
        )
        elapsed = time.time() - start
        
        # Should complete in reasonable time (e.g., < 30 seconds for 100 nodes)
        # This is a soft constraint - just for benchmarking
        self.assertLess(elapsed, 30.0,
                       f"Solver took {elapsed:.1f}s for 100 nodes - might be slow")
        
        # Print performance info
        print(f"\nPerformance benchmark:")
        print(f"  Nodes: {geometry.nodes}")
        print(f"  Time steps: {results['time_steps']}")
        print(f"  Computation time: {results['computation_time']:.2f}s")
        print(f"  Time per step: {results['computation_time']/results['time_steps']*1000:.1f}ms")
    
    def test_memory_usage(self):
        """Check memory usage doesn't explode"""
        # Test with moderate grid size
        water = MaterialProperties(
            k=0.6, c=4186, L=334000, rho=1000, T_melt=273.15, name="Water"
        )
        
        # Try with 200 nodes (our standard validation size)
        geometry = Geometry(shape='planar', L=0.01, nodes=200)
        solver = EnthalpySolver(water, geometry)
        
        # Just creating the solver shouldn't use excessive memory
        # This is a qualitative test - we just ensure it doesn't crash
        self.assertIsNotNone(solver)
        self.assertEqual(len(solver.T), 200)
        
        print(f"\nMemory test: Created solver with {geometry.nodes} nodes")


def run_numerical_tests():
    """Run all numerical solver tests"""
    # Create test suites
    basic_suite = unittest.TestLoader().loadTestsFromTestCase(TestEnthalpySolverBasics)
    physics_suite = unittest.TestLoader().loadTestsFromTestCase(TestSolverPhysics)
    stability_suite = unittest.TestLoader().loadTestsFromTestCase(TestSolverStability)
    validation_suite = unittest.TestLoader().loadTestsFromTestCase(TestValidationStudy)
    perf_suite = unittest.TestLoader().loadTestsFromTestCase(TestPerformanceBenchmarks)
    
    # Combine suites
    combined_suite = unittest.TestSuite([
        basic_suite,
        physics_suite,
        stability_suite,
        validation_suite,
        perf_suite
    ])
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(combined_suite)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    print("Running numerical solver tests...")
    print("="*60)
    
    success = run_numerical_tests()
    
    print("\n" + "="*60)
    if success:
        print("All tests passed! ✓")
    else:
        print("Some tests failed ✗")
    
    sys.exit(0 if success else 1)
