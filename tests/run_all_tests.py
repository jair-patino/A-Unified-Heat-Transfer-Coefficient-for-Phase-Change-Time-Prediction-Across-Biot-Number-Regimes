#!/usr/bin/env python3
"""
Run all tests for the unified heat transfer model.

This script runs both model prediction tests and numerical solver tests,
providing a comprehensive test report.

Usage:
    python run_all_tests.py
    python run_all_tests.py --verbose
    python run_all_tests.py --coverage
"""

import sys
import os
import argparse
import unittest

# Add parent directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


def run_all_tests(verbosity=2, coverage=False):
    """Run all test suites"""
    from test_model_predictions import run_model_tests
    from test_numerical_solver import run_numerical_tests
    
    print("="*70)
    print("COMPREHENSIVE TEST SUITE - UNIFIED HEAT TRANSFER MODEL")
    print("="*70)
    
    # Run model tests
    print("\n" + "="*70)
    print("MODEL PREDICTION TESTS")
    print("="*70)
    model_success = run_model_tests()
    
    # Run numerical tests
    print("\n" + "="*70)
    print("NUMERICAL SOLVER TESTS")
    print("="*70)
    numerical_success = run_numerical_tests()
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    print(f"Model Prediction Tests: {'PASS' if model_success else 'FAIL'}")
    print(f"Numerical Solver Tests: {'PASS' if numerical_success else 'FAIL'}")
    
    overall_success = model_success and numerical_success
    print(f"\nOverall Result: {'ALL TESTS PASSED ✓' if overall_success else 'SOME TESTS FAILED ✗'}")
    print("="*70)
    
    return overall_success


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run all model tests')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Increase verbosity')
    parser.add_argument('--coverage', action='store_true',
                       help='Generate coverage report (requires pytest-cov)')
    parser.add_argument('--quick', action='store_true',
                       help='Run only critical tests')
    
    args = parser.parse_args()
    
    verbosity = 2 if args.verbose else 1
    
    success = run_all_tests(verbosity=verbosity, coverage=args.coverage)
    
    sys.exit(0 if success else 1)
