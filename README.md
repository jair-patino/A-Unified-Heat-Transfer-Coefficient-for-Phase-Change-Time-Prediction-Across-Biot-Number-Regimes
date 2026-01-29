# A Unified Heat Transfer Coefficient for Phase Change Time Prediction Across Biot Number Regimes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![LaTeX](https://img.shields.io/badge/LaTeX-Article-brightgreen)

## Overview

This repository contains the complete code, data, and documentation for the research paper:

> **"A Unified Heat Transfer Coefficient for Phase Change Time Prediction Across Biot Number Regimes"**  
> **Author:** Jair Patiño B. (Independent Researcher)  
> **Date:** January 2026  
> **Status:** Submitted to Applied Physics Letters (APS)

The paper presents a novel analytical model for predicting phase change times in materials, valid across a wide range of Biot numbers (0.01 ≤ Bi ≤ 2). The model bridges the gap between lumped capacitance and moving boundary formulations, introducing a physically-motivated global heat transfer coefficient derived from thermal resistance analysis.

## Key Contributions

- **Unified Model**: Closed-form expression valid for 0.01 ≤ Bi ≤ 2
- **High Accuracy**: <5% error for Bi ≤ 1, <15% error at Bi = 2
- **Geometric Flexibility**: Derived Φ factors for planar, cylindrical, and spherical geometries
- **Supercooling Prediction**: Captures non-monotonic solidification time in water with minimum near -12°C
- **Rapid Optimization**: Enables parameter sweeps without numerical simulations

## Repository Structure
├── paper/ # LaTeX source for the manuscript
├── scripts/ # Python scripts for figures and simulations
├── data/ # Input data and simulation results
├── docs/ # Experimental protocol and supplementary material
└── tests/ # Unit tests for model validation
