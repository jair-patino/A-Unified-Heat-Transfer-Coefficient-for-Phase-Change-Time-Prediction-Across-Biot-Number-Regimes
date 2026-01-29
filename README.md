# Unified Phase Change Time Prediction Model

This repository contains the code and data for the paper:

**"A Unified Heat Transfer Coefficient for Phase Change Time Prediction Across Biot Number Regimes"**  
*Jair Patino B., 2026*  
*(Submitted to Physical Review Applied)*

## Abstract
Accurate prediction of phase change times governs efficiency in applications ranging from thermal energy storage to cryopreservation. While exact solutions require numerical methods, existing analytical approximations fail in the intermediate regime where internal and external thermal resistances are comparable (Bi ∼ 1). Here, we bridge this gap by deriving a closed-form expression valid for 0.01 ≤ Bi ≤ 2. Our formulation introduces a physically motivated global heat transfer coefficient U derived from first-principles thermal resistance analysis, with geometrically-derived Φ factors for planar, cylindrical, and spherical geometries.

## Features
- Unified analytical model bridging lumped capacitance and moving boundary formulations
- Closed-form solution with <5% error for Bi ≤ 1 and <15% error at Bi = 2
- Quantitative prediction of supercooling effects in water solidification
- Rapid parameter sweeps for system optimization without full numerical simulations

## Installation

### Using pip:
```bash
git clone https://github.com/jairpatino/phase-change-model.git
cd phase-change-model
pip install -r requirements.txt
