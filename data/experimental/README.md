# Experimental Data

This directory contains experimental data used for model validation.

## Files:

### 1. mueller_2015_supercooling.csv
Experimental data from:
Müller, K., Paschinger, H., & Fröba, A. P. (2015). 
Supercooling of water droplets in jet fuel. 
Fuel, 148, 16-24.

Columns:
- ambient_temp_C: Ambient temperature in °C
- exp_time_s: Experimental solidification time in seconds
- exp_time_error_s: Experimental uncertainty in seconds
- droplet_diameter_mm: Droplet diameter in millimeters
- purity: Water purity level (low/moderate/high)
- supercooling_degree_C: Degree of supercooling (T_pc - T_n) in °C
- source: Reference source

## Notes:
- Data points with moderate purity and 1.0 mm droplets were used for main validation.
- Supercooling degree of 8°C corresponds to nucleation at -8°C (T_n = -8°C).
- Minimum solidification time observed near -12°C.
