# PBH Lensing Optical Depth Analysis

This repository contains code to compute and visualize the gravitational lensing optical depth for primordial black holes (PBHs), using cosmological models and FRB-like time delay constraints.

## Structure

- `scripts/tau_vs_mass_plot.py`  
  Plots **f<sub>DM</sub> vs lens mass**, comparing against constraints from EROS, MACHO, FIRAS, etc.
  
- `scripts/integrated_tau_plot.py`  
  Plots **optical depth τ vs lens mass**, for various time delays.

- `DMconstraint/`  
  Contains constraint data from external sources (EROS, MACHO, WMAP3, etc.).

## Dependencies

Install dependencies with:

```bash
pip install -r requirements.txt
