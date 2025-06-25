# PBH Lensing Optical Depth Analysis with BURSTT

This repository contains code to compute and visualize the gravitational lensing optical depth for primordial black holes (PBHs), using cosmological models and FRB-like time delay constraints.

## Structure

- `scripts/tau_vs_mass_plot.py`  
  Plots **f<sub>DM</sub> vs lensing mass**, comparing against constraints from EROS, MACHO, FIRAS, etc.
  
- `scripts/integrated_tau_plot.py`  
  Plots **optical depth τ vs lensing mass**, for various time delays.

- `DMconstrain/`  
  Contains constraint data from external sources (EROS, MACHO, WMAP3, etc.).

## Dependencies

Install dependencies with:

```bash
pip install -r requirements.txt
```

## Usage
`python scripts/tau_vs_mass_plot.py` for the f<sub>DM</sub> vs lensing mass plot  

`python scripts/integrated_tau_plot.py` for the optical depth τ vs lensing mass plot

## Example plots
- f<sub>DM</sub> vs lens mass plot  
<img src="fdm_ml.png" alt="Batch 1" width="800">  

optical depth τ vs lens mass plot  
- <img src="" alt="Batch 2" width="800">

## Citation
Please cite [Ho et al. 2023, ApJ](https://iopscience.iop.org/article/10.3847/1538-4357/accb9e) if you use the code in your paper.

