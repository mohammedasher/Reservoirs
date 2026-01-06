# Reservoir-based lithium-ion battery degradation model (code + data)

This repository contains the code and supporting data for the paper:

**“Lithium-ion battery degradation: Introducing the concept of reservoirs to design for lifetime”**  
Mohammed Asheruddin Nazeeruddin, Ruihe Li, Simon E. J. O’Kane, Monica Marinescu, Gregory J. Offer

## Overview

The repository implements a DFN/P2D (physics-based) lithium-ion cell model in **PyBaMM** with **five coupled degradation mechanisms**, and uses a “reservoir” framing (finite internal resources that are progressively consumed) to study how tuning design variables (e.g., lithium inventory, porosity, electrolyte volume) shifts degradation pathways and service life.

**PyBaMM compatibility:** this codebase is compatible with **PyBaMM v22 (22.x series)**.

## Repository layout

All content is inside the `Reservoir/` directory.

### Core code
- `Full_5Exp_5Ts_AddLi.py`  
  Main driver script. Reads an input bundle, runs the PyBaMM experiment(s), and writes outputs (Excel/MAT/PKL/plots).  
  Includes an HPC/PBS execution mode.
- `Fun_P2.py`  
  Core model + pipeline: assembles the DFN/P2D model with 5 coupled degradation mechanisms, executes PyBaMM experiments, post-processes outputs, and performs validation/diagnostics against the included experimental data.
- `OKane2023.py`  
  Parameter file used by the model (imported by the pipeline).
- `Custom_Para_Func.py`  
  Custom parameter/utility functions used by the model and post-processing.
- `Reservoir.pbs`  
  Example PBS job script for HPC execution (edit for your cluster paths/modules/resources).

### Inputs
- `InputData/Full_8Exps_AddLi_Pore/Bundle_1.csv`  
  Example “bundle” defining the run configuration and parameter set(s) (single-row CSV in this repository).
- `InputData/Expt 2,2 - C-based Degradation 2/`  
  Experimental extracted data for cells **A–F** (at 10°C / 25°C / 40°C) plus `DMA Output/` used by the post-processing/validation.

### Example outputs included
- `Full_8Exps_AddLi_Pore_Case_1_1/`  
  Example outputs corresponding to the default bundle/configuration.
  - `Excel/` : summary workbooks
  - `Mats/`  : intermediate/results files (`.mat`, `.pkl`, reload saves)
  - `Plots/` : generated PNG diagnostics (LLI/LAM breakdowns, electrolyte fields, half-cell potentials, etc.)

## Requirements

- **PyBaMM:** v22 (22.x)
- Python and packages commonly used by the scripts:
  - `numpy`, `pandas`, `scipy`, `matplotlib`
  - `openpyxl` (Excel writing)
  - `pyDOE` (Latin hypercube sampling utilities used in the pipeline)
