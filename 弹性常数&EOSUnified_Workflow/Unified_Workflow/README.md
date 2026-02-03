# Unified EOS & Elastic Workflow

This directory contains a streamlined workflow for calculating EOS and Elastic constants using both VASP and LAMMPS.

## Directory Structure
- `VASP/`: Scripts and templates for VASP calculations.
- `LAMMPS/`: Scripts for LAMMPS calculations.
- `Analysis/`: Scripts for plotting and summarizing results.
- `config.json`: Master configuration file for both workflows.

## Prerequisites
1.  **Potentials**: Ensure `mlip.ini` and `current.mtp` are available (checked automatically in parent folders).
2.  **POTCAR**: For VASP, ensure `vaspkit` is configured or manually copy `POTCAR` to the run folders.
3.  **Python Env**: Requires `ase`, `numpy`, `scipy`, `pandas`, `matplotlib`.

## Usage Steps

### 1. Configuration
Edit `config.json` to set:
- Elements and Masses.
- EOS scale factors (applied to both VASP and LAMMPS).
- Paths to binaries (`vasp_std`, `lmp_...`).
- Queue settings.

### 2. VASP Workflow
```bash
cd VASP
python setup_vasp.py <input.xyz>
# This creates struct_0, struct_1... folders.
# Submit jobs (using generated scripts or submit_batch.py)
python submit_batch.py submission_scripts
```
Each VASP job will:
1.  Relax the structure (`INCAR_relax`).
2.  Calculate EOS points (`INCAR_eos`, `ISIF=2`).
3.  Calculate Elastic Constants (`INCAR_elastic`, `IBRION=6`).

### 3. LAMMPS Workflow
```bash
cd LAMMPS
python setup_lammps.py <input.xyz>
# Submit jobs
python submit_batch.py submission_scripts
```
Each LAMMPS job will:
1.  Relax the structure (`min_style cg`).
2.  Calculate EOS points (Static, Scaled).
3.  Calculate Elastic Constants (Stress-Strain).

### 4. Analysis & Plotting
After all jobs are finished:
```bash
cd Analysis
python plot_all.py
```
This will generate the following files in the `results/` folder:
- `eos_raw_data_vasp.csv`: VASP Raw Energy-Volume data.
- `eos_raw_data_lammps.csv`: LAMMPS Raw Energy-Volume data.
- `summary_results_vasp.csv`: VASP Fitted EOS parameters (E0, V0, B0), Elastic Moduli (Kv, Gv), and Stiffness Tensor (Cij).
- `summary_results_lammps.csv`: LAMMPS Fitted EOS parameters (E0, V0, B0), Elastic Moduli (Kv, Gv), and Stiffness Tensor (Cij).
- `plot_vasp_{struct_id}.png`: Individual VASP EOS curves for each structure.
- `plot_lammps_{struct_id}.png`: Individual LAMMPS EOS curves for each structure.
- `plot_comparison_{struct_id}.png`: VASP vs LAMMPS comparison overlay for each structure.
