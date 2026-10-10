import os, subprocess, re
import numpy as np
import matplotlib
matplotlib.use('Agg') # Use non-interactive backend
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.data import atomic_masses, atomic_numbers
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

# Set font to Times New Roman
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

# Configuration
XYZ_PATH = "153.xyz"       # Relative path, assumes running from the same directory
MODEL    = "mlip.ini"      # Using mlip.ini configuration
INFILE   = "in.sus2"
# Command to run LAMMPS. Can be overridden by LMP_BIN environment variable.
LMP_BIN  = os.environ.get("LMP_BIN", "lmp") 
OUTDIR   = "run_153"
ELEMENTS = ["Cu"]          # Define elements present in the system, order must match LAMMPS types

os.makedirs(OUTDIR, exist_ok=True)

print(f"Reading structures from {XYZ_PATH}...")
frames = read(XYZ_PATH, index=":")
print(f"Found {len(frames)} structures.")

import pandas as pd

def get_ref_fields(at):
    # Energy (eV)
    E = None
    # Try getting from calculator first (ASE standard for extxyz)
    try:
        E = at.get_potential_energy()
    except:
        # Fallback to info keys
        for k in ["energy", "E", "dft_energy", "ref_energy", "free_energy"]:
            if k in at.info:
                E = float(at.info[k])
                break
            
    # Forces (eV/A)
    F = None
    try:
        F = at.get_forces()
    except:
        for k in ["forces", "force", "dft_forces", "ref_forces"]:
            if k in at.arrays:
                F = np.array(at.arrays[k], float)
                break
            
    # Stress (Voigt: xx,yy,zz,xy,xz,yz)
    S = None
    # Try calculator stress first
    try:
        S = at.get_stress(voigt=False) # 3x3
        if S.shape == (3,3):
             voigt = np.array([S[0,0], S[1,1], S[2,2], S[0,1], S[0,2], S[1,2]], float)
             S = voigt
    except:
        pass

    if S is None:
        if "stress" in at.info:
            s = np.array(at.info["stress"], float)
            if s.shape == (3,3):
                voigt = np.array([s[0,0], s[1,1], s[2,2], s[0,1], s[0,2], s[1,2]], float)
            elif s.size == 9:
                s = s.reshape(3,3)
                voigt = np.array([s[0,0], s[1,1], s[2,2], s[0,1], s[0,2], s[1,2]], float)
            else:
                voigt = s
            S = voigt
        elif "virial" in at.info:
            v = np.array(at.info["virial"], float)
            if v.size == 9: # 3x3 flattened
                 v = v.reshape(3,3)
                 voigt = np.array([v[0,0], v[1,1], v[2,2], v[0,1], v[0,2], v[1,2]], float)
            elif v.size == 6:
                 voigt = v
            else:
                 voigt = None
            
            if voigt is not None:
                 vol = at.get_volume()
                 S = - voigt / vol 

    # Heuristic unit conversion for Stress
    if S is not None:
        if np.nanmedian(np.abs(S)) < 5.0: 
            S = S * 160.21766208
            
    return E, F, S

def write_lammps_data(i, at):
    data_path = os.path.join(OUTDIR, f"data_{i:04d}.lmp")
    # write LAMMPS data. By default ASE might not write masses if not explicitly asked,
    # or LAMMPS input needs them. We'll handle masses in the LAMMPS input script via variables.
    write(data_path, at, format="lammps-data", atom_style="atomic")
    return data_path

# Generate mass commands based on ELEMENTS
mass_commands = ""
for idx, el in enumerate(ELEMENTS, 1):
    mass = atomic_masses[atomic_numbers[el]]
    mass_commands += f"mass {idx} {mass}\n"

# Write a separate file for masses to be included in LAMMPS, or we can just pass them.
# It's cleaner to write a mass file and include it.
mass_file_path = os.path.join(OUTDIR, "masses.in")
with open(mass_file_path, "w") as f:
    f.write(mass_commands)

pred_E, ref_E = [], []
pred_F, ref_F = [], []
pred_S, ref_S = [], []

print("Starting LAMMPS calculations...")
for i, at in enumerate(frames, 1):
    data_path = write_lammps_data(i, at)
    outprefix = os.path.join(OUTDIR, f"conf_{i:04d}")
    
    # Run LAMMPS
    # Note: Ensure lmp is in your PATH or specify full path in LMP_BIN
    cmd = [
        LMP_BIN,
        "-var", "datafile", data_path,
        "-var", "model", MODEL,
        "-var", "outprefix", outprefix,
        "-var", "massfile", mass_file_path,
        "-in", INFILE
    ]
    
    # Check if output exists to skip calculation (Resume mode)
    if os.path.exists(f"{outprefix}_thermo.txt") and os.path.exists(f"{outprefix}_forces.dump"):
        # print(f"Frame {i} output exists, skipping calculation.")
        pass
    else:
        # We use shell=True if LMP_BIN contains spaces (e.g. "mpirun -np 4 lmp")
        # But subprocess.run with list args expects executable as first arg.
        # If LMP_BIN is complex, better to split it or use shell=True with string.
        # Here we split simplistic commands.
        if " " in LMP_BIN:
            cmd = LMP_BIN.split() + cmd[1:]
            
        try:
            # Run LAMMPS and capture output for debugging
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running LAMMPS for frame {i}: {e}")
            print("STDOUT:", e.stdout)
            print("STDERR:", e.stderr)
            continue
        except FileNotFoundError:
            print(f"Error: LAMMPS executable '{LMP_BIN}' not found.")
            break

    # Parse Thermo Output (Energy, Stress diagonal)
    try:
        with open(f"{outprefix}_thermo.txt") as f:
            line = f.read().strip()
        vals = dict(zip(line.split()[::2], line.split()[1::2]))
        Etot = float(vals["ETOT"])
        pxx = float(vals["PXX"])
        pyy = float(vals["PYY"])
        pzz = float(vals["PZZ"])
    except Exception as e:
        print(f"Error parsing thermo for frame {i}: {e}")
        continue

    # Parse Forces
    forces_file = f"{outprefix}_forces.dump"
    fx, fy, fz = [], [], []
    try:
        with open(forces_file) as f:
            for l in f:
                if l.strip().startswith("#") or l.strip().startswith("ITEM:"): 
                    continue
                parts = l.split()
                if len(parts) >= 8:
                    fx.append(float(parts[5]))
                    py = float(parts[6]) # Typo in list logic? No
                    fy.append(float(parts[6]))
                    fz.append(float(parts[7]))
        F_lmp = np.column_stack([fx, fy, fz])
    except Exception as e:
        print(f"Error parsing forces for frame {i}: {e}")
        continue

    # Get Reference
    E_ref, F_ref, S_ref = get_ref_fields(at)

    # Collect Data
    if E_ref is not None:
        pred_E.append(Etot / len(at))     # eV/atom
        ref_E.append(E_ref / len(at))
        
    if F_ref is not None and F_lmp.shape == F_ref.shape:
        pred_F.append(F_lmp.reshape(-1))
        ref_F.append(F_ref.reshape(-1))
        
    if S_ref is not None:
        # LAMMPS pxx, pyy, pzz are in bars usually? 
        # units metal -> pressure is in bars. 
        # Wait, check LAMMPS units metal: pressure = bars.
        # 1 bar = 0.0001 GPa.
        # But 'stress/atom' compute is in pressure*volume units (bar*A^3).
        # We need to divide by volume and convert to GPa.
        # Actually, the thermo output PXX PYY PZZ are global pressure components (in bars).
        # Pressure = (Pxx + Pyy + Pzz)/3.
        # Pxx in thermo is usually in bars for 'units metal'.
        # 10000 bar = 1 GPa.
        # So we need to divide by 10000 to get GPa.
        
        # NOTE: The provided OCR script used PXX directly and compared with Reference.
        # If Reference was converted to GPa (x160), and LAMMPS PXX is in bars...
        # We need to be careful.
        # Let's check 'compute stress/atom'. It produces stress*volume.
        # The script reads PXX from 'thermo'. 'thermo_style custom pxx'.
        # pxx in thermo is Pressure tensor component.
        # For 'units metal', pxx is in bars.
        # To match GPa, we divide by 10000.
        
        # Let's adjust the script to convert LAMMPS PXX (bars) to GPa.
        pxx_gpa = pxx / 10000.0
        pyy_gpa = pyy / 10000.0
        pzz_gpa = pzz / 10000.0
        
        # Note: LAMMPS pressure is -Stress. P = -S.
        # So Sxx = -Pxx.
        # We should invert the sign of PXX from LAMMPS to get Stress.
        
        if S_ref.size >= 3:
            # Taking only diagonal for simplicity as per OCR script
            pred_S.append(np.array([-pxx_gpa, -pyy_gpa, -pzz_gpa])) 
            ref_S.append(np.array([S_ref[0], S_ref[1], S_ref[2]]))
            
    if i % 10 == 0:
        print(f"Processed {i}/{len(frames)}")

print("Processing complete. Generating plots...")

# Convert to arrays
if pred_E:
    pred_E = np.array(pred_E)
    ref_E = np.array(ref_E)
if pred_F:
    pred_F = np.concatenate(pred_F)
    ref_F = np.concatenate(ref_F)
if pred_S:
    pred_S = np.array(pred_S)
    ref_S = np.array(ref_S)

# Plotting
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
plt.subplots_adjust(wspace=0.3, right=0.9, bottom=0.15, top=0.9)

def plot_pane(ax, x, y, xlabel, ylabel, label, unit_factor=1.0, unit_name=""):
    if len(x) == 0: return None
    
    # Range
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    xymin = min(np.min(x), np.min(y))
    # Extend range slightly
    margin = (xymax - xymin) * 0.1
    vmin = xymin - margin
    vmax = xymax + margin
    
    # y=x line
    ax.plot([vmin, vmax], [vmin, vmax], ls="--", c="0.6", lw=1.5)
    
    # Scatter plot (points)
    # Using simple color (e.g. blue) with alpha
    hb = ax.scatter(x, y, alpha=0.6, s=20, c='royalblue', edgecolors='none')
    
    # Calculate stats with unit conversion
    mae = mean_absolute_error(x, y) * unit_factor
    rmse = np.sqrt(mean_squared_error(x, y)) * unit_factor
    r2 = r2_score(x, y)
    
    # Format stats text
    stats_text = (f"MAE: {mae:.3f} {unit_name}\n"
                  f"RMSE: {rmse:.3f} {unit_name}\n"
                  f"R²: {r2:.3f}")
    
    # Annotation box
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='left', bbox=props)
            
    # Subplot label (a), (b), (c) - Bold, top-left inside
    ax.text(0.02, 0.98, label, transform=ax.transAxes, fontsize=14, fontweight='bold', 
            va='top', ha='left')
            
    # Axis labels and ticks
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(labelsize=10, direction='in', which='both', top=True, right=True)
    
    # Set limits to be square-ish
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    
    return hb

hb_last = None

if len(pred_E) > 0:
    # Energy: eV/atom -> meV/atom (x1000)
    hb_last = plot_pane(axes[0], ref_E, pred_E, 
              "DFT Energy (eV/atom)", "SUS Energy (eV/atom)", 
              "(a)", unit_factor=1000.0, unit_name="(meV/atom)")

if len(pred_F) > 0:
    # Force: eV/A -> meV/A (x1000)
    hb_last = plot_pane(axes[1], ref_F, pred_F, 
              "DFT Force (eV/$\AA$)", "SUS Force (eV/$\AA$)", 
              "(b)", unit_factor=1000.0, unit_name="(meV/$\AA$)")

if len(pred_S) > 0:
    # Stress: GPa -> GPa (x1)
    hb_last = plot_pane(axes[2], ref_S.reshape(-1), pred_S.reshape(-1), 
              "DFT Stress (GPa)", "SUS Stress (GPa)", 
              "(c)", unit_factor=1.0, unit_name="(GPa)")

# Add colorbar (only if density was used, but for scatter it's just blue)
# fig.colorbar(hb_last, cax=cbar_ax) # Commented out for scatter

plt.savefig(os.path.join(OUTDIR, "fig_error.png"), dpi=300, bbox_inches='tight')
print(f"Figure saved to {os.path.join(OUTDIR, 'fig_error.png')}")

# Save CSVs
print("Saving CSV data...")
if len(pred_E) > 0:
    df_E = pd.DataFrame({'DFT_Energy_eV_atom': ref_E, 'SUS_Energy_eV_atom': pred_E})
    df_E.to_csv(os.path.join(OUTDIR, "energy_data.csv"), index=False)

if len(pred_F) > 0:
    # F is flat, let's reshape to Nx3 for CSV clarity if possible, or just flat
    # ref_F is (N_atoms*3,)
    # Maybe 3 columns: Fx, Fy, Fz? No, it's concatenated.
    # Simple: Ref, Pred
    df_F = pd.DataFrame({'DFT_Force_eV_A': ref_F, 'SUS_Force_eV_A': pred_F})
    df_F.to_csv(os.path.join(OUTDIR, "force_data.csv"), index=False)

if len(pred_S) > 0:
    # S is (N_frames, 3) flattened? 
    # plot_pane uses reshape(-1), so it's 1D array of all components
    df_S = pd.DataFrame({'DFT_Stress_GPa': ref_S.reshape(-1), 'SUS_Stress_GPa': pred_S.reshape(-1)})
    df_S.to_csv(os.path.join(OUTDIR, "stress_data.csv"), index=False)

print("Data saved to csv files.")

np.savez(os.path.join(OUTDIR, "pred_ref.npz"),
         E_pred=pred_E, E_ref=ref_E,
         F_pred=pred_F, F_ref=ref_F,
         S_pred=pred_S, S_ref=ref_S)
print("Data saved to pred_ref.npz")
