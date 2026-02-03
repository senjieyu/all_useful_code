import os
import sys
import json
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# ==========================================
# EOS Fitting Function (Birch-Murnaghan 3rd)
# ==========================================
def bm3_eos(V, E0, V0, B0, Bp):
    """
    Birch-Murnaghan 3rd order EOS.
    """
    eta = (V0 / V)**(2/3)
    return E0 + 9*V0*B0/16 * ( (eta - 1)**3 * Bp + (eta - 1)**2 * (6 - 4*eta) )

def fit_eos(volumes, energies):
    """
    Fits E-V data to BM3 EOS.
    Returns params: E0, V0, B0 (eV/A^3), Bp
    """
    # Initial guess
    # Sort data
    data = sorted(zip(volumes, energies))
    V = np.array([x[0] for x in data])
    E = np.array([x[1] for x in data])
    
    min_idx = np.argmin(E)
    V0_guess = V[min_idx]
    E0_guess = E[min_idx]
    B0_guess = 1.0 # eV/A^3 ~ 160 GPa
    Bp_guess = 4.0
    
    try:
        popt, pcov = curve_fit(bm3_eos, V, E, p0=[E0_guess, V0_guess, B0_guess, Bp_guess], maxfev=5000)
        return popt
    except Exception as e:
        print(f"Fitting failed: {e}")
        return None

# ------------------------------
# Helpers for elastic properties
# ------------------------------
def get_val(entry, key):
    """
    Retrieve a scalar elastic property from an entry.
    Supports direct Cij components and elastic moduli stored in elastic_data.
    """
    if entry is None:
        return None
    # Elastic moduli from elastic_data if present
    if key in ["Kv", "Gv"]:
        ed = entry.get("elastic_data") or {}
        val = ed.get(f"{key}_GPa")
        return val
    # Cij from 6x6 stiffness matrix
    C = entry.get("Cij")
    if C is None:
        return None
    C = np.array(C)
    mapper = {
        "C11": (0, 0), "C22": (1, 1), "C33": (2, 2),
        "C12": (0, 1), "C13": (0, 2), "C23": (1, 2),
        "C44": (3, 3), "C55": (4, 4), "C66": (5, 5),
    }
    if key in mapper:
        r, c = mapper[key]
        return C[r, c]
    return None

def collect_data(method_dir, method_name):
    """
    Collects results from struct_* folders.
    Returns a list of dicts.
    """
    data_list = []
    # Find struct_* folders
    folders = sorted(glob.glob(os.path.join(method_dir, "struct_*")))
    
    for folder in folders:
        struct_id = os.path.basename(folder)
        
        # EOS
        eos_file = os.path.join(folder, "results_eos.json")
        eos_data = []
        if os.path.exists(eos_file):
            with open(eos_file, 'r') as f:
                eos_data = json.load(f)
                
        # Elastic
        elas_file = os.path.join(folder, "results_elastic.json")
        elas_data = {}
        if os.path.exists(elas_file):
            with open(elas_file, 'r') as f:
                elas_data = json.load(f)
        
        # Determine number of atoms
        # We need this for plotting per-atom quantities
        # Try to read POSCAR or start.lmp to get N
        natoms = 1
        try:
            # Check for POSCAR (VASP)
            poscar = os.path.join(folder, "POSCAR")
            if os.path.exists(poscar):
                from ase.io import read
                atoms = read(poscar)
                natoms = len(atoms)
            else:
                # Check for start.lmp (LAMMPS)
                start_lmp = os.path.join(folder, "start.lmp")
                if os.path.exists(start_lmp):
                    from ase.io import read
                    atoms = read(start_lmp, format='lammps-data', style='atomic')
                    natoms = len(atoms)
        except:
            print(f"Warning: Could not determine natoms for {struct_id}, assuming 1.")
            pass
            
        if not eos_data:
            print(f"Warning: No EOS data for {method_name} {struct_id}")
            continue
            
        # Fit EOS (Using per-atom values if desired? No, usually fit total, then divide)
        # But user wants plot in per-atom units.
        # Let's fit total energy first, then convert for plotting.
        vols = [d['volume'] for d in eos_data]
        ens = [d['energy'] for d in eos_data]
        popt = fit_eos(vols, ens)
        
        entry = {
            "struct_id": struct_id,
            "method": method_name,
            "natoms": natoms,
            "eos_data": eos_data, # List of {volume, energy}
            "elastic_data": elas_data,
            "fit_params": None
        }
        
        if popt is not None:
            entry["fit_params"] = {
                "E0": popt[0],
                "V0": popt[1],
                "B0_GPa": popt[2] * 160.21766208,
                "Bp": popt[3]
            }
            
        # Extract Cij if available
        if elas_data:
            cij = elas_data.get("Cij_GPa_Voigt") or elas_data.get("Cij_GPa")
            if cij:
                entry["Cij"] = cij
            else:
                entry["Cij"] = None
        else:
            entry["Cij"] = None
            
        data_list.append(entry)
        
    return data_list

def plot_eos(data_list, title, filename, compare=False):
    plt.figure(figsize=(8, 6))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(data_list)))
    
    for i, entry in enumerate(data_list):
        method = entry["method"]
        sid = entry["struct_id"]
        eos = entry["eos_data"]
        natoms = entry.get("natoms", 1)
        
        # Get raw data
        V_total = np.array([d['volume'] for d in eos])
        E_total = np.array([d['energy'] for d in eos])
        
        # Normalize to per-atom
        V_pa = V_total / natoms
        E_pa = E_total / natoms
        
        # Shift Energy
        # User requirement: minimum data point value as 0
        # We shift the per-atom energy so min is 0
        E_ref_pa = np.min(E_pa)
        E_shifted_pa = E_pa - E_ref_pa
        
        label = f"{method} {sid}"
        
        # Scatter points (Per Atom)
        plt.scatter(V_pa, E_shifted_pa, label=label, marker='o' if method=="VASP" else 'x')
        
        # Mark minimum point
        min_idx = np.argmin(E_shifted_pa)
        # VASP uses triangle, LAMMPS uses star
        min_marker = '^' if method == "VASP" else '*'
        plt.scatter(V_pa[min_idx], E_shifted_pa[min_idx], marker=min_marker, s=200, color='gold', edgecolors='black', zorder=10)
        
        # Plot fit line if available
        if entry["fit_params"]:
            # Fit was done on TOTAL Energy and Volume
            # fit_params: E0, V0, B0, Bp (Total quantities)
            
            # We need to generate curve for TOTAL, then normalize
            
            fit_V_total = np.linspace(min(V_total), max(V_total), 100)
            
            # Reconstruct popt for bm3_eos (B0 must be in eV/A^3)
            popt = [entry["fit_params"]["E0"], entry["fit_params"]["V0"], 
                    entry["fit_params"]["B0_GPa"]/160.21766208, entry["fit_params"]["Bp"]]
            
            fit_E_total = bm3_eos(fit_V_total, *popt)
            
            # Convert fit curve to per-atom
            fit_V_pa = fit_V_total / natoms
            fit_E_pa = fit_E_total / natoms
            
            # Shift fit curve by the same reference as data points
            # Note: E_ref_pa = min(E_data_pa)
            # So fit curve min might not be exactly 0, but data min is 0.
            
            plt.plot(fit_V_pa, fit_E_pa - E_ref_pa, linestyle='--' if method=="VASP" else '-')
            
    plt.xlabel("Volume/Atom ($\AA^3$)")
    plt.ylabel("$\Delta$E/Atom (eV)")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Saved {filename}")

def save_summary_csv(data_list, filename):
    rows = []
    for entry in data_list:
        row = {
            "Struct": entry["struct_id"],
            "Method": entry["method"],
        }
        if entry["fit_params"]:
            row.update(entry["fit_params"]) # E0, V0, B0_GPa, Bp
            
        if entry["elastic_data"]:
            row["Kv_GPa"] = entry["elastic_data"].get("Kv_GPa")
            row["Gv_GPa"] = entry["elastic_data"].get("Gv_GPa")
        
        if entry["Cij"]:
            cij = np.array(entry["Cij"])
            # Flatten Cij (6x6) -> C11, C12... C66
            for r in range(6):
                for c in range(6):
                    row[f"C{r+1}{c+1}"] = cij[r, c]
                    
        rows.append(row)
    
    if rows:
        pd.DataFrame(rows).to_csv(filename, index=False)
        print(f"Saved {filename}")
    else:
        print(f"No data to save for {filename}")

def save_analysis_report(data_list, filename):
    # Group by struct
    struct_map = {}
    for entry in data_list:
        sid = entry["struct_id"]
        if sid not in struct_map: struct_map[sid] = {}
        struct_map[sid][entry["method"]] = entry
        
    with open(filename, 'w') as f:
        for sid, methods in sorted(struct_map.items()):
            f.write(f"{'='*60}\n")
            f.write(f"Structure: {sid}\n")
            f.write(f"{'='*60}\n\n")
            
            # 1. EOS Parameters
            f.write("--- EOS Parameters (Birch-Murnaghan 3rd) ---\n")
            headers = ["Method", "E0/atom (eV)", "V0/atom (A^3)", "B0 (GPa)", "Bp", "a0 (A)"]
            row_fmt = "{:<10} {:<12} {:<12} {:<12} {:<10} {:<12}\n"
            f.write(row_fmt.format(*headers))
            f.write("-" * 75 + "\n")
            
            for method in ["VASP", "LAMMPS"]:
                if method in methods:
                    e = methods[method]
                    natoms = e.get("natoms", 1)
                    if e["fit_params"]:
                        fp = e["fit_params"]
                        # Approximate a0 assuming cubic
                        a0 = fp["V0"]**(1.0/3.0)
                        f.write(row_fmt.format(
                            method,
                            f"{fp['E0']/natoms:.6f}",
                            f"{fp['V0']/natoms:.4f}",
                            f"{fp['B0_GPa']:.2f}",
                            f"{fp['Bp']:.4f}",
                            f"{a0:.4f}"
                        ))
                    else:
                        f.write(f"{method:<10} Fit failed or no data\n")
            f.write("\n")
            
            # 2. Elastic Moduli & Cij
            f.write("--- Elastic Properties ---\n")
            
            # Print Matrices first
            for method in ["VASP", "LAMMPS"]:
                if method in methods and methods[method]["Cij"]:
                    f.write(f"[{method} Stiffness Tensor Cij (GPa)]:\n")
                    C = np.array(methods[method]["Cij"])
                    
                    # Clean up small values for display
                    C_display = C.copy()
                    C_display[np.abs(C_display) < 1] = 0.0
                    
                    # Print 6x6
                    for i in range(6):
                        f.write("  " + " ".join([f"{x:8.2f}" for x in C_display[i]]) + "\n")
                    f.write("\n")
            
            # Comparison Table
            if "VASP" in methods and "LAMMPS" in methods:
                v_data = methods["VASP"]
                l_data = methods["LAMMPS"]
                
                f.write("--- Comparison (VASP vs LAMMPS) ---\n")
                f.write(f"{'Property':<20} {'Unit':<6} {'VASP':<10} {'LAMMPS':<10} {'Diff':<10}\n")
                f.write("-" * 60 + "\n")
                
                # Calculate derived properties for LAMMPS
                if "LAMMPS" in methods and methods["LAMMPS"]["Cij"]:
                    C = np.array(methods["LAMMPS"]["Cij"])
                    # Assuming cubic for simplicity, but general formulas exist
                    # Reuss (lower bound)
                    # S = inv(C)
                    try:
                        S = np.linalg.inv(C)
                        Kr = 1.0 / (3.0 * (S[0,0] + 2*S[0,1])) # Approximate for cubic
                        # Gr = 15 / (4*(S[0,0]-S[0,1]) + 3*S[3,3]) # Approximate for cubic
                        # Use general Voigt-Reuss-Hill
                        # Voigt (already calculated in entry["elastic_data"] usually)
                        
                        # Let's use the stored Kv/Gv if available, else calc
                        l_entry = methods["LAMMPS"]
                        
                        # Recalculate full VRH set for table
                        # Voigt
                        Kv_v = (C[0,0] + C[1,1] + C[2,2] + 2*(C[0,1] + C[1,2] + C[2,0])) / 9.0
                        Gv_v = ((C[0,0] + C[1,1] + C[2,2]) - (C[0,1] + C[1,2] + C[2,0]) + 3*(C[3,3] + C[4,4] + C[5,5])) / 15.0
                        
                        # Reuss
                        Kr_v = 1.0 / ((S[0,0] + S[1,1] + S[2,2]) + 2*(S[0,1] + S[1,2] + S[2,0]))
                        Gr_v = 15.0 / (4*(S[0,0] + S[1,1] + S[2,2]) - 4*(S[0,1] + S[1,2] + S[2,0]) + 3*(S[3,3] + S[4,4] + S[5,5]))
                        
                        K_vrh = (Kv_v + Kr_v) / 2.0
                        G_vrh = (Gv_v + Gr_v) / 2.0
                        
                        # Update l_entry with these derived values for display
                        # We store them in a temporary dict for the loop
                        l_derived = {
                            "Kv": Kv_v, "Gv": Gv_v,
                            "Kr": Kr_v, "Gr": Gr_v,
                            "K_vrh": K_vrh, "G_vrh": G_vrh,
                            "Anisotropy": 2 * C[3,3] / (C[0,0] - C[0,1]) if (C[0,0] - C[0,1]) != 0 else 0,
                            "Poisson": (3*K_vrh - 2*G_vrh) / (2*(3*K_vrh + G_vrh))
                        }
                    except:
                        l_derived = {}
                else:
                    l_derived = {}

                # Do same for VASP if needed (omitted for brevity, assuming VASP has them or we compute on fly)
                
                def get_val_enhanced(entry, key, derived=None):
                    # Try derived first
                    if derived and key in derived:
                        return derived[key]
                    # Try standard lookup
                    return get_val(entry, key)

                compare_list = [
                    ("C11", "GPa"), ("C22", "GPa"), ("C33", "GPa"),
                    ("C12", "GPa"), ("C13", "GPa"), ("C23", "GPa"),
                    ("C44", "GPa"), ("C55", "GPa"), ("C66", "GPa"),
                    ("Kv", "GPa"), ("Gv", "GPa"),
                    ("Kr", "GPa"), ("Gr", "GPa"),
                    ("K_vrh", "GPa"), ("G_vrh", "GPa"),
                    ("Anisotropy", ""), ("Poisson", "")
                ]
                
                f.write(f"{'Property':<20} {'Unit':<6} {'VASP':<10} {'LAMMPS':<10} {'Diff %':<10}\n")
                f.write("-" * 65 + "\n")

                for key, unit in compare_list:
                    # Calculate VASP derived
                    v_derived = {}
                    if "VASP" in methods and methods["VASP"]["Cij"]:
                         try:
                             C = np.array(methods["VASP"]["Cij"])
                             S = np.linalg.inv(C)
                             Kv_v = (C[0,0] + C[1,1] + C[2,2] + 2*(C[0,1] + C[1,2] + C[2,0])) / 9.0
                             Gv_v = ((C[0,0] + C[1,1] + C[2,2]) - (C[0,1] + C[1,2] + C[2,0]) + 3*(C[3,3] + C[4,4] + C[5,5])) / 15.0
                             Kr_v = 1.0 / ((S[0,0] + S[1,1] + S[2,2]) + 2*(S[0,1] + S[1,2] + S[2,0]))
                             Gr_v = 15.0 / (4*(S[0,0] + S[1,1] + S[2,2]) - 4*(S[0,1] + S[1,2] + S[2,0]) + 3*(S[3,3] + S[4,4] + S[5,5]))
                             K_vrh = (Kv_v + Kr_v) / 2.0
                             G_vrh = (Gv_v + Gr_v) / 2.0
                             v_derived = {
                                "Kv": Kv_v, "Gv": Gv_v, "Kr": Kr_v, "Gr": Gr_v, "K_vrh": K_vrh, "G_vrh": G_vrh,
                                "Anisotropy": 2 * C[3,3] / (C[0,0] - C[0,1]) if (C[0,0] - C[0,1]) != 0 else 0,
                                "Poisson": (3*K_vrh - 2*G_vrh) / (2*(3*K_vrh + G_vrh))
                             }
                         except: pass

                    v_val = get_val_enhanced(v_data, key, v_derived)
                    l_val = get_val_enhanced(l_data, key, l_derived)
                    
                    v_str = f"{v_val:.2f}" if v_val is not None else "N/A"
                    l_str = f"{l_val:.2f}" if l_val is not None else "N/A"
                    
                    # Percent Diff: (LAMMPS - VASP) / VASP * 100
                    if v_val is not None and l_val is not None and v_val != 0:
                        diff_val = (l_val - v_val) / abs(v_val) * 100
                        diff_str = f"{diff_val:+.1f}%"
                    else:
                        diff_str = "-"
                    
                    f.write(f"{key:<20} {unit:<6} {v_str:<10} {l_str:<10} {diff_str:<10}\n")
                f.write("\n")
    
    print(f"Saved {filename}")

def main():
    # Paths
    base_dir = os.path.dirname(os.path.abspath(__file__)) # Unified_Workflow/Analysis
    root_dir = os.path.dirname(base_dir) # Unified_Workflow
    
    vasp_dir = os.path.join(root_dir, "VASP")
    lammps_dir = os.path.join(root_dir, "LAMMPS")
    
    # Create results directory
    results_dir = os.path.join(root_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    
    # Collect
    print("Collecting VASP data...")
    vasp_data = collect_data(vasp_dir, "VASP")
    print("Collecting LAMMPS data...")
    lammps_data = collect_data(lammps_dir, "LAMMPS")
    
    # Output Raw Data CSV
    def save_raw_csv(data_list, filename):
        rows = []
        for entry in data_list:
            for p in entry["eos_data"]:
                rows.append({
                    "Struct": entry["struct_id"],
                    "Method": entry["method"],
                    "Scale": p["scale"],
                    "Volume": p["volume"],
                    "Energy": p["energy"]
                })
        if rows:
            pd.DataFrame(rows).to_csv(filename, index=False)
            print(f"Saved {filename}")
        else:
            print(f"No data to save for {filename}")

    # Output Raw Data CSVs (Split by method)
    save_raw_csv(vasp_data, os.path.join(results_dir, "eos_raw_data_vasp.csv"))
    save_raw_csv(lammps_data, os.path.join(results_dir, "eos_raw_data_lammps.csv"))
    
    # Output Summary CSVs (Split by method)
    save_summary_csv(vasp_data, os.path.join(results_dir, "summary_results_vasp.csv"))
    save_summary_csv(lammps_data, os.path.join(results_dir, "summary_results_lammps.csv"))
    
    # Generate Analysis Report (TXT)
    save_analysis_report(vasp_data + lammps_data, os.path.join(results_dir, "analysis_report.txt"))
    
    # Organize data by structure
    struct_map = {}
    all_entries = vasp_data + lammps_data
    for entry in all_entries:
        sid = entry["struct_id"]
        if sid not in struct_map:
            struct_map[sid] = {"VASP": None, "LAMMPS": None}
        struct_map[sid][entry["method"]] = entry
        
    # Generate Plots per Structure
    for sid, methods in sorted(struct_map.items()):
        v_entry = methods["VASP"]
        l_entry = methods["LAMMPS"]
        
        # 1. VASP Only
        if v_entry:
            plot_eos([v_entry], f"VASP EOS - {sid}", os.path.join(results_dir, f"plot_vasp_{sid}.png"))
            
        # 2. LAMMPS Only
        if l_entry:
            plot_eos([l_entry], f"LAMMPS EOS - {sid}", os.path.join(results_dir, f"plot_lammps_{sid}.png"))
            
        # 3. Comparison
        if v_entry and l_entry:
            plot_eos([v_entry, l_entry], f"VASP vs LAMMPS - {sid}", os.path.join(results_dir, f"plot_comparison_{sid}.png"), compare=True)

if __name__ == "__main__":
    main()
