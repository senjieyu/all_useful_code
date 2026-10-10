import os
import sys
import json
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# ==========================================
# Helpers
# ==========================================
def load_style_info(root_dir):
    """
    Reads style.txt from root_dir to get structure types.
    Returns dict: {int(index): str(type)}
    """
    style_map = {}
    style_path = os.path.join(root_dir, "style.txt")
    if os.path.exists(style_path):
        try:
            with open(style_path, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            idx = int(parts[0])
                            stype = parts[1].upper()
                            style_map[idx] = stype
                        except:
                            continue
        except Exception as e:
            print(f"Error reading style.txt: {e}")
    return style_map

def get_structure_type(struct_id, style_map):
    """
    Parses struct_id (e.g. "struct_0") and returns type from map.
    """
    try:
        idx = int(struct_id.split('_')[1])
        return style_map.get(idx, "Unknown")
    except:
        return "Unknown"

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
    
    if len(V) < 4:
        return None

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

def collect_data(method_dir, method_name, style_map):
    """
    Collects results from struct_* folders.
    Returns a list of dicts.
    """
    data_list = []
    if not os.path.exists(method_dir):
        print(f"Warning: Directory {method_dir} does not exist.")
        return data_list
        
    # Find struct_* folders
    folders = sorted(glob.glob(os.path.join(method_dir, "struct_*")))
    
    for folder in folders:
        struct_id = os.path.basename(folder)
        stype = get_structure_type(struct_id, style_map)
        
        # EOS
        eos_file = os.path.join(folder, "results_eos.json")
        eos_data = []
        if os.path.exists(eos_file):
            try:
                with open(eos_file, 'r') as f:
                    eos_data = json.load(f)
            except: pass
                
        # Elastic
        elas_file = os.path.join(folder, "results_elastic.json")
        elas_data = {}
        if os.path.exists(elas_file):
            try:
                with open(elas_file, 'r') as f:
                    elas_data = json.load(f)
            except: pass
        
        # Determine number of atoms
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
            # print(f"Warning: Could not determine natoms for {struct_id}, assuming 1.")
            pass
            
        if not eos_data:
            # print(f"Warning: No EOS data for {method_name} {struct_id}")
            continue
            
        # Fit EOS
        vols = [d['volume'] for d in eos_data]
        ens = [d['energy'] for d in eos_data]
        popt = fit_eos(vols, ens)
        
        entry = {
            "struct_id": struct_id,
            "struct_type": stype,
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
    
    # Use distinct markers/styles
    markers = {'VASP': 'o', 'LAMMPS_SUS2': 's', 'LAMMPS_EAM': '^'}
    min_markers = {'VASP': '*', 'LAMMPS_SUS2': 'D', 'LAMMPS_EAM': 'P'}
    linestyles = {'VASP': '-', 'LAMMPS_SUS2': '--', 'LAMMPS_EAM': '-.'}
    colors = {'VASP': 'black', 'LAMMPS_SUS2': 'blue', 'LAMMPS_EAM': 'red'}
    
    # Fallback for unknown methods
    default_markers = ['o', 's', '^', 'D', 'x']
    default_min_markers = ['*', 'D', 'P', 'X', 'h']
    
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
        E_ref_pa = np.min(E_pa)
        E_shifted_pa = E_pa - E_ref_pa
        
        label = f"{method}"
        
        marker = markers.get(method, default_markers[i % len(default_markers)])
        color = colors.get(method, None)
        
        # Scatter points (Per Atom)
        plt.scatter(V_pa, E_shifted_pa, label=label, marker=marker, color=color)
        
        # Mark minimum point
        min_idx = np.argmin(E_shifted_pa)
        min_V = V_pa[min_idx]
        min_E = E_shifted_pa[min_idx]
        
        min_marker = min_markers.get(method, default_min_markers[i % len(default_min_markers)])
        
        # Plot min point with distinct shape and annotation
        plt.scatter(min_V, min_E, marker=min_marker, s=250, color='gold', edgecolors=color, linewidth=1.5, zorder=10)
        plt.annotate(f"{min_V:.2f}", (min_V, min_E), xytext=(0, 15), textcoords='offset points', ha='center', fontsize=9, color=color, fontweight='bold')

        # Plot fit line if available
        if entry["fit_params"]:
            fit_V_total = np.linspace(min(V_total), max(V_total), 100)
            
            # Reconstruct popt for bm3_eos (B0 must be in eV/A^3)
            popt = [entry["fit_params"]["E0"], entry["fit_params"]["V0"], 
                    entry["fit_params"]["B0_GPa"]/160.21766208, entry["fit_params"]["Bp"]]
            
            fit_E_total = bm3_eos(fit_V_total, *popt)
            
            # Convert fit curve to per-atom
            fit_V_pa = fit_V_total / natoms
            fit_E_pa = fit_E_total / natoms
            
            ls = linestyles.get(method, '-')
            plt.plot(fit_V_pa, fit_E_pa - E_ref_pa, linestyle=ls, color=color)
            
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
            "Type": entry.get("struct_type", "Unknown"),
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

def save_analysis_report(data_list, filename, methods_order):
    # Group by struct
    struct_map = {}
    for entry in data_list:
        sid = entry["struct_id"]
        if sid not in struct_map: struct_map[sid] = {}
        struct_map[sid][entry["method"]] = entry
        
    with open(filename, 'w') as f:
        for sid, methods in sorted(struct_map.items()):
            # Get structure type
            stype = "Unknown"
            for m in methods.values():
                stype = m.get("struct_type", "Unknown")
                if stype != "Unknown": break
            
            f.write(f"{'='*80}\n")
            f.write(f"Structure: {sid} ({stype})\n")
            f.write(f"{'='*80}\n\n")
            
            # 1. EOS Parameters
            f.write("--- EOS Parameters (Birch-Murnaghan 3rd) ---\n")
            headers = ["Method", "E0/atom (eV)", "V0/atom (A^3)", "B0 (GPa)", "Bp", "a0 (A)"]
            row_fmt = "{:<15} {:<12} {:<12} {:<12} {:<10} {:<12}\n"
            f.write(row_fmt.format(*headers))
            f.write("-" * 80 + "\n")
            
            for method in methods_order:
                if method in methods:
                    e = methods[method]
                    natoms = e.get("natoms", 1)
                    if e["fit_params"]:
                        fp = e["fit_params"]
                        # Approximate a0
                        v0_atom = fp["V0"]/natoms
                        a0 = 0.0
                        if stype in ["FCC"]:
                             a0 = (4 * v0_atom) ** (1.0/3.0)
                        elif stype in ["BCC"]:
                             a0 = (2 * v0_atom) ** (1.0/3.0)
                        elif stype in ["SC", "CUBIC"]:
                             a0 = v0_atom ** (1.0/3.0)
                        elif stype in ["HCP"]:
                             # Assuming ideal c/a for approximation
                             a0 = (2 * v0_atom / (np.sqrt(3) * 1.633)) ** (1.0/3.0)
                        else:
                             a0 = v0_atom ** (1.0/3.0)

                        f.write(row_fmt.format(
                            method,
                            f"{fp['E0']/natoms:.6f}",
                            f"{v0_atom:.4f}",
                            f"{fp['B0_GPa']:.2f}",
                            f"{fp['Bp']:.4f}",
                            f"{a0:.4f}"
                        ))
                    else:
                        f.write(f"{method:<15} Fit failed or no data\n")
            f.write("\n")
            
            # 2. Elastic Moduli & Cij
            f.write("--- Elastic Properties ---\n")
            
            # Comparison Table
            
            # Calculate derived properties
            def calc_derived(entry):
                if not entry or not entry.get("Cij"): return {}
                C = np.array(entry["Cij"])
                derived = {}
                try:
                    # Voigt (General)
                    Kv_v = (C[0,0] + C[1,1] + C[2,2] + 2*(C[0,1] + C[1,2] + C[2,0])) / 9.0
                    Gv_v = ((C[0,0] + C[1,1] + C[2,2]) - (C[0,1] + C[1,2] + C[2,0]) + 3*(C[3,3] + C[4,4] + C[5,5])) / 15.0
                    
                    derived.update({
                        "Kv": Kv_v, "Gv": Gv_v,
                    })
                    
                    if stype in ["FCC", "BCC", "SC", "CUBIC"]:
                        # Zener Anisotropy
                        A = 2 * C[3,3] / (C[0,0] - C[0,1]) if (C[0,0] - C[0,1]) != 0 else 0
                        derived["Anisotropy"] = A
                except: pass
                return derived

            derived_maps = {}
            for method in methods_order:
                if method in methods:
                    derived_maps[method] = calc_derived(methods[method])

            def get_val_enhanced(entry, key, derived=None):
                if derived and key in derived: return derived[key]
                return get_val(entry, key)

            compare_list = [
                ("C11", "GPa"), ("C12", "GPa"), ("C44", "GPa"),
                ("Kv", "GPa"), ("Gv", "GPa"),
                ("Anisotropy", "")
            ]
            
            # Build Header
            header = f"{'Property':<20} {'Unit':<6}"
            for method in methods_order:
                header += f" {method:<15}"
            f.write(header + "\n")
            f.write("-" * len(header) + "\n")

            for key, unit in compare_list:
                line = f"{key:<20} {unit:<6}"
                has_data = False
                for method in methods_order:
                    val_str = "N/A"
                    if method in methods:
                        val = get_val_enhanced(methods[method], key, derived_maps.get(method))
                        if val is not None:
                            val_str = f"{val:.2f}"
                            has_data = True
                    line += f" {val_str:<15}"
                
                if has_data:
                    f.write(line + "\n")
            f.write("\n")
    
    print(f"Saved {filename}")

def main():
    # Paths
    base_dir = os.path.dirname(os.path.abspath(__file__)) # Unified_Workflow/Analysis
    root_dir = os.path.dirname(base_dir) # Unified_Workflow
    
    vasp_dir = os.path.join(root_dir, "VASP")
    lammps_sus2_dir = os.path.join(root_dir, "LAMMPS_SUS2")
    lammps_eam_dir = os.path.join(root_dir, "LAMMPS_EAM")
    
    # Create results directory
    results_dir = os.path.join(root_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    
    style_map = load_style_info(root_dir)
    
    # Collect
    print("Collecting VASP data...")
    vasp_data = collect_data(vasp_dir, "VASP", style_map)
    
    print("Collecting LAMMPS_SUS2 data...")
    sus2_data = collect_data(lammps_sus2_dir, "LAMMPS_SUS2", style_map)
    
    print("Collecting LAMMPS_EAM data...")
    eam_data = collect_data(lammps_eam_dir, "LAMMPS_EAM", style_map)
    
    all_data = vasp_data + sus2_data + eam_data
    methods_order = ["VASP", "LAMMPS_SUS2", "LAMMPS_EAM"]
    
    # Output Summary CSVs
    save_summary_csv(vasp_data, os.path.join(results_dir, "summary_results_vasp.csv"))
    save_summary_csv(sus2_data, os.path.join(results_dir, "summary_results_sus2.csv"))
    save_summary_csv(eam_data, os.path.join(results_dir, "summary_results_eam.csv"))
    
    # Generate Analysis Report (TXT)
    save_analysis_report(all_data, os.path.join(results_dir, "analysis_report.txt"), methods_order)
    
    # Organize data by structure
    struct_map = {}
    for entry in all_data:
        sid = entry["struct_id"]
        if sid not in struct_map:
            struct_map[sid] = []
        struct_map[sid].append(entry)
        
    # Generate Plots per Structure
    for sid, entries in sorted(struct_map.items()):
        # Filter entries to sort by method order
        sorted_entries = []
        for m in methods_order:
            for e in entries:
                if e["method"] == m:
                    sorted_entries.append(e)
                    break
        
        if sorted_entries:
             plot_eos(sorted_entries, f"EOS Comparison - {sid}", os.path.join(results_dir, f"plot_comparison_{sid}.png"), compare=True)

if __name__ == "__main__":
    main()
