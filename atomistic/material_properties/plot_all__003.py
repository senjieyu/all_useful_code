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

def vinet_eos(V, E0, V0, B0, Bp):
    """
    Vinet (Universal) EOS.
    """
    # x = (V / V0) ** (1/3)
    # eta = 1.5 * (Bp - 1)
    # E = E0 + (9 * B0 * V0 / (eta**2)) * (1 + (eta * (1 - x) - 1) * np.exp(eta * (1 - x)))
    
    # To avoid division by zero if Bp=1, though unlikely for solids (Bp ~ 3-5)
    # And ensure V/V0 is positive
    
    V_safe = np.maximum(V, 1e-8) # Avoid zero or negative volume
    x = (V_safe / V0) ** (1.0/3.0)
    eta = 1.5 * (Bp - 1)
    
    term = 1 + (eta * (1 - x) - 1) * np.exp(eta * (1 - x))
    return E0 + (9 * B0 * V0 / (eta**2)) * term

def fit_eos(volumes, energies):
    """
    Fits E-V data to Vinet EOS.
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
        popt, pcov = curve_fit(vinet_eos, V, E, p0=[E0_guess, V0_guess, B0_guess, Bp_guess], maxfev=10000)
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
            
        # Fallback if natoms is 1 but structure type implies more
        if natoms == 1:
            if stype == "FCC": natoms = 4
            elif stype == "BCC": natoms = 2
            elif stype == "HCP": natoms = 2
            
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

def get_lattice_constant(vol_atom, struct_type):
    """
    Converts volume per atom to lattice constant based on structure type.
    """
    stype = struct_type.upper()
    if stype in ["FCC"]:
        # V = a^3 / 4  => a = (4V)^(1/3)
        return (4 * vol_atom) ** (1.0/3.0)
    elif stype in ["BCC"]:
        # V = a^3 / 2 => a = (2V)^(1/3)
        return (2 * vol_atom) ** (1.0/3.0)
    elif stype in ["SC", "CUBIC"]:
        # V = a^3
        return vol_atom ** (1.0/3.0)
    elif stype in ["DIAMOND"]:
        # V = a^3 / 8
        return (8 * vol_atom) ** (1.0/3.0)
    elif stype in ["HCP"]:
        # Ideal HCP: V_pa = a^3 / sqrt(2) => a = (sqrt(2) * V)^(1/3)
        # Note: This assumes ideal c/a ratio.
        return (np.sqrt(2) * vol_atom) ** (1.0/3.0)
    else:
        # Default fallback to cubic
        return vol_atom ** (1.0/3.0)

def plot_eos(data_list, title, filename, compare=False, x_range=None):
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
        stype = entry.get("struct_type", "Unknown")
        eos = entry["eos_data"]
        natoms = entry.get("natoms", 1)
        
        # Get raw data
        V_total = np.array([d['volume'] for d in eos])
        E_total = np.array([d['energy'] for d in eos])
        
        # Normalize to per-atom
        V_pa = V_total / natoms
        E_pa = E_total / natoms
        
        # Convert V to a
        a_vals = get_lattice_constant(V_pa, stype)
        
        # Shift Energy
        E_ref_pa = np.min(E_pa)
        E_shifted_pa = E_pa - E_ref_pa
        
        # Filter by x_range if provided
        if x_range:
            mask = (a_vals >= x_range[0]) & (a_vals <= x_range[1])
            # If no points in range, skip plotting data but might still plot fit if relevant? 
            # Better to just filter what we plot.
            if not np.any(mask):
                continue
            
            # Apply mask for plotting
            plot_a = a_vals[mask]
            plot_E = E_shifted_pa[mask]
        else:
            plot_a = a_vals
            plot_E = E_shifted_pa
        
        label = f"{method}"
        
        marker = markers.get(method, default_markers[i % len(default_markers)])
        color = colors.get(method, None)
        
        # Scatter points (Lattice Constant)
        plt.scatter(plot_a, plot_E, label=label, marker=marker, color=color)
        
        # Mark minimum point (Always mark the global minimum of the dataset, 
        # but check if it's within range to plot the marker)
        min_idx = np.argmin(E_shifted_pa)
        min_a = a_vals[min_idx]
        min_E = E_shifted_pa[min_idx]
        
        # Only plot min marker if it's within range (or if no range specified)
        if x_range is None or (x_range[0] <= min_a <= x_range[1]):
            min_marker = min_markers.get(method, default_min_markers[i % len(default_min_markers)])
            
            # Plot min point with distinct shape and annotation
            plt.scatter(min_a, min_E, marker=min_marker, s=250, color='gold', edgecolors=color, linewidth=1.5, zorder=10)
            
            # Dynamic offset to avoid overlap
            # Offsets are (x, y) in points
            offsets = [
                (0, 15),    # Top
                (0, -25),   # Bottom
                (25, 0),    # Right
                (-25, 0),   # Left
                (20, 20),   # Top-Right
                (-20, 20)   # Top-Left
            ]
            offset = offsets[i % len(offsets)]
            
            # Adjust alignment based on offset
            ha = 'center'
            if offset[0] > 10: ha = 'left'
            elif offset[0] < -10: ha = 'right'
                
            va = 'bottom'
            if offset[1] < 0: va = 'top'
            elif offset[1] > 10: va = 'bottom'
            else: va = 'center'
            
            plt.annotate(f"{min_a:.3f}", (min_a, min_E), xytext=offset, textcoords='offset points', 
                        ha=ha, va=va, fontsize=9, color=color, fontweight='bold',
                        arrowprops=dict(arrowstyle='-', color=color, lw=0.5))

        # Plot fit line if available
        if entry["fit_params"]:
            # Generate fit curve over the plotted range
            if x_range:
                # Extend slightly beyond range for smooth look
                plot_min_a, plot_max_a = x_range
            else:
                plot_min_a, plot_max_a = min(a_vals), max(a_vals)
                
            # We need to convert a range back to V to generate the curve
            # This is tricky because V->a depends on structure. 
            # We can generate V range corresponding to a range.
            
            # a = (C * V)^(1/3) => a^3 = C * V => V = a^3 / C
            # Factor C depends on structure.
            # get_lattice_constant logic:
            # FCC: a = (4V)^(1/3) => V = a^3 / 4
            # BCC: a = (2V)^(1/3) => V = a^3 / 2
            # SC: a = V^(1/3) => V = a^3
            # Diamond: a = (8V)^(1/3) => V = a^3 / 8
            # HCP: a = (sqrt(2)V)^(1/3) => V = a^3 / sqrt(2)
            
            def get_vol_from_a(a, stype):
                stype = stype.upper()
                if stype == "FCC": return (a**3)/4.0
                elif stype == "BCC": return (a**3)/2.0
                elif stype == "SC": return a**3
                elif stype == "DIAMOND": return (a**3)/8.0
                elif stype == "HCP": return (a**3)/np.sqrt(2)
                else: return a**3
            
            # Generate points in a-space
            fit_a_vals = np.linspace(plot_min_a, plot_max_a, 100)
            fit_V_pa = get_vol_from_a(fit_a_vals, stype)
            fit_V_total = fit_V_pa * natoms
            
            # Reconstruct popt for vinet_eos (B0 must be in eV/A^3)
            popt = [entry["fit_params"]["E0"], entry["fit_params"]["V0"], 
                    entry["fit_params"]["B0_GPa"]/160.21766208, entry["fit_params"]["Bp"]]
            
            fit_E_total = vinet_eos(fit_V_total, *popt)
            
            # Convert fit curve to per-atom
            fit_E_pa = fit_E_total / natoms
            
            ls = linestyles.get(method, '-')
            plt.plot(fit_a_vals, fit_E_pa - E_ref_pa, linestyle=ls, color=color)
            
    plt.xlabel("Lattice Constant ($\AA$)")
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
            # Get structure type (header)
            stype_header = "Unknown"
            for m in methods.values():
                t = m.get("struct_type", "Unknown")
                if t != "Unknown": 
                    stype_header = t
                    break
            
            f.write(f"{'='*80}\n")
            f.write(f"Structure: {sid} ({stype_header})\n")
            f.write(f"{'='*80}\n\n")
            
            # 1. EOS Parameters
            f.write("--- EOS Parameters (Vinet / Universal) ---\n")
            headers = ["Method", "E0/atom (eV)", "V0/atom (A^3)", "B0 (GPa)", "Bp", "a0 (A)"]
            row_fmt = "{:<15} {:<12} {:<12} {:<12} {:<10} {:<12}\n"
            f.write(row_fmt.format(*headers))
            f.write("-" * 80 + "\n")
            
            for method in methods_order:
                if method in methods:
                    e = methods[method]
                    natoms = e.get("natoms", 1)
                    stype = e.get("struct_type", "Unknown")
                    
                    if e["fit_params"]:
                        fp = e["fit_params"]
                        # Approximate a0
                        v0_atom = fp["V0"]/natoms
                        
                        a0 = get_lattice_constant(v0_atom, stype)

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
        
        # Determine x_range based on structure type
        # Auto-detect reasonable range based on data or structure type
        # Or user can hardcode here.
        
        x_range = None
        
        # Heuristic: Find min and max of optimal lattice constants found so far
        # and add a buffer. Or simply exclude extreme outliers.
        # But user specifically asked for "points within range".
        
        # Let's set a broad reasonable range for Cu
        # struct_0 (BCC) ~ 2.9 A
        # struct_1 (FCC) ~ 3.6 A
        # struct_2 (HCP) ~ 2.5 A (a)
        
        # If we want to filter out the huge expansions (like scale 8.0)
        # We can set a limit like [2.0, 5.5]
        x_range = (2.0, 5.5)
        
        if sorted_entries:
             plot_eos(sorted_entries, f"EOS Comparison - {sid}", os.path.join(results_dir, f"plot_comparison_{sid}.png"), compare=True, x_range=x_range)

if __name__ == "__main__":
    main()
