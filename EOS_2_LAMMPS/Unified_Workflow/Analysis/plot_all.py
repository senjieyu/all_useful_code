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
                
        if not eos_data:
            print(f"Warning: No EOS data for {method_name} {struct_id}")
            continue
            
        # Fit EOS
        vols = [d['volume'] for d in eos_data]
        ens = [d['energy'] for d in eos_data]
        popt = fit_eos(vols, ens)
        
        entry = {
            "struct_id": struct_id,
            "method": method_name,
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
        
        V = np.array([d['volume'] for d in eos])
        E = np.array([d['energy'] for d in eos])
        
        # Shift Energy
        # User requirement: minimum value 0 point as curve lowest point
        # If fit exists, use E0. If not, use min(E).
        if entry["fit_params"]:
            E_ref = entry["fit_params"]["E0"]
        else:
            E_ref = np.min(E)
            
        E_shifted = E - E_ref
        
        label = f"{method} {sid}"
        
        # Scatter points
        plt.scatter(V, E_shifted, label=label, marker='o' if method=="VASP" else 'x')
        
        # Plot fit line if available
        if entry["fit_params"]:
            # fit_params["E0"] is the absolute E0.
            # We want to plot (Fit - E_ref).
            # If E_ref == E0, then the minimum of the curve is at 0.
            
            fit_V = np.linspace(min(V), max(V), 100)
            # Reconstruct popt for bm3_eos (B0 must be in eV/A^3)
            popt = [entry["fit_params"]["E0"], entry["fit_params"]["V0"], 
                    entry["fit_params"]["B0_GPa"]/160.21766208, entry["fit_params"]["Bp"]]
            
            fit_E_curve = bm3_eos(fit_V, *popt)
            plt.plot(fit_V, fit_E_curve - E_ref, linestyle='--' if method=="VASP" else '-')
            
    plt.xlabel("Volume ($\AA^3$)")
    plt.ylabel("Energy (eV)")
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
