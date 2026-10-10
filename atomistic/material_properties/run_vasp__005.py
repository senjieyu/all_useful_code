import os
import sys
import shutil
import json
import subprocess
import re
import numpy as np
from ase.io import read, write

# ==========================================
# 1. Configuration
# ==========================================
def load_config():
    if os.path.exists("config.json"):
        with open("config.json", 'r') as f: return json.load(f)
    elif os.path.exists("../config.json"):
        with open("../config.json", 'r') as f: return json.load(f)
    return {}

CONFIG = load_config()
EOS_SCALES = CONFIG.get("eos_scale_factors", [0.94, 1.08])
VAS_CONF = CONFIG.get("vasp", {})
VASP_BIN = VAS_CONF.get("binary", "vasp_std")
QUEUE = CONFIG.get("queue", {"nproc": 40})
NPROC = QUEUE["nproc"]

# ==========================================
# 2. Helpers
# ==========================================
def run_vasp(folder="."):
    """Runs VASP in the specified folder."""
    cwd = os.getcwd()
    try:
        os.chdir(folder)
        # Check files
        if not os.path.exists("INCAR") or not os.path.exists("POSCAR") or not os.path.exists("POTCAR"):
            print(f"Error: Missing VASP input files in {folder}")
            return False
            
        print(f"Running VASP in {folder} ...")
        cmd = f"mpirun -np {NPROC} {VASP_BIN}"
        with open("vasp_run.log", "w") as f:
            subprocess.check_call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
        return True
    except subprocess.CalledProcessError as e:
        print(f"VASP failed in {folder}: {e}")
        return False
    finally:
        os.chdir(cwd)

def get_energy_volume(folder="."):
    """Parses Energy and Volume from OUTCAR."""
    outcar = os.path.join(folder, "OUTCAR")
    if not os.path.exists(outcar): return None, None
    
    energy = None
    volume = None
    
    with open(outcar, 'r') as f:
        lines = f.readlines()
        
    # Get last energy and volume
    for line in lines:
        if "volume of cell :" in line:
            parts = line.split()
            # line: "volume of cell : 123.45"
            try: volume = float(parts[-1])
            except: pass
        if "free  energy   TOTEN" in line:
            parts = line.split()
            # line: "free  energy   TOTEN  = -123.456 eV"
            try: energy = float(parts[-2])
            except: pass
            
    return energy, volume

def get_elastic_tensor(folder="."):
    """Parses Elastic Tensor from OUTCAR (TOTAL ELASTIC MODULI)."""
    outcar = os.path.join(folder, "OUTCAR")
    if not os.path.exists(outcar): return None
    
    C = np.zeros((6, 6))
    found = False
    
    with open(outcar, 'r') as f:
        lines = f.readlines()
        
    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI (kBar)" in line:
            found = True
            # Next line is header
            # Then 6 lines of matrix
            start = i + 3
            for r in range(6):
                parts = lines[start+r].split()
                # format: XX YY ZZ XY YZ XZ
                # VASP order: XX YY ZZ XY YZ XZ
                # We need to map to 0..5
                # VASP columns: 1:XX, 2:YY, 3:ZZ, 4:XY, 5:YZ, 6:XZ
                # Values are in kBar. 1 kBar = 0.1 GPa.
                vals = [float(x) for x in parts[1:7]]
                C[r, :] = vals
            break
            
    if found:
        return C * 0.1 # Convert kBar to GPa
    return None

# ==========================================
# 3. Workflow
# ==========================================
def main():
    # Check POTCAR
    if not os.path.exists("POTCAR"):
        print("POTCAR not found. Trying to generate with vaspkit (103)...")
        
        # Get vaspkit binary from config
        vaspkit_bin = CONFIG.get("vaspkit", {}).get("binary", "vaspkit")
        vaspkit_bin = os.path.expanduser(vaspkit_bin)
        
        # Using subprocess to capture potential errors
        try:
            # Check if vaspkit is available
            subprocess.check_call(f"{vaspkit_bin} -h > /dev/null 2>&1", shell=True)
            
            # Run vaspkit 103
            cmd = f"(echo 103) | {vaspkit_bin}"
            with open("vaspkit_potcar.log", "w") as log:
                subprocess.check_call(cmd, shell=True, stdout=log, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            print("Error: vaspkit execution failed. Check vaspkit_potcar.log or ensure vaspkit is in PATH.")
            sys.exit(1)
            
        if not os.path.exists("POTCAR"):
            print("Error: POTCAR generation failed. Please provide POTCAR.")
            sys.exit(1)

    # --- Step 1: Relaxation ---
    print("--- 1. Relaxation ---")
    if not os.path.exists("INCAR_relax"):
        print("Error: INCAR_relax missing.")
        sys.exit(1)
        
    shutil.copy("INCAR_relax", "INCAR")
    if run_vasp("."):
        if os.path.exists("CONTCAR") and os.path.getsize("CONTCAR") > 0:
            shutil.copy("CONTCAR", "POSCAR_relaxed")
            print("Relaxation Done.")
        else:
            print("Error: CONTCAR empty or missing.")
            sys.exit(1)
    else:
        sys.exit(1)
        
    # Load relaxed structure
    atoms = read("POSCAR_relaxed")
    
    # --- Step 2: EOS ---
    print("--- 2. EOS Calculation ---")
    eos_results = []
    cell0 = atoms.get_cell()
    spos = atoms.get_scaled_positions()
    
    for scale in EOS_SCALES:
        sub_dir = f"eos_{scale:.4f}"
        os.makedirs(sub_dir, exist_ok=True)
        
        # Prepare Inputs
        shutil.copy("INCAR_eos", os.path.join(sub_dir, "INCAR"))
        shutil.copy("POTCAR", os.path.join(sub_dir, "POTCAR"))
        # No KPOINTS (using KSPACING)
        
        # Scale Structure
        ls = scale ** (1.0/3.0)
        scaled = atoms.copy()
        scaled.set_cell(cell0 * ls, scale_atoms=False)
        scaled.set_scaled_positions(spos)
        write(os.path.join(sub_dir, "POSCAR"), scaled, format='vasp', direct=True)
        
        # Run
        print(f"Running EOS scale {scale}...")
        run_vasp(sub_dir)
        
        # Collect
        en, vol = get_energy_volume(sub_dir)
        if en is not None:
            eos_results.append({
                "scale": scale,
                "volume": vol,
                "energy": en,
                "energy_per_atom": en / len(atoms)
            })
            print(f"  V={vol:.2f} E={en:.4f}")
            
    with open("results_eos.json", 'w') as f:
        json.dump(eos_results, f, indent=4)
        
    # --- Step 3: Elastic ---
    print("--- 3. Elastic Calculation ---")
    elas_dir = "elastic_calc"
    os.makedirs(elas_dir, exist_ok=True)
    
    shutil.copy("INCAR_elastic", os.path.join(elas_dir, "INCAR"))
    shutil.copy("POTCAR", os.path.join(elas_dir, "POTCAR"))
    shutil.copy("POSCAR_relaxed", os.path.join(elas_dir, "POSCAR"))
    
    print("Running Elastic Calculation...")
    run_vasp(elas_dir)
    
    C = get_elastic_tensor(elas_dir)
    if C is not None:
        # Calculate Moduli
        # VASP C matrix: 1:XX, 2:YY, 3:ZZ, 4:XY, 5:YZ, 6:XZ
        # Voigt standard: XX, YY, ZZ, YZ, XZ, XY (indices 0,1,2,3,4,5)
        # VASP indices 0,1,2 map to Voigt 0,1,2
        # VASP index 3 (XY) maps to Voigt 5
        # VASP index 4 (YZ) maps to Voigt 3
        # VASP index 5 (XZ) maps to Voigt 4
        
        # Reorder C to standard Voigt
        C_voigt = np.zeros((6,6))
        # Mapping: VASP 0,1,2,3,4,5 -> Voigt 0,1,2,5,3,4
        v_map = [0, 1, 2, 5, 3, 4]
        
        for r in range(6):
            for c in range(6):
                C_voigt[v_map[r], v_map[c]] = C[r, c]
                
        # Calculate Kv, Gv
        # ... (Same formula as before)
        Kv = (C_voigt[0,0]+C_voigt[1,1]+C_voigt[2,2]+2*(C_voigt[0,1]+C_voigt[1,2]+C_voigt[2,0]))/9.0
        Gv = ((C_voigt[0,0]+C_voigt[1,1]+C_voigt[2,2])-(C_voigt[0,1]+C_voigt[1,2]+C_voigt[2,0])+3*(C_voigt[3,3]+C_voigt[4,4]+C_voigt[5,5]))/15.0
        
        res = {
            "Cij_GPa_VASP_Raw": C.tolist(),
            "Cij_GPa_Voigt": C_voigt.tolist(),
            "Kv_GPa": Kv,
            "Gv_GPa": Gv
        }
        with open("results_elastic.json", 'w') as f:
            json.dump(res, f, indent=4)
        print("Elastic Results saved.")
    else:
        print("Failed to parse Elastic Constants.")

if __name__ == "__main__":
    main()
