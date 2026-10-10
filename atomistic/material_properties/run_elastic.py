import os
import sys
import shutil
import json
import subprocess
import re
import numpy as np
from ase.io import read, write

# ==========================================
# 1. Configuration & Helpers (Copied from run_vasp.py)
# ==========================================
def load_config():
    if os.path.exists("config.json"):
        with open("config.json", 'r') as f:
            return json.load(f)
    elif os.path.exists("../config.json"):
        with open("../config.json", 'r') as f:
            return json.load(f)
    return {}

CONFIG = load_config()
VAS_CONF = CONFIG.get("vasp", {})
VASP_BIN = VAS_CONF.get("binary", "vasp_std")
QUEUE = CONFIG.get("queue", {"nproc": 40})
NPROC = QUEUE["nproc"]
INCAR_PATCH = CONFIG.get("incar_patch", {})

DEFAULT_ELASTIC = {
    "ISPIN": 1,
    "LREAL": False,
    "ALGO": "Normal",
    "LWAVE": False,
    "LCHARG": False,
    "ADDGRID": True,
    "LASPH": True,
    "PREC": "Accurate",
    "ISMEAR": 1,
    "SIGMA": 0.1,
    "NELM": 120,
    "EDIFF": 1e-07,
    "NCORE": 4,
    "ENCUT": 500,
    "IBRION": 6,
    "ISIF": 3,
    "NFREE": 2,
}

def get_structure_type():
    # Look for style.txt in parent directories
    paths = ["style.txt", "../style.txt", "../../style.txt"]
    style_file = None
    for p in paths:
        if os.path.exists(p):
            style_file = p
            break
    
    if not style_file:
        print("Warning: style.txt not found. Using default parameters.")
        return None

    cwd = os.getcwd()
    dirname = os.path.basename(cwd)
    idx = None
    if dirname.startswith("struct_"):
        try:
            idx = int(dirname.split("_")[1])
        except:
            pass
            
    if idx is None:
        return None

    try:
        with open(style_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        if int(parts[0]) == idx:
                            stype = parts[1].upper()
                            return stype
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error parsing style.txt: {e}")
    return None

def generate_incar(mode, output_filename="INCAR"):
    vaspkit_bin = CONFIG.get("vaspkit", {}).get("binary", "vaspkit")
    vaspkit_bin = os.path.expanduser(vaspkit_bin)
    inp_str = f"101\n{mode}\n"
    
    print(f"Generating INCAR ({mode}) using vaspkit...")
    try:
        if os.path.exists("INCAR"): os.remove("INCAR")
        cmd = f"{vaspkit_bin}"
        process = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate(input=inp_str.encode())
        
        if os.path.exists("INCAR"):
            with open("INCAR", "a") as f:
                f.write("\nNCORE = 4\n")
            if output_filename != "INCAR":
                shutil.move("INCAR", output_filename)
            return True
        else:
            print("Error: vaspkit did not generate INCAR.")
            return False
    except Exception as e:
        print(f"Exception running vaspkit: {e}")
        return False

def patch_incar(incar_path, overrides):
    if not os.path.exists(incar_path):
        raise FileNotFoundError(f"INCAR not found: {incar_path}")

    with open(incar_path, "r") as f:
        lines = f.readlines()

    key_last_idx = {}
    key_pattern = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=")
    for i, ln in enumerate(lines):
        m = key_pattern.match(ln)
        if m:
            key_last_idx[m.group(1).upper()] = i

    def fmt_value(v):
        if isinstance(v, bool):
            return ".TRUE." if v else ".FALSE."
        if isinstance(v, (int, float)):
            return f"{v}"
        return str(v)

    for k, v in overrides.items():
        ku = k.upper()
        new_line = f"{ku:<8} = {fmt_value(v)}\n"
        if ku in key_last_idx:
            lines[key_last_idx[ku]] = new_line
        else:
            lines.append(new_line)

    with open(incar_path, "w") as f:
        f.writelines(lines)

def remove_incar_tags(incar_path, tags):
    if not os.path.exists(incar_path): return
    with open(incar_path, "r") as f: lines = f.readlines()
    new_lines = []
    for line in lines:
        is_removed = False
        for tag in tags:
            if re.match(r"^\s*" + re.escape(tag) + r"\s*=", line, re.IGNORECASE):
                is_removed = True
                break
        if not is_removed: new_lines.append(line)
    with open(incar_path, "w") as f: f.writelines(new_lines)

import textwrap
def write_kpoints(folder, content):
    dedented_content = textwrap.dedent(content).strip()
    with open(os.path.join(folder, "KPOINTS"), "w") as f:
        f.write(dedented_content + "\n")

def run_vasp(folder="."):
    cwd = os.getcwd()
    try:
        os.chdir(folder)
        if (not os.path.exists("INCAR") or not os.path.exists("POSCAR") or not os.path.exists("POTCAR")):
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

def get_elastic_tensor(folder="."):
    outcar = os.path.join(folder, "OUTCAR")
    if not os.path.exists(outcar): return None
    C = np.zeros((6, 6))
    found = False
    with open(outcar, 'r') as f: lines = f.readlines()
    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI (kBar)" in line:
            found = True
            start = i + 3
            for r in range(6):
                parts = lines[start + r].split()
                vals = [float(x) for x in parts[1:7]]
                C[r, :] = vals
            break
    if found: return C * 0.1
    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI" in line and "kBar" in line:
            for j in range(1, 10):
                if "XX" in lines[i+j]:
                    start = i + j + 2
                    try:
                        for r in range(6):
                            parts = lines[start + r].split()
                            vals = [float(x) for x in parts[1:7]]
                            C[r, :] = vals
                        return C * 0.1
                    except: pass
            break
    return None

# ==========================================
# 3. Main Workflow (Elastic Only)
# ==========================================
def main():
    print("--- Elastic Calculation Only ---")
    
    # Check if POSCAR_relaxed exists
    if not os.path.exists("POSCAR_relaxed"):
        if os.path.exists("relaxation/CONTCAR") and os.path.getsize("relaxation/CONTCAR") > 0:
             print("Found relaxation/CONTCAR, copying to POSCAR_relaxed...")
             shutil.copy("relaxation/CONTCAR", "POSCAR_relaxed")
        else:
             print("Error: POSCAR_relaxed not found. Cannot run elastic calculation without relaxed structure.")
             sys.exit(1)
             
    struct_type = get_structure_type()
    if struct_type:
        print(f"Structure Type Detected: {struct_type}")

    elas_dir = "elastic_calc"
    # Back up previous run if exists? Optional.
    if os.path.exists(elas_dir):
        print(f"Warning: {elas_dir} exists. Overwriting...")
    os.makedirs(elas_dir, exist_ok=True)

    # Generate/Prepare INCAR
    if not generate_incar("ST", "INCAR_elastic_template"):
        print("Failed to generate INCAR for Elastic.")
        sys.exit(1)

    elas_overrides = DEFAULT_ELASTIC.copy()
    elas_overrides.update(INCAR_PATCH.get("elastic", {}))
    
    # Remove NCORE/NPAR/KSPACING/KGAMMA
    if "NCORE" in elas_overrides: del elas_overrides["NCORE"]
    if "NPAR" in elas_overrides: del elas_overrides["NPAR"]
    if "KSPACING" in elas_overrides: del elas_overrides["KSPACING"]
    if "KGAMMA" in elas_overrides: del elas_overrides["KGAMMA"]
    
    remove_incar_tags("INCAR_elastic_template", ["NCORE", "NPAR", "KSPACING", "KGAMMA"])
    
    # KPOINTS Logic
    if struct_type == "HCP":
        kpoints_content = """HCP Cu elastic constants (recommended) 
0 
Gamma 
24 24 15 
0 0 0"""
    else:
        # Default for FCC/BCC (struct_0, struct_1)
        kpoints_content = """Elastic constants KPOINTS 
0 
Gamma 
24 24 24 
0 0 0"""

    write_kpoints(elas_dir, kpoints_content)
    patch_incar("INCAR_elastic_template", elas_overrides)

    # Setup Directory
    shutil.copy("INCAR_elastic_template", os.path.join(elas_dir, "INCAR"))
    if os.path.exists("POTCAR"):
        shutil.copy("POTCAR", os.path.join(elas_dir, "POTCAR"))
    else:
        print("Error: POTCAR not found in current directory.")
        sys.exit(1)
        
    shutil.copy("POSCAR_relaxed", os.path.join(elas_dir, "POSCAR"))

    # Run VASP
    print("Running Elastic Calculation...")
    run_vasp(elas_dir)

    # Parse Results
    C = get_elastic_tensor(elas_dir)
    if C is not None:
        C_voigt = np.zeros((6, 6))
        vasp_to_voigt_map = {0: 0, 1: 1, 2: 2, 3: 5, 4: 3, 5: 4}

        for r in range(6):
            for c in range(6):
                v_r = vasp_to_voigt_map[r]
                v_c = vasp_to_voigt_map[c]
                C_voigt[v_r, v_c] = C[r, c]

        Kv = (C_voigt[0, 0] + C_voigt[1, 1] + C_voigt[2, 2]
              + 2 * (C_voigt[0, 1] + C_voigt[1, 2] + C_voigt[2, 0])) / 9.0

        Gv = ((C_voigt[0, 0] + C_voigt[1, 1] + C_voigt[2, 2])
              - (C_voigt[0, 1] + C_voigt[1, 2] + C_voigt[2, 0])
              + 3 * (C_voigt[3, 3] + C_voigt[4, 4] + C_voigt[5, 5])) / 15.0

        res = {
            "Cij_GPa_VASP_Raw": C.tolist(),
            "Cij_GPa_Voigt": C_voigt.tolist(),
            "Kv_GPa": Kv,
            "Gv_GPa": Gv
        }
        with open("results_elastic.json", 'w') as f:
            json.dump(res, f, indent=4)
        print("Elastic Results saved to results_elastic.json.")
    else:
        print("Failed to parse Elastic Constants.")

if __name__ == "__main__":
    main()