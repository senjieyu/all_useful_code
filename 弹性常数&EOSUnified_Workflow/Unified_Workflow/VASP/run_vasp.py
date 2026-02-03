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
        with open("config.json", 'r') as f:
            return json.load(f)
    elif os.path.exists("../config.json"):
        with open("../config.json", 'r') as f:
            return json.load(f)
    return {}

CONFIG = load_config()
EOS_SCALES = CONFIG.get("eos_scale_factors", [0.94, 1.08])
VAS_CONF = CONFIG.get("vasp", {})
VASP_BIN = VAS_CONF.get("binary", "vasp_std")
QUEUE = CONFIG.get("queue", {"nproc": 40})
NPROC = QUEUE["nproc"]

# Optional per-step INCAR overrides from config.json
# Example:
# "incar_patch": {
#   "relax":   {"ISMEAR": 1, "SIGMA": 0.2},
#   "eos":     {"KSPACING": 0.15, "KGAMMA": true, "ISMEAR": 1, "SIGMA": 0.2},
#   "elastic": {"KSPACING": 0.12, "KGAMMA": true, "ISMEAR": 1, "SIGMA": 0.2}
# }
INCAR_PATCH = CONFIG.get("incar_patch", {})

# Robust defaults for Cu / metals
DEFAULT_METAL = {
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "LORBIT": 0,  # Disable LORBIT
    "NEDOS": 301, # Default NEDOS
}

DEFAULT_EOS = {
    **DEFAULT_METAL,
    "KSPACING": 0.15,   # EOS energy differences: 0.12~0.18 typical
    "KGAMMA": True,
}

DEFAULT_ELASTIC = {
    **DEFAULT_METAL,
    "KSPACING": 0.1,    # Denser k-mesh for stress/elastic
    "KGAMMA": True,
    "PREC": "Accurate",
    "ENCUT": 500,       # High cutoff required for accurate stress/elastic
    "EDIFF": 1e-7,
    "LREAL": False,
    "ADDGRID": True,
    "LASPH": True,
    "ISYM": 0,          # Symmetry off often safer for finite differences
}

# ==========================================
# 2. Helpers
# ==========================================
def generate_incar(mode, output_filename="INCAR"):
    """
    Generates INCAR using vaspkit (101).
    mode: 'SR' (Relax), 'ST' (Static/EOS), etc.
    """
    vaspkit_bin = CONFIG.get("vaspkit", {}).get("binary", "vaspkit")
    vaspkit_bin = os.path.expanduser(vaspkit_bin)
    
    # Input for vaspkit: 101 -> mode
    # e.g. "101\nSR\n"
    inp_str = f"101\n{mode}\n"
    
    print(f"Generating INCAR ({mode}) using vaspkit...")
    try:
        # vaspkit generates 'INCAR' in current dir
        if os.path.exists("INCAR"): os.remove("INCAR")
        
        cmd = f"{vaspkit_bin}"
        # Use shell=True for piping
        # Note: vaspkit might need interactive input
        process = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate(input=inp_str.encode())
        
        if process.returncode != 0:
            print(f"Error running vaspkit: {err.decode()}")
            return False
            
        if os.path.exists("INCAR"):
            # Add NCORE for parallel efficiency if not present
            with open("INCAR", "a") as f:
                f.write("\nNCORE = 4\n")
                
            if output_filename != "INCAR":
                shutil.move("INCAR", output_filename)
            return True
        else:
            print("Error: vaspkit did not generate INCAR.")
            print(out.decode())
            return False
    except Exception as e:
        print(f"Exception running vaspkit: {e}")
        return False


def patch_incar(incar_path, overrides):
    """
    Patch an INCAR file by forcing specific tags (overrides).
    - If a key exists: replace its last occurrence.
    - If missing: append at end.
    Keeps vaspkit defaults for everything else.
    """
    if not os.path.exists(incar_path):
        raise FileNotFoundError(f"INCAR not found: {incar_path}")

    with open(incar_path, "r") as f:
        lines = f.readlines()

    # Record last occurrence of each key
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


def run_vasp(folder="."):
    """Runs VASP in the specified folder."""
    cwd = os.getcwd()
    try:
        os.chdir(folder)

        # Check files
        if (not os.path.exists("INCAR")
                or not os.path.exists("POSCAR")
                or not os.path.exists("POTCAR")):
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
    if not os.path.exists(outcar):
        return None, None

    energy = None
    volume = None

    with open(outcar, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if "volume of cell :" in line:
            parts = line.split()
            try:
                volume = float(parts[-1])
            except Exception:
                pass
        if "free  energy   TOTEN" in line:
            parts = line.split()
            try:
                energy = float(parts[-2])
            except Exception:
                pass

    return energy, volume


def get_elastic_tensor(folder="."):
    """Parses Elastic Tensor from OUTCAR (TOTAL ELASTIC MODULI)."""
    outcar = os.path.join(folder, "OUTCAR")
    if not os.path.exists(outcar):
        return None

    C = np.zeros((6, 6))
    found = False

    with open(outcar, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI (kBar)" in line:
            found = True
            start = i + 3
            for r in range(6):
                parts = lines[start + r].split()
                vals = [float(x) for x in parts[1:7]]
                C[r, :] = vals
            break
            
    if found:
        return C * 0.1  # kBar -> GPa
        
    # Retry with looser parsing if exact match failed
    # Sometimes VASP 5.4.4 OUTCAR format differs slightly
    # Look for the section and parse the 6x6 matrix
    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI" in line and "kBar" in line:
            # Found header, but maybe offset is different?
            # Look ahead for "XX" line
            for j in range(1, 10):
                if "XX" in lines[i+j]:
                    start = i + j + 2 # Skip header line and dashed line
                    try:
                        for r in range(6):
                            parts = lines[start + r].split()
                            # parts: Direction, XX, YY, ZZ, XY, YZ, ZX
                            # We want indices 1 to 6
                            vals = [float(x) for x in parts[1:7]]
                            C[r, :] = vals
                        return C * 0.1
                    except:
                        pass
            break
            
    return None


# ==========================================
# 3. Workflow
# ==========================================
def main():
    # Check POTCAR
    if not os.path.exists("POTCAR"):
        print("POTCAR not found. Trying to generate with vaspkit (103)...")

        vaspkit_bin = CONFIG.get("vaspkit", {}).get("binary", "vaspkit")
        vaspkit_bin = os.path.expanduser(vaspkit_bin)

        try:
            subprocess.check_call(f"{vaspkit_bin} -h > /dev/null 2>&1", shell=True)
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
    
    relax_dir = "relaxation"
    os.makedirs(relax_dir, exist_ok=True)

    if not generate_incar("SR", os.path.join(relax_dir, "INCAR")):
        print("Failed to generate INCAR for relaxation.")
        sys.exit(1)
        
    # Patch Relax INCAR (keep smearing consistent)
    relax_overrides = DEFAULT_METAL.copy()
    relax_overrides.update(INCAR_PATCH.get("relax", {}))
    # Ensure ISIF=3 for full relaxation
    relax_overrides["ISIF"] = 3
    relax_overrides["NCORE"] = 4
    patch_incar(os.path.join(relax_dir, "INCAR"), relax_overrides)
    
    shutil.copy("POSCAR", os.path.join(relax_dir, "POSCAR"))
    shutil.copy("POTCAR", os.path.join(relax_dir, "POTCAR"))

    print("Running Relaxation...")
    if run_vasp(relax_dir):
        contcar = os.path.join(relax_dir, "CONTCAR")
        if os.path.exists(contcar) and os.path.getsize(contcar) > 0:
            shutil.copy(contcar, "POSCAR_relaxed")
            print("Relaxation Done.")
        else:
            print("Error: CONTCAR empty or missing.")
            sys.exit(1)
    else:
        sys.exit(1)

    atoms = read("POSCAR_relaxed")

    # --- Step 2: EOS ---
    print("--- 2. EOS Calculation ---")

    if not generate_incar("ST", "INCAR_eos_template"):
        print("Failed to generate INCAR for EOS.")
        sys.exit(1)

    eos_overrides = DEFAULT_EOS.copy()
    eos_overrides.update(INCAR_PATCH.get("eos", {}))
    patch_incar("INCAR_eos_template", eos_overrides)

    eos_results = []
    cell0 = atoms.get_cell()
    spos = atoms.get_scaled_positions()

    for scale in EOS_SCALES:
        sub_dir = f"eos_{scale:.4f}"
        os.makedirs(sub_dir, exist_ok=True)

        shutil.copy("INCAR_eos_template", os.path.join(sub_dir, "INCAR"))
        shutil.copy("POTCAR", os.path.join(sub_dir, "POTCAR"))
        # No KPOINTS (using KSPACING in INCAR)

        ls = scale ** (1.0 / 3.0)
        scaled = atoms.copy()
        scaled.set_cell(cell0 * ls, scale_atoms=False)
        scaled.set_scaled_positions(spos)
        write(os.path.join(sub_dir, "POSCAR"), scaled, format='vasp', direct=True)

        print(f"Running EOS scale {scale}...")
        run_vasp(sub_dir)

        en, vol = get_energy_volume(sub_dir)
        if en is not None:
            eos_results.append({
                "scale": scale,
                "volume": vol,
                "energy": en,
                "energy_per_atom": en / len(atoms)
            })
            print(f"  V={vol:.2f} E={en:.6f}")

    with open("results_eos.json", 'w') as f:
        json.dump(eos_results, f, indent=4)

    # --- Step 3: Elastic ---
    print("--- 3. Elastic Calculation ---")
    elas_dir = "elastic_calc"
    os.makedirs(elas_dir, exist_ok=True)

    if not generate_incar("ST", "INCAR_elastic_template"):
        print("Failed to generate INCAR for Elastic.")
        sys.exit(1)

    elas_overrides = DEFAULT_ELASTIC.copy()
    elas_overrides.update(INCAR_PATCH.get("elastic", {}))
    patch_incar("INCAR_elastic_template", elas_overrides)

    # Append Elastic driver tags (robust)
    with open("INCAR_elastic_template", "a") as f:
        f.write("\n# Elastic Calculation Driver\n")
        f.write("IBRION  = 6\n")
        f.write("ISIF    = 3\n")
        f.write("NFREE   = 2\n")
        f.write("POTIM   = 0.01\n")   # Smaller displacement for linear regime
        f.write("NSW     = 1\n")
        # f.write("NCORE   = 4\n")    # Commented out NCORE for IBRION=6 stability

    shutil.copy("INCAR_elastic_template", os.path.join(elas_dir, "INCAR"))
    shutil.copy("POTCAR", os.path.join(elas_dir, "POTCAR"))
    shutil.copy("POSCAR_relaxed", os.path.join(elas_dir, "POSCAR"))

    print("Running Elastic Calculation...")
    run_vasp(elas_dir)

    C = get_elastic_tensor(elas_dir)
    if C is not None:
        # Reorder C to standard Voigt (11,22,33,23,13,12)
        C_voigt = np.zeros((6, 6))

        # VASP order: XX, YY, ZZ, XY, YZ, ZX
        # Voigt:     11, 22, 33, 23, 13, 12
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
        print("Elastic Results saved.")
    else:
        print("Failed to parse Elastic Constants.")


if __name__ == "__main__":
    main()
