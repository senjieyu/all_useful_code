import os
import sys
import shutil
import json
import numpy as np
from ase.io import read, write

# ==========================================
# 1. Configuration & Constants
# ==========================================
def load_config():
    # Look for config.json in current dir or parent dir
    if os.path.exists("config.json"):
        with open("config.json", 'r') as f: return json.load(f)
    elif os.path.exists("../config.json"):
        with open("../config.json", 'r') as f: return json.load(f)
    else:
        print("Error: config.json not found.")
        sys.exit(1)

CONFIG = load_config()
ELEMENTS = CONFIG.get("elements", ["Cu"])
MASSES = CONFIG.get("masses", {"Cu": 63.546})
EOS_SCALES = CONFIG.get("eos_scale_factors", [0.94, 1.08])
LMP_BIN = CONFIG["lammps"]["binary"]
# Check for local potential files first
if os.path.exists("mlip.ini"):
    MLIP_INI = "mlip.ini"
else:
    MLIP_INI = CONFIG["lammps"]["potential_files"]["mlip_ini"]

if os.path.exists("current.mtp"):
    MTP_FILE = "current.mtp"
else:
    MTP_FILE = CONFIG["lammps"]["potential_files"]["mtp_file"]

QUEUE = CONFIG["queue"]

# ==========================================
# 2. LAMMPS Helper Functions
# ==========================================
def write_lammps_input(label, atoms, run_type="static", strain_matrix=None):
    """
    Generates LAMMPS input and data file.
    run_type: 'relax', 'static'
    """
    data_fname = f"data_{label}.lmp"
    input_fname = f"in_{label}.lmp"
    log_fname = f"log_{label}.lmp"
    
    # Write data file
    write(data_fname, atoms, format='lammps-data', atom_style='atomic', specorder=ELEMENTS)
    
    with open(input_fname, 'w') as f:
        f.write("clear\n")
        f.write("units metal\n")
        f.write("dimension 3\n")
        f.write("boundary p p p\n")
        f.write("atom_style atomic\n")
        f.write(f"read_data {data_fname}\n")
        
        # Masses
        for i, el in enumerate(ELEMENTS):
            m = MASSES.get(el, 1.0)
            f.write(f"mass {i+1} {m}\n")
            
        # Potentials
        f.write(f"pair_style mlip {MLIP_INI}\n")
        f.write("pair_coeff * *\n")
        
        # Output setup
        f.write("thermo 10\n")
        f.write("thermo_style custom step pe press pxx pyy pzz pyz pxz pxy vol\n")
        
        if run_type == "relax":
            # ISIF 3: Relax ions + cell volume + shape
            f.write("min_style cg\n")
            f.write("minimize 1.0e-10 1.0e-10 1000 10000\n")
            f.write("fix 1 all box/relax iso 0.0 vmax 0.001\n")
            f.write("minimize 1.0e-10 1.0e-10 1000 10000\n")
            f.write("write_data relaxed.lmp\n")
            
        elif run_type == "static":
            f.write("run 0\n")
            
    return input_fname, log_fname

def run_lammps(input_fname, log_fname):
    nproc = QUEUE.get("nproc", 1)
    if nproc > 1:
        cmd = f"mpirun -np {nproc} {LMP_BIN} -in {input_fname} -log {log_fname}"
    else:
        cmd = f"{LMP_BIN} -in {input_fname} -log {log_fname}"
    
    # Run
    ret = os.system(f"{cmd} > /dev/null 2>&1")
    if ret != 0:
        print(f"Error running LAMMPS: {cmd}")
        sys.exit(1)

def parse_log(log_fname):
    """Parses PE, Volume, Stress from log file."""
    # Simple parser for the last thermo line
    # Thermo style: step pe press pxx pyy pzz pyz pxz pxy vol
    import re
    with open(log_fname, 'r') as f:
        lines = f.readlines()
        
    # Find last numeric line
    data = None
    for line in reversed(lines):
        if "Loop time" in line: continue
        try:
            parts = list(map(float, line.split()))
            if len(parts) >= 10: # step + 8 props
                data = parts
                break
        except:
            continue
            
    if data is None: return None
    
    # Mapping based on thermo_style
    # 0:step, 1:pe, 2:press, 3:pxx, 4:pyy, 5:pzz, 6:pyz, 7:pxz, 8:pxy, 9:vol
    pe = data[1]
    press = data[2] # bars
    stress = np.array(data[3:9]) # bars, pxx, pyy...
    vol = data[9]
    
    # Convert Stress (bar) to GPa
    # 1 bar = 0.0001 GPa
    # Note: LAMMPS outputs Pressure, Stress = -Pressure
    # But for pxx etc, these are pressure tensor components.
    # Sigma_xx = - Pxx
    sigma_gpa = -stress * 0.0001
    
    return {"pe": pe, "vol": vol, "stress_gpa": sigma_gpa}

# ==========================================
# 3. Workflows
# ==========================================

def run_relaxation():
    print("--- 1. Relaxation ---")
    atoms = read("start.lmp", format='lammps-data', style='atomic')
    # Set species for mass writing
    if len(ELEMENTS) == 1:
        atoms.symbols = [ELEMENTS[0]] * len(atoms)
        
    inp, log = write_lammps_input("relax", atoms, "relax")
    run_lammps(inp, log)
    
    # Load relaxed structure
    if not os.path.exists("relaxed.lmp"):
        print("Error: Relaxation failed, relaxed.lmp not found.")
        sys.exit(1)
        
    relaxed_atoms = read("relaxed.lmp", format='lammps-data', style='atomic')
    if len(ELEMENTS) == 1:
        relaxed_atoms.symbols = [ELEMENTS[0]] * len(relaxed_atoms)
        
    print("Relaxation Done.")
    return relaxed_atoms

def run_eos(atoms):
    print("--- 2. EOS Calculation ---")
    results = []
    
    cell0 = atoms.get_cell()
    spos = atoms.get_scaled_positions()
    
    for scale in EOS_SCALES:
        # Scale volume: V = V0 * scale
        # L = L0 * scale^(1/3)
        ls = scale ** (1.0/3.0)
        
        scaled_atoms = atoms.copy()
        scaled_atoms.set_cell(cell0 * ls, scale_atoms=False) # Only cell, positions are scaled
        scaled_atoms.set_scaled_positions(spos)
        
        label = f"eos_{scale:.4f}"
        inp, log = write_lammps_input(label, scaled_atoms, "static")
        run_lammps(inp, log)
        
        data = parse_log(log)
        results.append({
            "scale": scale,
            "volume": data["vol"],
            "energy": data["pe"],
            "energy_per_atom": data["pe"] / len(atoms)
        })
        print(f"  V_scale={scale:.2f}  E={data['pe']:.4f} eV")
        
    # Save results
    with open("results_eos.json", 'w') as f:
        json.dump(results, f, indent=4)
    print("EOS Results saved to results_eos.json")

def run_elastic(atoms):
    print("--- 3. Elastic Calculation ---")
    delta = 0.01
    
    # 0:xx, 1:yy, 2:zz, 3:yz, 4:xz, 5:xy
    C = np.zeros((6, 6))
    
    # Reference stress (should be near 0)
    inp0, log0 = write_lammps_input("elas_0", atoms, "minimize")
    run_lammps(inp0, log0)
    # s0 = parse_log(log0)["stress_gpa"] # Use if not zero
    
    cell0 = atoms.get_cell()
    
    for j in range(6):
        # Apply +delta
        D_plus = np.eye(3)
        if j == 0: D_plus[0,0] += delta
        elif j == 1: D_plus[1,1] += delta
        elif j == 2: D_plus[2,2] += delta
        elif j == 3: D_plus[1,2] += delta
        elif j == 4: D_plus[0,2] += delta
        elif j == 5: D_plus[0,1] += delta
        
        atoms_plus = atoms.copy()
        atoms_plus.set_cell(D_plus @ cell0, scale_atoms=True)
        inp_p, log_p = write_lammps_input(f"elas_p{j}", atoms_plus, "static")
        run_lammps(inp_p, log_p)
        s_plus = parse_log(log_p)["stress_gpa"]
        
        # Apply -delta
        D_minus = np.eye(3)
        if j == 0: D_minus[0,0] -= delta
        elif j == 1: D_minus[1,1] -= delta
        elif j == 2: D_minus[2,2] -= delta
        elif j == 3: D_minus[1,2] -= delta
        elif j == 4: D_minus[0,2] -= delta
        elif j == 5: D_minus[0,1] -= delta
        
        atoms_minus = atoms.copy()
        atoms_minus.set_cell(D_minus @ cell0, scale_atoms=True)
        inp_m, log_m = write_lammps_input(f"elas_m{j}", atoms_minus, "static")
        run_lammps(inp_m, log_m)
        s_minus = parse_log(log_m)["stress_gpa"]
        
        # Cij = (sigma+ - sigma-) / (2*delta)
        C[:, j] = (s_plus - s_minus) / (2 * delta)
        
    # Summarize
    # Moduli (Voigt)
    Kv = (C[0,0]+C[1,1]+C[2,2]+2*(C[0,1]+C[1,2]+C[2,0]))/9.0
    Gv = ((C[0,0]+C[1,1]+C[2,2])-(C[0,1]+C[1,2]+C[2,0])+3*(C[3,3]+C[4,4]+C[5,5]))/15.0
    
    results = {
        "Cij_GPa": C.tolist(),
        "Bulk_Modulus_Kv_GPa": Kv,
        "Shear_Modulus_Gv_GPa": Gv
    }
    
    with open("results_elastic.json", 'w') as f:
        json.dump(results, f, indent=4)
    print("Elastic Results saved to results_elastic.json")
    print("Cij Matrix (GPa):")
    print(C)

def main():
    if not os.path.exists("start.lmp"):
        print("Error: start.lmp not found.")
        sys.exit(1)
        
    # 1. Relax
    relaxed_atoms = run_relaxation()
    
    # 2. EOS
    run_eos(relaxed_atoms)
    
    # 3. Elastic
    run_elastic(relaxed_atoms)

if __name__ == "__main__":
    main()
