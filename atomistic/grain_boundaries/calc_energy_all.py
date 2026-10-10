import os
import sys

def parse_vasp(outcar_path):
    if not os.path.exists(outcar_path):
        return None, None
    n_ions = 0
    toten = 0.0
    found_nions = False
    found_toten = False
    with open(outcar_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if "NIONS" in line and not found_nions:
                parts = line.split("=")
                if len(parts) > 1:
                    try:
                        n_ions = int(parts[-1].strip().split()[0])
                        found_nions = True
                    except ValueError:
                        pass
            if "free  energy   TOTEN" in line:
                parts = line.split("=")
                if len(parts) > 1:
                    val_str = parts[-1].split("eV")[0].strip()
                    try:
                        toten = float(val_str)
                        found_toten = True
                    except ValueError:
                        pass
    if found_nions and found_toten:
        return n_ions, toten
    return None, None

def parse_lammps(log_path):
    if not os.path.exists(log_path):
        return None, None
    n_ions = 0
    toten = 0.0
    found_nions = False
    found_toten = False
    with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # We want to match the actual output line, not the print command echo
            # e.g., "Final_Energy = -14.1406046129829"
            if "Num_Atoms =" in line:
                try:
                    val = line.split("=")[1].strip()
                    n_ions = int(val)
                    found_nions = True
                except ValueError:
                    pass
            elif "Final_Energy =" in line:
                try:
                    val = line.split("=")[1].strip()
                    toten = float(val)
                    found_toten = True
                except ValueError:
                    pass
    if found_nions and found_toten:
        return n_ions, toten
    return None, None

def print_result(method, n_ions, toten):
    print(f"=== {method} ===")
    if n_ions is None or toten is None:
        print("Error: Could not extract data or data incomplete.\n")
        return
    avg_energy = toten / n_ions
    print(f"Total atoms (NIONS)        = {n_ions}")
    print(f"Total energy (TOTEN)       = {toten:.6f} eV")
    print(f"Average energy per atom    = {avg_energy:.6f} eV/atom\n")

if __name__ == "__main__":
    vasp_n, vasp_e = parse_vasp("VASP/OUTCAR")
    eam_n, eam_e = parse_lammps("LAMMPS_EAM/log.lammps")
    sus_n, sus_e = parse_lammps("LAMMPS_SUS2/log.lammps")

    print_result("VASP", vasp_n, vasp_e)
    print_result("LAMMPS_EAM", eam_n, eam_e)
    print_result("LAMMPS_SUS2", sus_n, sus_e)
