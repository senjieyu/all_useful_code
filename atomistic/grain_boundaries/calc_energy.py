import os
import sys

def calculate_energy(outcar_path="OUTCAR"):
    if not os.path.exists(outcar_path):
        print(f"Error: {outcar_path} does not exist.")
        sys.exit(1)

    n_ions = 0
    toten = 0.0
    found_nions = False
    found_toten = False

    with open(outcar_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if "NIONS" in line and not found_nions:
                # e.g., "   number of ions     NIONS =      4"
                parts = line.split("=")
                if len(parts) > 1:
                    try:
                        n_ions = int(parts[-1].strip().split()[0])
                        found_nions = True
                    except ValueError:
                        pass
            if "free  energy   TOTEN" in line:
                # e.g., "  free  energy   TOTEN  =       -14.32001234 eV"
                parts = line.split("=")
                if len(parts) > 1:
                    val_str = parts[-1].split("eV")[0].strip()
                    try:
                        toten = float(val_str)
                        found_toten = True
                    except ValueError:
                        pass

    if not found_nions:
        print("Error: Could not find number of atoms (NIONS) in OUTCAR.")
        sys.exit(1)
    if not found_toten:
        print("Error: Could not find total energy (TOTEN) in OUTCAR.")
        sys.exit(1)

    avg_energy = toten / n_ions

    print(f"Total atoms (NIONS)        = {n_ions}")
    print(f"Total energy (TOTEN)       = {toten:.6f} eV")
    print(f"Average energy per atom    = {avg_energy:.6f} eV/atom")

if __name__ == "__main__":
    calculate_energy()
