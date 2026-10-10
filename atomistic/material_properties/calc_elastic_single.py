import numpy as np
import sys
import os
import subprocess
import re

# Constants
EV_A3_TO_GPA = 160.21766208

def read_lammps_stress(log_file):
    """
    Reads the LAST pressure tensor from LAMMPS log file.
    Assumes thermo output: pxx pyy pzz pyz pxz pxy
    Note: LAMMPS 'pressure' is usually -stress. 
    However, we want the stress tensor (sigma). 
    Sigma = -Pressure.
    """
    # Pattern to find thermo lines. 
    # We look for the line after the header. 
    # Or just grab the last numeric line that matches the thermo format.
    # But reliable way: look for "Step Pxx Pyy Pzz Pyz Pxz Pxy" header
    
    last_stress = None
    
    with open(log_file, 'r') as f:
        lines = f.readlines()
        
    # Reverse search for the data
    # Typical log:
    # Step Pxx Pyy Pzz Pyz Pxz Pxy
    # 100 -1.2 -1.2 -1.2 0 0 0
    # Loop time ...
    
    for i in range(len(lines)-1, -1, -1):
        line = lines[i].strip()
        if "Loop time" in line:
            # The data is before this
            continue
        
        # Try to parse 7 floats (Step + 6 stress comps) or just 6 if we customized
        # In our script, we will define thermo_style custom step pxx pyy pzz pyz pxz pxy
        try:
            parts = list(map(float, line.split()))
            if len(parts) >= 7: # step + 6 components
                # pxx, pyy, pzz, pyz, pxz, pxy
                # LAMMPS order: xx, yy, zz, yz, xz, xy
                # Stress = -Pressure
                press_tensor = np.array(parts[1:7])
                stress_tensor = -press_tensor # Convert P to Sigma
                
                # Voigt order for Cij calculation typically: xx, yy, zz, yz, xz, xy
                # (matches LAMMPS order exactly)
                return stress_tensor
        except ValueError:
            continue
            
    return None

def write_lammps_input(fname, data_file, strain_tensor=None):
    """
    Writes a LAMMPS input file.
    If strain_tensor is provided (Voigt: e1..e6 or 3x3), apply deformation.
    """
    
    # Strain application string
    # We use 'change_box' or 'deform'
    # But 'change_box' with 'remap x' is good for static.
    # However, simply modifying the box vectors in the data file or 
    # using 'change_box' commands is tricky if triclinic.
    # Easier: use 'change_box' command explicitly.
    
    deform_cmd = ""
    if strain_tensor is not None:
        # Strain tensor is likely [exx, eyy, ezz, eyz, exz, exy]
        # LAMMPS change_box uses: x final ... y final ... xy final ...
        # For small strains, we can assume:
        # Lx' = Lx * (1 + exx)
        # Ly' = Ly * (1 + eyy)
        # Lz' = Lz * (1 + ezz)
        # Tilt factors change too.
        # This is complicated to script purely in LAMMPS input for generic triclinic.
        pass

    # Instead of complex Python-generation of inputs, 
    # we will use a simpler approach: 
    # We already have a relaxed structure. 
    # We will let LAMMPS handle the deformation if we pass the right command.
    # OR, we modify the box in Python (ASE) and write a new data file.
    # This script assumes 'ASE' is available to modify the cell.
    pass

# We will rely on ASE to handle file I/O and cell deformation 
# because it handles the box vector math correctly.
from ase.io import read, write
from ase.constraints import UnitCellFilter

def get_stress_from_lammps(atoms, label, lmp_bin, pot_file, mlip_ini):
    """
    Runs LAMMPS for a given atoms object and returns stress (Voigt 6-vector, GPa).
    """
    data_fname = f"data_{label}.lmp"
    input_fname = f"in_{label}.lmp"
    log_fname = f"log_{label}.lmp"
    
    # Write data file
    # Important: atom_style atomic for MLIP
    write(data_fname, atoms, format='lammps-data', atom_style='atomic', specorder=['Cu'])
    
    # Write input file
    with open(input_fname, 'w') as f:
        f.write("clear\n")
        f.write("units metal\n")
        f.write("dimension 3\n")
        f.write("boundary p p p\n")
        f.write("atom_style atomic\n")
        f.write(f"read_data {data_fname}\n")
        f.write("mass 1 63.546\n") # Cu mass
        
        f.write(f"pair_style mlip {mlip_ini}\n")
        f.write("pair_coeff * *\n")
        
        f.write("neighbor 2.0 bin\n")
        f.write("neigh_modify delay 0 every 1 check yes\n")
        
        # Static calculation
        f.write("thermo 1\n")
        f.write("thermo_style custom step pxx pyy pzz pyz pxz pxy\n")
        f.write("run 0\n")
    
    # Run LAMMPS
    cmd = f"{lmp_bin} -in {input_fname} -log {log_fname}"
    # Suppress stdout
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL)
    
    # Parse stress
    stress = read_lammps_stress(log_fname)
    if stress is None:
        raise RuntimeError(f"Could not read stress from {log_fname}")
        
    # LAMMPS metal units: pressure in bars.
    # Convert bars to GPa. 1 bar = 1e-4 GPa.
    # Wait, EV_A3_TO_GPA is 160.21... 
    # LAMMPS 'metal' units pressure is [bars].
    # We want GPa.
    # 1 bar = 0.0001 GPa
    
    stress_gpa = stress * 0.0001 
    return stress_gpa

def main():
    if len(sys.argv) < 2:
        print("Usage: python calc_elastic_single.py <structure_file> [lmp_bin]")
        sys.exit(1)
        
    struct_file = sys.argv[1]
    lmp_bin = sys.argv[2] if len(sys.argv) > 2 else "lmp_intel_cpu_intelmpi"
    
    # Ensure config files exist
    if not os.path.exists("mlip.ini") or not os.path.exists("current.mtp"):
        print("Error: mlip.ini or current.mtp missing.")
        sys.exit(1)
        
    # 1. Load and Relax Structure
    # We do the relaxation in LAMMPS to be consistent
    # Or we can just trust the input is close? 
    # Better to relax first to zero stress.
    
    print("Relaxing structure...")
    # Explicitly specify format since extension .lmp is ambiguous for ASE
    # It is a LAMMPS data file.
    atoms0 = read(struct_file, format='lammps-data', style='atomic')
    
    # ASE's lammps-data reader (style=atomic) assigns species as 'H' (atomic number 1)
    # or 'He' etc based on type ID if it doesn't know.
    # We must enforce Cu species for writing, because we use specorder=['Cu'].
    # If atoms.symbols are 'H', but specorder is ['Cu'], it fails.
    atoms0.symbols = ['Cu'] * len(atoms0)
    
    # Relax using a dedicated run
    relax_in = "relax_pre.in"
    relax_data = "relax_pre.data"
    relax_final = "relaxed.data"
    
    write(relax_data, atoms0, format='lammps-data', atom_style='atomic', specorder=['Cu'])
    
    with open(relax_in, 'w') as f:
        f.write("clear\n")
        f.write("units metal\n")
        f.write("atom_style atomic\n")
        f.write(f"read_data {relax_data}\n")
        f.write("mass 1 63.546\n")
        f.write("pair_style mlip mlip.ini\n")
        f.write("pair_coeff * *\n")
        f.write("min_style cg\n")
        f.write("minimize 1.0e-10 1.0e-10 1000 10000\n") # Relax atoms
        f.write("fix 1 all box/relax iso 0.0 vmax 0.001\n") # Relax cell
        f.write("minimize 1.0e-10 1.0e-10 1000 10000\n")
        f.write(f"write_data {relax_final}\n")
        
    cmd = f"{lmp_bin} -in {relax_in} -log relax.log"
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL)
    
    # Load relaxed structure
    atoms = read(relax_final, format='lammps-data', style='atomic')
    # Assign species manually since lammps-data doesn't store it for atomic
    atoms.symbols = ['Cu'] * len(atoms) 
    
    # Get reference volume
    vol = atoms.get_volume()
    
    # 2. Calculate Elastic Constants
    # Method: Apply +/- delta strain to each of 6 components
    # Calculate stress difference
    # C_ij approx dSigma_i / dStrain_j
    
    delta = 0.01 # 1% strain
    
    # Voigt map: 0:xx, 1:yy, 2:zz, 3:yz, 4:xz, 5:xy
    
    stress0 = get_stress_from_lammps(atoms, "0", lmp_bin, "current.mtp", "mlip.ini")
    
    C = np.zeros((6, 6))
    
    print("Calculating elastic constants...")
    
    # We iterate over 6 strain components
    for j in range(6):
        # Apply +delta
        atoms_plus = atoms.copy()
        cell0 = atoms.get_cell()
        
        # Construct strain matrix
        strain_mat = np.eye(3)
        
        if j == 0: strain_mat[0,0] += delta
        elif j == 1: strain_mat[1,1] += delta
        elif j == 2: strain_mat[2,2] += delta
        elif j == 3: 
            strain_mat[1,2] += delta
            strain_mat[2,1] += delta # Symmetric pure shear? 
            # Voigt strain usually defines gamma_yz = 2 * epsilon_yz
            # If we strain by delta in Voigt notation:
            # gamma_yz = delta.
            # So epsilon_yz = delta / 2.
            # Let's stick to engineering strain definitions for Cij
            # Or standard definition: sigma = C * epsilon
            # If we apply strain epsilon_j, C_ij = sigma_i / epsilon_j
            pass
        elif j == 4: 
            strain_mat[0,2] += delta
            strain_mat[2,0] += delta
        elif j == 5: 
            strain_mat[0,1] += delta
            strain_mat[1,0] += delta
            
        # Wait, strictly speaking:
        # Normal strains (0,1,2): box length change.
        # Shear strains (3,4,5): tilt.
        # ASE set_cell with scale_atoms=True handles this if we provide the new cell matrix.
        # New Cell = Strain * Old Cell
        
        # Correct strain matrix for Voigt (j=3,4,5 are shears)
        # For shear, we apply gamma. 
        # epsilon_tensor = 0.5 * gamma
        # Deformation gradient F = I + epsilon_tensor (small strain)
        # Actually F = I + strain_matrix_displacement_gradient
        
        # Let's use a simpler "box deformation" logic:
        # e_xx: x -> x(1+d)
        # e_yz: y -> y + z*d/2 ? No.
        
        # Standard approach:
        # Apply deformation matrix D. new_cell = D @ old_cell
        D = np.eye(3)
        if j == 0: D[0,0] += delta
        elif j == 1: D[1,1] += delta
        elif j == 2: D[2,2] += delta
        elif j == 3: D[1,2] += delta # Shear yz (tilt y along z) -> this is gamma_yz = delta
        elif j == 4: D[0,2] += delta # Shear xz
        elif j == 5: D[0,1] += delta # Shear xy
        
        # Note: applying simple shear like D[0,1] corresponds to gamma_xy = delta
        # and epsilon_xy = delta/2.
        # The stress-strain relation is sigma_i = C_ij * epsilon_j (Voigt notation).
        # In Voigt notation, strain vector e = [exx, eyy, ezz, 2exy, 2exz, 2eyz] ??
        # Or [exx, eyy, ezz, gyz, gxz, gxy].
        # Usually C matrix relates stress (xx,yy,zz,yz,xz,xy) to strain (xx,yy,zz,yz,xz,xy) ?
        # Be careful with factors of 2.
        # Common convention: Sigma_voigt = C_voigt * Epsilon_voigt
        # Where Epsilon_voigt has gammas (2*epsilon) for shears.
        # So if we apply D[0,1] = delta, then gamma_xy = delta.
        # The strain is exactly delta in the Voigt vector index 5.
        
        atoms_plus.set_cell(D @ cell0, scale_atoms=True)
        s_plus = get_stress_from_lammps(atoms_plus, f"p{j}", lmp_bin, "current.mtp", "mlip.ini")
        
        # Apply -delta
        atoms_minus = atoms.copy()
        D = np.eye(3)
        if j == 0: D[0,0] -= delta
        elif j == 1: D[1,1] -= delta
        elif j == 2: D[2,2] -= delta
        elif j == 3: D[1,2] -= delta
        elif j == 4: D[0,2] -= delta
        elif j == 5: D[0,1] -= delta
        
        atoms_minus.set_cell(D @ cell0, scale_atoms=True)
        s_minus = get_stress_from_lammps(atoms_minus, f"m{j}", lmp_bin, "current.mtp", "mlip.ini")
        
        # Finite difference
        # C_ij = (sigma_i(+) - sigma_i(-)) / (2 * delta)
        C[:, j] = (s_plus - s_minus) / (2 * delta)

    # Output Results
    print("Elastic Constants (GPa):")
    print(C)
    
    # Calculate Moduli
    # Bulk Modulus (Voigt)
    Kv = (C[0,0] + C[1,1] + C[2,2] + 2*(C[0,1] + C[1,2] + C[2,0])) / 9.0
    
    # Shear Modulus (Voigt)
    Gv = ((C[0,0] + C[1,1] + C[2,2]) - (C[0,1] + C[1,2] + C[2,0]) + 3*(C[3,3] + C[4,4] + C[5,5])) / 15.0
    
    # Reuss bounds (Compliance S = C^-1)
    try:
        S = np.linalg.inv(C)
        Kr = 1.0 / ((S[0,0] + S[1,1] + S[2,2]) + 2*(S[0,1] + S[1,2] + S[2,0]))
        Gr = 15.0 / (4*(S[0,0] + S[1,1] + S[2,2]) - 4*(S[0,1] + S[1,2] + S[2,0]) + 3*(S[3,3] + S[4,4] + S[5,5]))
    except np.linalg.LinAlgError:
        Kr = 0.0
        Gr = 0.0
        
    Kvrh = (Kv + Kr) / 2.0
    Gvrh = (Gv + Gr) / 2.0
    
    # Poisson's ratio
    # nu = (3K - 2G) / (2(3K + G))
    if (3*Kvrh + Gvrh) != 0:
        nu = (3*Kvrh - 2*Gvrh) / (2*(3*Kvrh + Gvrh))
    else:
        nu = 0.0
        
    # Anisotropy (Universal)
    # Au = 5 * (Gv/Gr) + (Kv/Kr) - 6
    if Gr != 0 and Kr != 0:
        Au = 5 * (Gv/Gr) + (Kv/Kr) - 6
    else:
        Au = 0.0
        
    # Save to file
    out_txt = "elastic_properties.txt"
    with open(out_txt, 'w') as f:
        f.write("Stiffness Tensor Cij (GPa)\n")
        np.savetxt(f, C, fmt="%10.2f")
        f.write("\n")
        f.write(f"Bulk Modulus Kv:   {Kv:.2f} GPa\n")
        f.write(f"Bulk Modulus Kr:   {Kr:.2f} GPa\n")
        f.write(f"Bulk Modulus Kvrh: {Kvrh:.2f} GPa\n")
        f.write("\n")
        f.write(f"Shear Modulus Gv:   {Gv:.2f} GPa\n")
        f.write(f"Shear Modulus Gr:   {Gr:.2f} GPa\n")
        f.write(f"Shear Modulus Gvrh: {Gvrh:.2f} GPa\n")
        f.write("\n")
        f.write(f"Elastic Anisotropy: {Au:.2f}\n")
        f.write(f"Poisson's Ratio:    {nu:.2f}\n")
        
    print(f"Results saved to {out_txt}")

if __name__ == "__main__":
    main()
