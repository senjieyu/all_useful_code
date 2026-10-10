#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate orthorhombic FCC Cu structures with [100], [110], [111] z-orientation.
Apply tensile strain along z up to 50% (step 2%).
Generate VASP and LAMMPS input files for relaxation (fix z, relax xy).
"""

import os
import shutil
import numpy as np
import argparse
from ase.build import bulk, make_supercell
from ase.io import read, write

# --- Parameters ---
A0_VASP = 3.631
A0_EAM = 3.615
A0_SUS2 = 3.633
EPS_MAX = 0.50
STEP = 0.02
STRAINS = np.arange(0.0, EPS_MAX + STEP/2, STEP) # Include endpoint

# VASP settings
VASP_INCAR_TEMPLATE = """SYSTEM = Cu {orientation} strain {strain:.2f}
ISTART = 0
ICHARG = 2
ENCUT = 500
ISMEAR = 1
SIGMA = 0.1
ISIF = 3
IBRION = 2
NSW = 100
NELM = 120
NELMIN = 6
EDIFF = 1E-6
EDIFFG = -0.002
KSPACING = 0.20
KGAMMA = .FALSE.
LREAL = Auto
NCORE = 4
LWAVE = .FALSE.
LCHARG = .FALSE.
"""

# Template for strain 00
VASP_SUB_SH_RELAX_00 = """#!/bin/bash
#BSUB -J Cu_{orientation}_00
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH

echo "Starting VASP Relaxation (ISIF=3, vasp_std_str_c)..."
# Link the special binary
if [ -f ../../../vasp_std_str_c ]; then
    ln -sf ../../../vasp_std_str_c .
fi

COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP ./vasp_std_str_c"
$COMMAND_std > vasp_relax_00.log

if [ -f CONTCAR ]; then
    echo "Relaxation done."
fi
"""

# Template for strain > 0
VASP_SUB_SH_STRAIN = """#!/bin/bash
#BSUB -J Cu_{orientation}_{strain_pct:02d}
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH

# Link the special binary
if [ -f ../../../vasp_std_str_c ]; then
    ln -sf ../../../vasp_std_str_c .
fi

COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP ./vasp_std_str_c"

echo "Starting VASP relaxation (ISIF=3, vasp_std_str_c)..."
$COMMAND_std > vasp_relax.log
"""

RELAX_XY_SCRIPT = """import os
import shutil
import numpy as np
from ase.io import read, write

# VASP run command (will be executed by os.system)
VASP_CMD = "./run_vasp.sh" 

# Tolerance for stress (kB). 10 kB = 1 GPa.
# We want to be closer to 0 pressure in lateral directions.
TOL = 2.0 

# Bulk modulus guess (kB) for updating cell
# Cu bulk modulus ~ 140 GPa = 1400 kB.
B_MOD = 1400.0 
# Mixing factor (small to be safe)
ALPHA = 0.5 

MAX_STEPS = 15

def get_stress_from_outcar(filename="OUTCAR"):
    s_xx = None
    s_yy = None
    if not os.path.exists(filename):
        return None, None
        
    with open(filename, "r") as f:
        lines = f.readlines()
        
    # Find the LAST occurrence of "in kB"
    for line in reversed(lines):
        if "in kB" in line:
            parts = line.split()
            # Format: in kB  XX YY ZZ XY YZ ZX
            if len(parts) >= 8 and parts[0] == "in" and parts[1] == "kB":
                try:
                    s_xx = float(parts[2])
                    s_yy = float(parts[3])
                    return s_xx, s_yy
                except ValueError:
                    pass
    return None, None

def main():
    print("Starting external relaxation loop (ISIF=2, adjusting XY)...")
    
    for step in range(MAX_STEPS):
        print(f"--- Step {step+1}/{MAX_STEPS} ---")
        
        # 1. Run VASP
        exit_code = os.system(VASP_CMD)
        if exit_code != 0:
            print("Error: VASP execution failed.")
            break
            
        # 2. Check Stress
        s_xx, s_yy = get_stress_from_outcar("OUTCAR")
        if s_xx is None:
            print("Error: Could not read stress from OUTCAR.")
            break
            
        print(f"Current Stress (kB): XX={s_xx:.2f}, YY={s_yy:.2f}")
        
        # 3. Check Convergence
        if abs(s_xx) < TOL and abs(s_yy) < TOL:
            print(f"Converged! |XX| and |YY| < {TOL} kB.")
            # Backup final
            shutil.copy("CONTCAR", "CONTCAR_final")
            shutil.copy("OUTCAR", "OUTCAR_final")
            break
            
        # 4. Update Cell
        # Hooke's Law: delta_L / L = - Stress / Modulus
        eps_x = - s_xx / B_MOD * ALPHA
        eps_y = - s_yy / B_MOD * ALPHA
        
        print(f"Updating cell: eps_x={eps_x:.6f}, eps_y={eps_y:.6f}")
        
        # Read latest structure
        # Use CONTCAR for next step
        if not os.path.exists("CONTCAR") or os.path.getsize("CONTCAR") == 0:
             print("Error: CONTCAR not found or empty.")
             break
             
        atoms = read("CONTCAR")
        cell = atoms.get_cell()
        
        # Scale X and Y vectors
        # Note: Assuming cell is aligned to axes (which we ensured in build.py)
        # cell[0] is X vector, cell[1] is Y vector
        cell[0] *= (1.0 + eps_x)
        cell[1] *= (1.0 + eps_y)
        
        atoms.set_cell(cell, scale_atoms=True)
        
        # Write new POSCAR
        write("POSCAR", atoms, format="vasp", vasp5=True, direct=True)
        
        # Cleanup? (Optional)
        # shutil.copy("OUTCAR", f"OUTCAR_step{step}")

if __name__ == "__main__":
    main()
"""

K_SPACING = 60.0

VASP_KPOINTS_TEMPLATE = """K-Points
0
Gamma
 {kx} {ky} {kz}
 0 0 0
"""

VASP_OPTCELL = """1 1 0
0 0 0
0 0 0
"""

VASP_SUB_SH = """#!/bin/bash
#BSUB -J Cu_{orientation}_{strain_pct:02d}
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH

# Create the run wrapper script used by python
echo '#!/bin/bash' > run_vasp.sh
echo "mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std > vasp_std.log" >> run_vasp.sh
chmod +x run_vasp.sh

echo "Starting VASP external relaxation loop..."
python3 relax_xy.py > relax_xy.log
"""

# LAMMPS settings
LAMMPS_EAM_TEMPLATE = """units metal
dimension 3
boundary p p p
atom_style atomic
read_data data.lammps_eam

mass 1 63.546

# Potential
pair_style eam/alloy
pair_coeff * * ../../lammps_eam/Cu_mishin1.eam.alloy Cu

# Relaxation (Fix Z, Relax XY)
# For strain 00, relax everything
{fix_cmd}
minimize 1.0e-10 1.0e-10 1000 10000

# Update thermo
run 0

# Calculate stress
variable pzz equal pzz
variable stress_zz equal -v_pzz/10000.0
print "ZZ_STRESS_GPA ${{stress_zz}}" file stress_eam.txt
"""

LAMMPS_SUS2_TEMPLATE = """units metal
dimension 3
boundary p p p
atom_style atomic
read_data data.lammps_sus2

mass 1 63.546

# Potential
pair_style mlip mlip.ini
pair_coeff * *

# Relaxation (Fix Z, Relax XY)
# For strain 00, relax everything
{fix_cmd}
minimize 1.0e-10 1.0e-10 1000 10000

# Update thermo
run 0

# Calculate stress
variable pzz equal pzz
variable stress_zz equal -v_pzz/10000.0
print "ZZ_STRESS_GPA ${{stress_zz}}" file stress_sus2.txt
"""

def is_orthorhombic(atoms, tol_deg=1e-3):
    ang = np.array(atoms.get_cell().angles(), dtype=float)
    return np.max(np.abs(ang - 90.0)) < tol_deg

def build_oriented_orth_cell(symbol, a0, x, y, z, tol_deg=1e-3):
    base = bulk(symbol, "fcc", a=a0, cubic=True)
    x = np.array(x, dtype=int)
    y = np.array(y, dtype=int)
    z = np.array(z, dtype=int)
    if (np.dot(x, y) != 0) or (np.dot(x, z) != 0) or (np.dot(y, z) != 0):
        raise ValueError(f"Directions not orthogonal:\n x={x}, y={y}, z={z}")
    P = np.column_stack([x, y, z]).astype(int)
    cand1 = make_supercell(base, P)
    if is_orthorhombic(cand1, tol_deg=tol_deg):
        cand1.set_cell(np.diag(cand1.cell.lengths()), scale_atoms=True)
        return cand1
    cand2 = make_supercell(base, P.T)
    if is_orthorhombic(cand2, tol_deg=tol_deg):
        cand2.set_cell(np.diag(cand2.cell.lengths()), scale_atoms=True)
        return cand2
    raise RuntimeError("Failed to build orthorhombic cell.")

def get_base_structure(orientation, a0):
    if orientation == '100':
        # x=[010], y=[001], z=[100] -> Z is [100]
        # Repeat to get ~100 atoms -> Changed to unit cell
        base = build_oriented_orth_cell('Cu', a0, [0,1,0], [0,0,1], [1,0,0])
        return base #.repeat((3,3,3)) # 4*27 = 108
    elif orientation == '110':
        # x=[1-10], y=[001], z=[110]
        base = build_oriented_orth_cell('Cu', a0, [1,-1,0], [0,0,1], [1,1,0])
        return base #.repeat((3,2,2)) # 8*12 = 96
    elif orientation == '111':
        # x=[1-10], y=[11-2], z=[111]
        base = build_oriented_orth_cell('Cu', a0, [1,-1,0], [1,1,-2], [1,1,1])
        return base #.repeat((2,1,2)) # 24*4 = 96
    return None

def get_kpoints(atoms, k_spacing=K_SPACING):
    cell = atoms.get_cell()
    lengths = np.linalg.norm(cell, axis=1)
    kpts = []
    for l in lengths:
        k = int(round(k_spacing / l))
        kpts.append(max(1, k))
    return kpts

def generate_structure_and_write(ori, strain, a0, filename, format):
    base_atoms = get_base_structure(ori, a0)
    atoms = base_atoms.copy()
    cell = atoms.get_cell()
    cell[2, 2] *= (1.0 + strain)
    atoms.set_cell(cell, scale_atoms=True)
    write(filename, atoms, format=format)
    return atoms

def main():
    # Removed argument parser, back to one-shot generation
    orientations = ['100', '110', '111']
    
    for ori in orientations:
        work_dir = f"work_{ori}"
        os.makedirs(work_dir, exist_ok=True)
        
        # Determine coupling for box/relax
        couple_style = "xy" if ori == '100' else "none"

        for strain in STRAINS:
            strain_pct = int(round(strain * 100))
            folder = os.path.join(work_dir, f"strain_{strain_pct:02d}")
            os.makedirs(folder, exist_ok=True)
            
            # 1. Generate VASP structure (POSCAR) using A0_VASP
            atoms_vasp = generate_structure_and_write(ori, strain, A0_VASP, 
                                        os.path.join(folder, "POSCAR"), "vasp")
            
            # 2. Generate EAM structure (data.lammps_eam) using A0_EAM
            generate_structure_and_write(ori, strain, A0_EAM, 
                                        os.path.join(folder, "data.lammps_eam"), "lammps-data")
            
            # 3. Generate SUS2 structure (data.lammps_sus2) using A0_SUS2
            generate_structure_and_write(ori, strain, A0_SUS2, 
                                        os.path.join(folder, "data.lammps_sus2"), "lammps-data")
            
            # Write VASP inputs
            with open(os.path.join(folder, "INCAR"), "w") as f:
                f.write(VASP_INCAR_TEMPLATE.format(orientation=ori, strain=strain))
            
            # Generate KPOINTS based on cell size (use VASP atoms)
            # Using KSPACING in INCAR, so KPOINTS file is not strictly needed or can be minimal/automatic
            # But since we provide KSPACING, VASP ignores KPOINTS usually or we can skip writing it.
            # To be safe, we can write a dummy KPOINTS or remove it.
            # However, if we write KPOINTS, VASP might use it instead of KSPACING depending on version.
            # Standard practice: If KSPACING is present, KPOINTS file is optional or should be 'Automatic'.
            # Let's keep writing KPOINTS for compatibility but note KSPACING takes precedence in newer VASP if KPOINTS is not found?
            # Actually, if KSPACING is set, VASP generates KPOINTS internally. We should probably DELETE KPOINTS file if it exists
            # or just not write it.
            kpoints_path = os.path.join(folder, "KPOINTS")
            if os.path.exists(kpoints_path):
                os.remove(kpoints_path)
            # kx, ky, kz = get_kpoints(atoms_vasp)
            # with open(os.path.join(folder, "KPOINTS"), "w") as f:
            #     f.write(VASP_KPOINTS_TEMPLATE.format(kx=kx, ky=ky, kz=kz))
            
            with open(os.path.join(folder, "OPTCELL"), "w") as f:
                f.write(VASP_OPTCELL)
            
            # Write sub.sh
            if strain_pct == 0:
                 with open(os.path.join(folder, "sub.sh"), "w") as f:
                    f.write(VASP_SUB_SH_RELAX_00.format(orientation=ori))
            else:
                with open(os.path.join(folder, "sub.sh"), "w") as f:
                    f.write(VASP_SUB_SH_STRAIN.format(orientation=ori, strain_pct=strain_pct))
                # relax_xy.py removed

            # Copy POTCAR if exists
            potcar_src = os.path.abspath(os.path.join("vasp", "POTCAR"))
            if os.path.exists(potcar_src):
                shutil.copy(potcar_src, os.path.join(folder, "POTCAR"))
            
            # EAM Input
            if strain_pct == 0:
                fix_cmd = "fix 1 all box/relax aniso 0.0"
            else:
                fix_cmd = f"fix 1 all box/relax x 0.0 y 0.0 couple {couple_style}"

            with open(os.path.join(folder, "in.eam"), "w") as f:
                f.write(LAMMPS_EAM_TEMPLATE.format(fix_cmd=fix_cmd))
                
            # SUS2 Input
            with open(os.path.join(folder, "in.sus2"), "w") as f:
                f.write(LAMMPS_SUS2_TEMPLATE.format(fix_cmd=fix_cmd))
            
            # Copy MLIP files
            mlip_src = os.path.abspath(os.path.join("lammps_sus2", "mlip.ini"))
            mtp_src = os.path.abspath(os.path.join("lammps_sus2", "current.mtp"))
            if os.path.exists(mlip_src):
                shutil.copy(mlip_src, os.path.join(folder, "mlip.ini"))
            if os.path.exists(mtp_src):
                shutil.copy(mtp_src, os.path.join(folder, "current.mtp"))

    print("All structures and inputs generated in work_100/, work_110/, work_111/.")

if __name__ == "__main__":
    main()
