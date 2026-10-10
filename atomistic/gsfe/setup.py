import os
import glob
import shutil
import numpy as np

def main():
    stru_files = sorted(glob.glob("stru/*.vasp"))
    if not stru_files:
        print("No .vasp files found in stru/ directory.")
        return

    # Create VASP submission scripts
    for vasp_file in stru_files:
        basename = os.path.basename(vasp_file)
        disp_str = basename.split("disp")[-1].replace(".vasp", "")
        
        vasp_dir = os.path.join("vasp", f"calc_{disp_str}")
        os.makedirs(vasp_dir, exist_ok=True)
        shutil.copy(vasp_file, os.path.join(vasp_dir, "POSCAR"))
        
        if os.path.exists("vasp/POTCAR"):
            shutil.copy("vasp/POTCAR", os.path.join(vasp_dir, "POTCAR"))
        if os.path.exists("vasp/INCAR_relax"):
            shutil.copy("vasp/INCAR_relax", os.path.join(vasp_dir, "INCAR"))
        if os.path.exists("vasp/INCAR_single"):
            shutil.copy("vasp/INCAR_single", os.path.join(vasp_dir, "INCAR_single"))
        
        if os.path.exists("vasp/sub.sh"):
            with open("vasp/sub.sh", "r") as f:
                sub_lines = f.read()
            
            vasp_script = sub_lines + f"""
echo "Starting VASP relaxation..."
$COMMAND_std > vasp_relax.log
cp CONTCAR POSCAR
cp INCAR_single INCAR
echo "Starting VASP single point..."
$COMMAND_std > vasp_single.log
"""
            with open(os.path.join(vasp_dir, "sub.sh"), "w") as f:
                f.write(vasp_script)

    # Setup LAMMPS scripts
    # Calculate exact Z boundaries for 24 layers
    a = 3.615
    z_layer = a / np.sqrt(3)
    z_mid = 11.5 * z_layer
    z_top = 24.5 * z_layer

    in_eam = f"""# LAMMPS EAM GSFE Calculation
label SFE
variable x_disp index 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0

units metal
dimension 3
boundary p p f
atom_style atomic

variable latparam1 equal 3.615
variable xdim equal ${{latparam1}}*sqrt(6)/2*10
variable ydim equal ${{latparam1}}*sqrt(2)/2*10

lattice fcc ${{latparam1}} orient x 1 1 2 orient y -1 1 0 orient z -1 -1 1

# EXACT 24 LAYERS: 12 bottom, 12 top (vacuum added)
# Box from -5 to 15 in Z.
# Fixed boundary 'f' prevents lost atoms better than 's' for high energy relaxation,
# providing we have enough vacuum (which we do: -5 to 1.9 and 10.1 to 15).
region whole block 0 5 0 5 -5 15 units lattice
create_box 1 whole
mass 1 63.546

# Shift regions up by 2 units to be safely inside
region 1 block -0.1 5.1 -0.1 5.1 1.9 5.9 units lattice
create_atoms 1 region 1

region 2 block -0.1 5.1 -0.1 5.1 5.9 10.1 units lattice
create_atoms 1 region 2

pair_style eam/alloy
pair_coeff * * Cu_mishin1.eam.alloy Cu

group bot region 1
group top region 2

compute peratom all pe/atom
compute eatoms all reduce sum c_peratom

min_style cg
minimize 1e-10 1e-10 1000 1000

thermo 1
thermo_style custom step pe c_eatoms
run 0
variable E equal "c_eatoms"
variable Eo equal ${{E}}

variable A equal (${{xdim}}*${{ydim}})*1e-20
variable Cf equal 1.60217657e-16

displace_atoms top move -${{x_disp}} 0.0 0.0 units box
fix 1 all setforce 0.0 0.0 NULL

min_style cg
minimize 1e-10 1e-10 1000 1000

thermo 1
thermo_style custom step pe c_eatoms
run 0

variable Ef equal "c_eatoms"
variable SFE equal ((${{Ef}}-${{Eo}})*${{Cf}})/(${{A}})

print "DISP: ${{x_disp}} SFE: ${{SFE}} mJ/m^2"

clear
next x_disp
jump in.eam SFE
"""
    eam_dir = "lammps_eam"
    os.makedirs(eam_dir, exist_ok=True)
    with open(os.path.join(eam_dir, "in.eam"), "w") as f:
        f.write(in_eam)

    in_sus2 = f"""# LAMMPS SUS2 GSFE Calculation
label SFE
variable x_disp index 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0

units metal
dimension 3
boundary p p f
atom_style atomic

variable latparam1 equal 3.615
variable xdim equal ${{latparam1}}*sqrt(6)/2*10
variable ydim equal ${{latparam1}}*sqrt(2)/2*10

lattice fcc ${{latparam1}} orient x 1 1 2 orient y -1 1 0 orient z -1 -1 1

# EXACT 24 LAYERS: 12 bottom, 12 top (vacuum added)
# Box from -5 to 15 in Z.
# Fixed boundary 'f' prevents lost atoms better than 's' for high energy relaxation,
# providing we have enough vacuum (which we do: -5 to 1.9 and 10.1 to 15).
region whole block 0 5 0 5 -5 15 units lattice
create_box 1 whole
mass 1 63.546

# Shift regions up by 2 units to be safely inside
region 1 block -0.1 5.1 -0.1 5.1 1.9 5.9 units lattice
create_atoms 1 region 1

region 2 block -0.1 5.1 -0.1 5.1 5.9 10.1 units lattice
create_atoms 1 region 2

pair_style mlip mlip.ini
pair_coeff * *

group bot region 1
group top region 2

compute peratom all pe/atom
compute eatoms all reduce sum c_peratom

min_style cg
minimize 1e-10 1e-10 1000 1000

thermo 1
thermo_style custom step pe c_eatoms
run 0
variable E equal "c_eatoms"
variable Eo equal ${{E}}

variable A equal (${{xdim}}*${{ydim}})*1e-20
variable Cf equal 1.60217657e-16

displace_atoms top move -${{x_disp}} 0.0 0.0 units box
fix 1 all setforce 0.0 0.0 NULL

min_style cg
minimize 1e-10 1e-10 1000 1000

thermo 1
thermo_style custom step pe c_eatoms
run 0

variable Ef equal "c_eatoms"
variable SFE equal ((${{Ef}}-${{Eo}})*${{Cf}})/(${{A}})

print "DISP: ${{x_disp}} SFE: ${{SFE}} mJ/m^2"

clear
next x_disp
jump in.sus2 SFE
"""
    sus2_dir = "lammps_sus2"
    os.makedirs(sus2_dir, exist_ok=True)
    with open(os.path.join(sus2_dir, "in.sus2"), "w") as f:
        f.write(in_sus2)

    # Note about potential file location fix
    eam_script_path = os.path.join(eam_dir, "in.eam")
    with open(eam_script_path, 'r') as f:
        content = f.read()
    content = content.replace("../Cu_mishin1.eam.alloy", "Cu_mishin1.eam.alloy")
    with open(eam_script_path, 'w') as f:
        f.write(content)

    print("Setup completed. VASP uses individual job folders. LAMMPS uses a single in-script loop (no calc_* folders).")

if __name__ == "__main__":
    main()
