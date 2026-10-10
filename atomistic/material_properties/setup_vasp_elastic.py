import os
import sys
import shutil
from ase.io import read, write

def main():
    if len(sys.argv) < 3:
        print("Usage: python setup_vasp_elastic.py <input.xyz> <output_folder>")
        sys.exit(1)
        
    input_xyz = sys.argv[1]
    output_folder = sys.argv[2]
    
    # Absolute paths
    cwd = os.getcwd()
    input_xyz = os.path.abspath(input_xyz)
    output_folder = os.path.abspath(output_folder)
    
    # Check dependencies in script dir
    script_dir = os.path.dirname(os.path.abspath(__file__))
    incar_relax = os.path.join(script_dir, "INCAR_relax")
    incar_elastic = os.path.join(script_dir, "INCAR_elastic")
    
    # Fallback checks for INCAR files
    if not os.path.exists(incar_relax) or not os.path.exists(incar_elastic):
        if os.path.exists("INCAR_relax"): incar_relax = "INCAR_relax"
        if os.path.exists("INCAR_elastic"): incar_elastic = "INCAR_elastic"
        
    if not os.path.exists(incar_relax) or not os.path.exists(incar_elastic):
        print("Error: INCAR_relax or INCAR_elastic not found.")
        sys.exit(1)
            
    # Read XYZ frames
    print(f"Reading {input_xyz}...")
    frames = read(input_xyz, index=':')
    print(f"Found {len(frames)} structures.")
    
    # Create base output folder
    os.makedirs(output_folder, exist_ok=True)
    
    # Submission scripts folder
    scripts_dir = os.path.join(output_folder, "submission_scripts")
    os.makedirs(scripts_dir, exist_ok=True)
    
    # ==========================================
    # FIXED SUBMIT TEMPLATE (Corrected Environment)
    # ==========================================
    submit_template = """#!/bin/bash
#BSUB -J vasp_elas_{idx}
#BSUB -q 33
#BSUB -n 40
#BSUB -R "span[ptile=40]"
#BSUB -o %J.out
#BSUB -e %J.err

set -euo pipefail

# Environment Setup
export I_MPI_HYDRA_BOOTSTRAP=lsf

# 1. Load Intel Compiler and MPI (Matches your .bashrc)
module purge
module load mpi/2021.6.0 compiler/2022.1.0 mkl/2022.2.0

# 2. Set VASP 5.4.4 Path (Manually specified)
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH

echo "Host: $(hostname)"
echo "Start: $(date)"

# Debug: check if vasp_std is found
which vasp_std || echo "WARNING: vasp_std not found in PATH"

# 1. Relaxation Step
if [ ! -f "RELAX_DONE" ]; then
    echo "Starting Relaxation..."
    cp INCAR_relax INCAR
    
    # Check if POTCAR exists (Vaspkit should have generated it)
    if [ ! -f "POTCAR" ]; then
        echo "Error: POTCAR missing! Vaspkit generation failed or file not copied."
        exit 1
    fi
    
    # Run vasp_std
    mpirun -np 40 vasp_std > vasp_relax.log
    
    cp CONTCAR POSCAR
    touch RELAX_DONE
    echo "Relaxation Done."
fi

# 2. Elastic Calculation Step
echo "Starting Elastic Calculation..."
cp INCAR_elastic INCAR
mpirun -np 40 vasp_std > vasp_elastic.log

echo "End: $(date)"
"""

    for i, atoms in enumerate(frames):
        folder_name = f"struct_{i}"
        work_dir = os.path.join(output_folder, folder_name)
        os.makedirs(work_dir, exist_ok=True)
        
        # 1. Write POSCAR
        write(os.path.join(work_dir, "POSCAR"), atoms, format='vasp', direct=True, sort=True)
        
        # 2. Copy INCARs
        shutil.copy(incar_relax, os.path.join(work_dir, "INCAR_relax"))
        shutil.copy(incar_elastic, os.path.join(work_dir, "INCAR_elastic"))
        
        # 3. Generate POTCAR using vaspkit
        try:
            cwd_temp = os.getcwd()
            os.chdir(work_dir)
            # Run vaspkit option 103 (Generate POTCAR)
            # Redirect output to prevent clutter
            ret = os.system(r'(echo 103) | vaspkit > /dev/null 2>&1')
            if ret != 0:
                print(f"Warning: vaspkit returned error code in {folder_name}")
            os.chdir(cwd_temp)
        except Exception as e:
            print(f"Warning: Failed to run vaspkit in {work_dir}: {e}")

        # 4. Write submit.sh
        submit_content = submit_template.format(idx=i)
        with open(os.path.join(work_dir, "submit.sh"), 'w') as f:
            f.write(submit_content)
            
        # 5. Write wrapper script (for easy bulk submission)
        wrapper_path = os.path.join(scripts_dir, f"submit_{folder_name}.sh")
        abs_work_dir = work_dir.replace(os.sep, '/') 
        
        with open(wrapper_path, 'w', newline='\n') as f:
            f.write("#!/bin/bash\n")
            f.write(f"cd \"{abs_work_dir}\"\n")
            f.write("bsub < submit.sh\n")
            f.write(f"echo \"Submitted vasp_elas_{i}\"\n")
            
    print(f"Setup complete. {len(frames)} folders created in {output_folder}.")
    print(f"Submission scripts in {scripts_dir}.")
    print("NOTE: Vaspkit was run to generate POTCARs.")
    print("      Please verify that 'POTCAR' exists in the subfolders before submitting.")

if __name__ == "__main__":
    main()

