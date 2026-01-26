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
    
    if not os.path.exists(incar_relax) or not os.path.exists(incar_elastic):
        # Fallback to current dir
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
    
    # VASP submit template (LSF)
    submit_template = """#!/bin/bash
#BSUB -J vasp_elas_{idx}
#BSUB -q 33
#BSUB -n 40
#BSUB -R "span[ptile=40]"
#BSUB -o %J.out
#BSUB -e %J.err

set -euo pipefail

# Environment
# Adjust these paths for your cluster
export I_MPI_HYDRA_BOOTSTRAP=lsf
module load vasp/5.4.4 2>/dev/null || true
# Or export PATH directly if module not available
# export PATH=/path/to/vasp_std:$PATH

echo "Host: $(hostname)"
echo "Start: $(date)"

# 1. Relaxation Step
if [ ! -f "RELAX_DONE" ]; then
    echo "Starting Relaxation..."
    cp INCAR_relax INCAR
    # Ensure POTCAR exists (user must provide POTCAR in setup or link it)
    if [ ! -f "POTCAR" ]; then
        echo "Error: POTCAR missing!"
        exit 1
    fi
    
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

    # POTCAR warning
    print("WARNING: You must manually place or link a POTCAR file into the output folders,")
    print("         or modify this script to copy a POTCAR from a known location.")
    
    for i, atoms in enumerate(frames):
        folder_name = f"struct_{i}"
        work_dir = os.path.join(output_folder, folder_name)
        os.makedirs(work_dir, exist_ok=True)
        
        # 1. Write POSCAR
        write(os.path.join(work_dir, "POSCAR"), atoms, format='vasp', direct=True, sort=True)
        
        # 2. Copy INCARs
        shutil.copy(incar_relax, os.path.join(work_dir, "INCAR_relax"))
        shutil.copy(incar_elastic, os.path.join(work_dir, "INCAR_elastic"))
        
        # 3. Generate POTCAR using vaspkit (User confirmed they have license/vaspkit installed)
        # Using the logic provided by user: (echo 103) | vaspkit
        try:
            cwd = os.getcwd()
            os.chdir(work_dir)
            # Suppress output to keep terminal clean
            os.system(r'(echo 103) | vaspkit > /dev/null 2>&1')
            os.chdir(cwd)
        except Exception as e:
            print(f"Warning: Failed to run vaspkit in {work_dir}: {e}")

        # 4. Write submit.sh
        submit_content = submit_template.format(idx=i)
        with open(os.path.join(work_dir, "submit.sh"), 'w') as f:
            f.write(submit_content)
            
        # 4. Write wrapper script
        wrapper_path = os.path.join(scripts_dir, f"submit_{folder_name}.sh")
        abs_work_dir = work_dir.replace(os.sep, '/') 
        
        with open(wrapper_path, 'w', newline='\n') as f:
            f.write("#!/bin/bash\n")
            f.write(f"cd \"{abs_work_dir}\"\n")
            f.write("bsub < submit.sh\n")
            f.write(f"echo \"Submitted vasp_elas_{i}\"\n")
            
    print(f"Setup complete. {len(frames)} folders created in {output_folder}.")
    print(f"Submission scripts in {scripts_dir}.")
    print("IMPORTANT: Copy your POTCAR to all subfolders before submitting!")
    print("Example: for d in " + output_folder + "/struct_*; do cp /path/to/POTCAR $d/; done")

if __name__ == "__main__":
    main()
