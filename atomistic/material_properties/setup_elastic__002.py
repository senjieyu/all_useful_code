import os
import sys
import shutil
import glob
from ase.io import read, write

def main():
    if len(sys.argv) < 3:
        print("Usage: python setup_elastic.py <input.xyz> <output_folder>")
        sys.exit(1)
        
    input_xyz = sys.argv[1]
    output_folder = sys.argv[2]
    
    # Absolute paths
    cwd = os.getcwd()
    input_xyz = os.path.abspath(input_xyz)
    output_folder = os.path.abspath(output_folder)
    
    # Check dependencies
    # We look for files in the same directory as this script (if it's in the bundle)
    # or in the current working directory.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Files to look for: mlip.ini, current.mtp, calc_elastic_single.py
    # We prioritize files in script_dir, then CWD.
    
    def get_path(fname):
        p1 = os.path.join(script_dir, fname)
        if os.path.exists(p1): return p1
        p2 = os.path.join(cwd, fname)
        if os.path.exists(p2): return p2
        return None

    mlip_path = get_path("mlip.ini")
    mtp_path = get_path("current.mtp")
    calc_path = get_path("calc_elastic_single.py")

    if not mlip_path or not mtp_path or not calc_path:
        print("Error: Missing dependencies (mlip.ini, current.mtp, or calc_elastic_single.py).")
        print(f"Searched in {script_dir} and {cwd}")
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
    
    # Template for submission
    # Using the same LMP path as before
    LMP_BIN = "/work/phy-luoyf/apps/lammps_sus2/interface-lammps-mlip-2/lmp_intel_cpu_intelmpi"
    
    submit_template = """#!/bin/bash
#BSUB -J elastic_{idx}
#BSUB -q 33
#BSUB -n 1
#BSUB -R "span[ptile=1]"
#BSUB -o %J.out
#BSUB -e %J.err

set -euo pipefail

# Environment
export I_MPI_HYDRA_BOOTSTRAP=lsf
LMP_BIN="{lmp_bin}"

echo "Host: $(hostname)"
echo "Start: $(date)"

# Run calculation
# We use python to drive the calculation.
# It requires ase and numpy.
python calc_elastic_single.py start.lmp "$LMP_BIN"

echo "End: $(date)"
"""

    for i, atoms in enumerate(frames):
        folder_name = f"struct_{i}"
        work_dir = os.path.join(output_folder, folder_name)
        os.makedirs(work_dir, exist_ok=True)
        
        # 1. Write structure file
        write(os.path.join(work_dir, "start.lmp"), atoms, format='lammps-data', atom_style='atomic', specorder=['Cu'])
        
        # 2. Copy scripts and potentials
        shutil.copy(mlip_path, work_dir)
        shutil.copy(mtp_path, work_dir)
        shutil.copy(calc_path, work_dir)
        
        # 3. Write submit.sh
        submit_content = submit_template.format(idx=i, lmp_bin=LMP_BIN)
        with open(os.path.join(work_dir, "submit.sh"), 'w') as f:
            f.write(submit_content)
            
        # 4. Write wrapper script
        wrapper_path = os.path.join(scripts_dir, f"submit_{folder_name}.sh")
        abs_work_dir = work_dir.replace(os.sep, '/') # Posix path for bash
        
        with open(wrapper_path, 'w', newline='\n') as f:
            f.write("#!/bin/bash\n")
            f.write(f"cd \"{abs_work_dir}\"\n")
            f.write("bsub < submit.sh\n")
            f.write(f"echo \"Submitted elastic_{i}\"\n")
            
    print(f"Setup complete. {len(frames)} folders created in {output_folder}.")
    print(f"Submission scripts in {scripts_dir}.")
    print("Run 'python submit_batch.py " + scripts_dir + "' to submit.")

if __name__ == "__main__":
    main()
