import os
import sys
import shutil
import json
from ase.io import read, write

def load_config():
    if os.path.exists("../config.json"):
        with open("../config.json", 'r') as f: return json.load(f)
    if os.path.exists("config.json"):
        with open("config.json", 'r') as f: return json.load(f)
    return None

def main():
    if len(sys.argv) < 2:
        print("Usage: python setup_lammps.py <input.xyz>")
        sys.exit(1)
        
    input_xyz = sys.argv[1]
    config = load_config()
    if not config:
        print("Error: ../config.json not found.")
        sys.exit(1)
        
    # Find config file path to copy
    config_path = "config.json"
    if os.path.exists("../config.json"):
        config_path = "../config.json"
    elif os.path.exists("config.json"):
        config_path = "config.json"
    else:
        # Should be covered by load_config check, but for safety
        config_path = None
        
    # Paths
    cwd = os.getcwd()
    # Locate potentials
    # Config has filenames, but we need absolute paths to copy them
    # Assuming they are in the directory where setup is run, or user provided full path
    # Actually, let's search in parent directories if not found
    pot_files = config["lammps"]["potential_files"]
    mlip_ini = pot_files["mlip_ini"]
    mtp_file = pot_files["mtp_file"]
    
    def find_file(fname):
        if os.path.exists(fname): return os.path.abspath(fname)
        # Try parent of setup script (Unified_Workflow/LAMMPS/..)
        p1 = os.path.join(os.path.dirname(cwd), fname)
        if os.path.exists(p1): return p1
        # Try grandparent
        p2 = os.path.join(os.path.dirname(os.path.dirname(cwd)), fname)
        if os.path.exists(p2): return p2
        # Try specific known location from previous context
        p3 = os.path.join(os.path.dirname(os.path.dirname(cwd)), "EOS_2_LAMMPS", fname)
        if os.path.exists(p3): return p3
        return None
        
    abs_mlip = find_file(mlip_ini)
    abs_mtp = find_file(mtp_file)
    
    if not abs_mlip or not abs_mtp:
        print(f"Error: Could not find {mlip_ini} or {mtp_file}.")
        sys.exit(1)
        
    print(f"Using potentials: {abs_mlip}, {abs_mtp}")
    
    # Run Script
    run_script = os.path.abspath("run_lammps.py")
    if not os.path.exists(run_script):
        print("Error: run_lammps.py not found in current directory.")
        sys.exit(1)
        
    # Read XYZ
    atoms_list = read(input_xyz, index=':')
    print(f"Found {len(atoms_list)} structures.")
    
    # Queue settings
    queue = config["queue"]
    
    submit_template = f"""#!/bin/bash
#BSUB -J lmp_workflow_{{idx}}
#BSUB -q {queue['queue_name']}
#BSUB -n {queue['nproc']}
#BSUB -R "span[ptile={queue['ptile']}]"
#BSUB -o %J.out
#BSUB -e %J.err

set -euo pipefail

# Environment
export I_MPI_HYDRA_BOOTSTRAP=lsf
# Add any module loads here if needed

echo "Host: $(hostname)"
echo "Start: $(date)"

# Run Workflow
python run_lammps.py

echo "End: $(date)"
"""

    scripts_dir = "submission_scripts"
    os.makedirs(scripts_dir, exist_ok=True)
    
    for i, atoms in enumerate(atoms_list):
        folder = f"struct_{i}"
        os.makedirs(folder, exist_ok=True)
        
        # Write start.lmp
        # Use atomic style and ensure elements match config
        elements = config["elements"]
        if len(elements) == 1:
            atoms.symbols = [elements[0]] * len(atoms)
            
        write(os.path.join(folder, "start.lmp"), atoms, format='lammps-data', atom_style='atomic', specorder=elements)
        
        # Copy files
        shutil.copy(abs_mlip, os.path.join(folder, "mlip.ini"))
        shutil.copy(abs_mtp, os.path.join(folder, "current.mtp"))
        shutil.copy(run_script, os.path.join(folder, "run_lammps.py"))
        if config_path:
            shutil.copy(config_path, os.path.join(folder, "config.json"))
        
        # Write submit.sh
        with open(os.path.join(folder, "submit.sh"), 'w') as f:
            f.write(submit_template.format(idx=i))
            
        # Write wrapper
        wrapper = os.path.join(scripts_dir, f"submit_{folder}.sh")
        abs_folder = os.path.abspath(folder).replace(os.sep, '/')
        with open(wrapper, 'w', newline='\n') as f:
            f.write("#!/bin/bash\n")
            f.write(f"cd \"{abs_folder}\"\n")
            f.write("bsub < submit.sh\n")
            f.write(f"echo \"Submitted {folder}\"\n")
            
    print(f"Setup complete. Created {len(atoms_list)} folders.")
    print(f"To submit: python submit_batch.py {scripts_dir} (You need to copy submit_batch.py first or run the shell scripts manually)")

if __name__ == "__main__":
    main()
