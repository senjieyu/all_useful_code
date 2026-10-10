import os
import sys
import shutil
import json
from ase.io import read, write

def load_config():
    # Prioritize local config.json
    if os.path.exists("config.json"):
        print(f"Loading config from {os.path.abspath('config.json')}")
        with open("config.json", 'r') as f: return json.load(f)
    
    # Fallback to parent config.json
    if os.path.exists("../config.json"):
        print(f"Loading config from {os.path.abspath('../config.json')}")
        with open("../config.json", 'r') as f: return json.load(f)
        
    return None

def main():
    if len(sys.argv) < 2:
        print("Usage: python setup_lammps.py <input.xyz>")
        sys.exit(1)
        
    input_xyz = sys.argv[1]
    config = load_config()
    if not config:
        print("Error: config.json not found in current or parent directory.")
        sys.exit(1)
        
    # Find config file path to copy (to put in struct folders)
    # Must match the one we loaded to ensure consistency
    config_path = None
    if os.path.exists("config.json"):
        config_path = "config.json"
    elif os.path.exists("../config.json"):
        config_path = "../config.json"
        
    # Paths
    cwd = os.getcwd()
    
    # Locate potentials
    if "potential_files" not in config["lammps"]:
        print("Error: 'potential_files' missing in config['lammps'].")
        sys.exit(1)
        
    pot_files = config["lammps"]["potential_files"]
    eam_file = pot_files.get("eam_file")
    
    if not eam_file:
         print("Error: 'eam_file' not specified in config['lammps']['potential_files'].")
         print("Please ensure you are using the correct config.json (local one for EAM).")
         sys.exit(1)
    
    def find_file(fname):
        if os.path.exists(fname): return os.path.abspath(fname)
        # Try parent of setup script (Unified_Workflow/LAMMPS/..)
        p1 = os.path.join(os.path.dirname(cwd), fname)
        if os.path.exists(p1): return p1
        # Try grandparent
        p2 = os.path.join(os.path.dirname(os.path.dirname(cwd)), fname)
        if os.path.exists(p2): return p2
        return None
        
    abs_eam = find_file(eam_file)
    
    if not abs_eam:
        print(f"Error: Could not find {eam_file}.")
        sys.exit(1)
        
    print(f"Using potential: {abs_eam}")
    
    # Run Script
    run_script = os.path.abspath("run_lammps.py")
    if not os.path.exists(run_script):
        print("Error: run_lammps.py not found in current directory.")
        sys.exit(1)
        
    # Read XYZ
    atoms_list = read(input_xyz, index=':')
    print(f"Found {len(atoms_list)} structures.")
    
    # Find style.txt to copy (optional but recommended)
    style_path = "style.txt"
    if os.path.exists("../style.txt"):
        style_path = "../style.txt"
    elif os.path.exists("style.txt"):
        style_path = "style.txt"
    else:
        style_path = None

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
        shutil.copy(abs_eam, os.path.join(folder, os.path.basename(eam_file)))
        shutil.copy(run_script, os.path.join(folder, "run_lammps.py"))
        if config_path:
            shutil.copy(config_path, os.path.join(folder, "config.json"))
        if style_path:
            shutil.copy(style_path, os.path.join(folder, "style.txt"))
        
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
