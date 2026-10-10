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
        print("Usage: python setup_vasp.py <input.xyz>")
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
        config_path = None

    # Queue settings
    queue = config["queue"]
    
    # Get vaspkit path from config or use default
    vaspkit_bin = config.get("vaspkit", {}).get("binary", "vaspkit")
    # Expand user path (~) if necessary
    vaspkit_bin = os.path.expanduser(vaspkit_bin)
    
    # Read XYZ
    print(f"Reading {input_xyz}...")
    atoms_list = read(input_xyz, index=':')
    print(f"Found {len(atoms_list)} structures.")
    
    submit_template = f"""#!/bin/bash
#BSUB -J vasp_workflow_{{idx}}
#BSUB -q {queue['queue_name']}
#BSUB -n {queue['nproc']}
#BSUB -R "span[ptile={queue['ptile']}]"
#BSUB -o %J.out
#BSUB -e %J.err

set -euo pipefail

# Environment
export I_MPI_HYDRA_BOOTSTRAP=lsf
{config.get("vasp", {}).get("module_load", "# module load vasp")}
export PATH=$PATH:{os.path.dirname(config.get("vasp", {}).get("binary", ""))}

echo "Host: $(hostname)"
echo "Start: $(date)"

# Run Workflow
python3 run_vasp.py

echo "End: $(date)"
"""

    scripts_dir = "submission_scripts"
    os.makedirs(scripts_dir, exist_ok=True)
    
    files_to_copy = ["run_vasp.py", "INCAR_relax", "INCAR_eos", "INCAR_elastic", "POTCAR"]
    # Check if they exist
    existing_files = []
    for f in files_to_copy:
        if os.path.exists(f):
            existing_files.append(f)
        else:
            # Special handling for POTCAR: might not exist if user wants vaspkit to generate it
            if f == "POTCAR":
                print(f"Warning: {f} not found. Ensure vaspkit is configured or POTCAR is generated later.")
                continue
            print(f"Error: {f} not found in current directory.")
            sys.exit(1)
            
    for i, atoms in enumerate(atoms_list):
        folder = f"struct_{i}"
        os.makedirs(folder, exist_ok=True)
        
        # Write POSCAR (initial)
        write(os.path.join(folder, "POSCAR"), atoms, format='vasp', direct=True, sort=True)
        
        # Copy files
        for f in existing_files:
            shutil.copy(f, os.path.join(folder, f))
            
        if config_path:
            shutil.copy(config_path, os.path.join(folder, "config.json"))

        # If POTCAR is missing, try to generate it with vaspkit
        if "POTCAR" not in existing_files:
            print(f"Generating POTCAR for {folder} using vaspkit...")
            try:
                import subprocess
                # vaspkit requires POSCAR in the directory
                # We run vaspkit inside the folder
                # Command: echo 103 | /path/to/vaspkit
                cmd = f"echo 103 | {vaspkit_bin}"
                # Use shell=True for piping
                process = subprocess.Popen(cmd, shell=True, cwd=folder, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = process.communicate()
                
                if process.returncode != 0:
                    print(f"  Warning: vaspkit failed with code {process.returncode}")
                    print(f"  Error output: {err.decode().strip()}")
                elif os.path.exists(os.path.join(folder, "POTCAR")):
                    print(f"  POTCAR generated successfully.")
                else:
                    print(f"  Warning: vaspkit ran but POTCAR not found.")
                    print(f"  Output: {out.decode().strip()}")
            except Exception as e:
                print(f"  Warning: Failed to run vaspkit: {e}")
            
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
    print(f"To submit: python submit_batch.py {scripts_dir}")

if __name__ == "__main__":
    main()
