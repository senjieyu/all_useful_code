import os
import sys
import shutil
import json

def load_config():
    if os.path.exists("../config.json"):
        with open("../config.json", 'r') as f: return json.load(f)
    if os.path.exists("config.json"):
        with open("config.json", 'r') as f: return json.load(f)
    return None

def main():
    # 1. Load Config
    config = load_config()
    if not config:
        print("Error: config.json not found.")
        sys.exit(1)
        
    queue = config["queue"]
    
    # 2. Prepare submission scripts directory
    scripts_dir = "submission_scripts_elastic"
    os.makedirs(scripts_dir, exist_ok=True)
    
    # 3. Define Submission Template
    submit_template = f"""#!/bin/bash
#BSUB -J elastic_workflow_{{idx}}
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

# Run Elastic Workflow
python3 run_elastic.py

echo "End: $(date)"
"""

    # 4. Find all struct_X folders
    folders = [f for f in os.listdir(".") if os.path.isdir(f) and f.startswith("struct_")]
    # Sort them by index
    try:
        folders.sort(key=lambda x: int(x.split('_')[1]))
    except:
        folders.sort()
        
    print(f"Found {len(folders)} structure folders.")
    
    # 5. Process each folder
    count = 0
    for folder in folders:
        try:
            idx = int(folder.split('_')[1])
        except:
            idx = 999
            
        print(f"Setting up elastic submission for {folder}...")
        
        # Copy run_elastic.py
        if os.path.exists("run_elastic.py"):
            shutil.copy("run_elastic.py", os.path.join(folder, "run_elastic.py"))
        else:
            print("Error: run_elastic.py not found in current directory.")
            sys.exit(1)
            
        # Write submit_elastic.sh
        with open(os.path.join(folder, "submit_elastic.sh"), 'w', newline='\n') as f:
            f.write(submit_template.format(idx=idx))
            
        # Write wrapper script
        wrapper_name = f"submit_elastic_{folder}.sh"
        wrapper_path = os.path.join(scripts_dir, wrapper_name)
        abs_folder = os.path.abspath(folder).replace(os.sep, '/')
        
        with open(wrapper_path, 'w', newline='\n') as f:
            f.write("#!/bin/bash\n")
            f.write(f"cd \"{abs_folder}\"\n")
            f.write("bsub < submit_elastic.sh\n")
            f.write(f"echo \"Submitted Elastic Calculation for {folder}\"\n")
            
        count += 1
        
    print(f"Successfully set up {count} elastic calculation submissions.")
    print(f"To submit all, run: python submit_batch.py {scripts_dir}")

if __name__ == "__main__":
    main()