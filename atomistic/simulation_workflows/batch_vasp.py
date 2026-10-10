import os
import glob
import math

# Configuration
VASP_BIN_PATH = "/work/phy-huangj/app/vasp.5.4.4/bin"
TASKS_PER_NODE = 3
CORES_PER_NODE = 40

def submit_vasp_jobs():
    # Find all strain directories under work_* folders
    strain_dirs = sorted([d for d in glob.glob("work_*/strain_*") if os.path.isdir(d)])
    
    if not strain_dirs:
        print("No strain directories found! (Did you run 'python stru/build.py'?)")
        return

    # Chunk directories
    chunks = [strain_dirs[i:i + TASKS_PER_NODE] for i in range(0, len(strain_dirs), TASKS_PER_NODE)]
    
    print(f"Found {len(strain_dirs)} directories. Creating {len(chunks)} packed jobs.")

    for i, chunk in enumerate(chunks):
        job_name = f"vasp_pack_{i}"
        script_name = f"sub_vasp_pack_{i}.sh"
        
        with open(script_name, 'w', newline='\n') as f:
            f.write("#!/bin/bash\n")
            f.write(f"#BSUB -J {job_name}\n")
            f.write("#BSUB -q 33\n")
            f.write(f"#BSUB -n {CORES_PER_NODE}\n")
            f.write("#BSUB -e %J.err\n")
            f.write("#BSUB -o %J.out\n")
            f.write(f"#BSUB -R \"span[ptile={CORES_PER_NODE}]\"\n\n")
            
            f.write("hostfile=`echo $LSB_DJOB_HOSTFILE`\n")
            f.write("NP=`cat $hostfile | wc -l`\n")
            f.write(f"export PATH={VASP_BIN_PATH}:$PATH\n")
            f.write("cd $LS_SUBCWD\n")
            f.write("BASE_DIR=$LS_SUBCWD\n")
            # Get absolute path of the special VASP binary
            f.write("VASP_EXE=$BASE_DIR/vasp_std_str_c\n\n")
            f.write("echo \"Starting packed VASP job (Serial execution)...\"\n\n")
            
            # Write tasks
            for folder in chunk:
                # Ensure we use Linux-style paths if running on Windows
                folder_linux = folder.replace(os.sep, '/')
                
                f.write(f"# --- Task: {folder_linux} ---\n")
                f.write(f"echo \"Processing {folder_linux}...\"\n")
                f.write(f"cd $BASE_DIR/{folder_linux}\n")
                
                # Read the actual sub.sh to execute its commands
                # This ensures we are always in sync with build.py
                sub_sh_path = os.path.join(folder, "sub.sh")
                if os.path.exists(sub_sh_path):
                    with open(sub_sh_path, 'r') as subf:
                        lines = subf.readlines()
                        for line in lines:
                            # Skip scheduler directives and directory changes that break batch mode
                            if line.strip().startswith("#BSUB"):
                                continue
                            if "cd $LS_SUBCWD" in line:
                                continue
                            if line.strip().startswith("#!"): # Skip shebang
                                continue
                            
                            # Replace relative execution with absolute path
                            # This fixes the "execvp error" if symlinks fail
                            if "./vasp_std_str_c" in line:
                                line = line.replace("./vasp_std_str_c", "$VASP_EXE")
                                
                            f.write(line)
                else:
                    f.write(f"echo \"Error: sub.sh not found in {folder_linux}\"\n")
                
                f.write(f"\necho \"Finished {folder_linux}\"\n")
                f.write("sleep 1\n\n") # Brief pause
            
            f.write("echo \"All tasks in this packed job finished.\"\n")

        # Submit the script
        print(f"Generated {script_name}. Submitting...")
        if os.system("which bsub > /dev/null 2>&1") == 0:
            os.system(f"bsub < {script_name}")
        else:
            print(f"  'bsub' command not found. Please submit {script_name} manually.")

if __name__ == "__main__":
    submit_vasp_jobs()
