import os
import sys
import shutil
import subprocess
import glob

def main():
    if len(sys.argv) < 3:
        print("Usage: python setup_batch.py <input.xyz> <output_folder>")
        print("Example: python setup_batch.py CU331.xyz batch_work")
        sys.exit(1)

    input_xyz = sys.argv[1]
    output_folder = sys.argv[2]
    
    # Get absolute paths
    cwd = os.getcwd()
    input_xyz = os.path.abspath(input_xyz)
    output_folder = os.path.abspath(output_folder)
    
    # 1. Run rand_dis_battery_data.py to generate folders and data
    print("--- Step 1: Generating EOS data folders ---")
    script_path = os.path.join(cwd, "rand_dis_battery_data.py")
    
    if not os.path.exists(script_path):
        print(f"Error: {script_path} not found.")
        sys.exit(1)

    # Call the script
    cmd = [sys.executable, script_path, input_xyz, output_folder]
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print("Error generating EOS data.")
        sys.exit(1)
    
    # 2. Prepare files to copy
    # User specified: mlip.ini, relax_isif.in (assuming relax_isif2.in based on file list), current.mtp
    files_to_copy = ["mlip.ini", "relax_isif2.in", "current.mtp"]
    submit_template = "submit_eos_isif2.sh"
    
    # Check existence
    for f in files_to_copy + [submit_template]:
        if not os.path.exists(f):
            print(f"Error: Required file {f} not found in current directory.")
            sys.exit(1)
            
    # Read submit template
    with open(submit_template, 'r', encoding='utf-8') as f:
        submit_content = f.read()
        
    # Create submission_scripts folder (Folder ②)
    scripts_dir = os.path.join(output_folder, "submission_scripts")
    os.makedirs(scripts_dir, exist_ok=True)
    
    # 3. Iterate over generated struct folders
    # rand_dis_battery_data.py generates folders named "struct_X"
    struct_dirs = sorted(glob.glob(os.path.join(output_folder, "struct_*")))
    
    if not struct_dirs:
        print("No struct directories found! Check if input xyz had valid frames.")
        sys.exit(1)
        
    print(f"--- Step 2: Setting up {len(struct_dirs)} job folders ---")
    
    for s_dir in struct_dirs:
        folder_name = os.path.basename(s_dir) # e.g. struct_0
        
        print(f"Setting up {folder_name}...")
        
        # Copy files to the structure folder
        for f in files_to_copy:
            shutil.copy(f, os.path.join(s_dir, f))
            
        # Customize and write submit script
        # Change job name to include folder name
        job_name = f"eos_{folder_name}"
        new_content = submit_content.replace("#BSUB -J eos_mlip_isif2", f"#BSUB -J {job_name}")
        
        # Write submit.sh inside the task folder
        submit_path = os.path.join(s_dir, "submit.sh")
        with open(submit_path, 'w', encoding='utf-8') as f:
            f.write(new_content)
            
        # Create wrapper script in submission_scripts folder
        # This script cds to the task folder and submits the job
        wrapper_path = os.path.join(scripts_dir, f"submit_{folder_name}.sh")
        
        # Use absolute path for s_dir to ensure it works from anywhere
        abs_s_dir = os.path.abspath(s_dir)
        # On Windows, path separators might be backslashes, but bash needs forward slashes or escaped.
        # Assuming this runs in a bash-like environment (WSL or Linux cluster), we should ensure paths are POSIX-like if possible.
        # However, Python on Windows produces backslashes. 
        # If the user is submitting from Windows to a Linux cluster via some mechanism, paths might be tricky.
        # But usually 'bsub' implies a Linux environment.
        # I'll try to convert to forward slashes just in case.
        abs_s_dir_posix = abs_s_dir.replace(os.sep, '/')
        
        with open(wrapper_path, 'w', encoding='utf-8', newline='\n') as f:
            f.write("#!/bin/bash\n")
            f.write(f"cd \"{abs_s_dir_posix}\"\n")
            f.write("bsub < submit.sh\n")
            f.write("echo \"Submitted job for %s\"\n" % folder_name)
            
    print("--- Setup Complete ---")
    print(f"Job folders prepared in: {output_folder}")
    print(f"Submission scripts generated in: {scripts_dir}")
    print("You can now run 'python submit_batch.py <submission_scripts_folder>' to submit all jobs.")

if __name__ == "__main__":
    main()
