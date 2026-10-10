import os
import sys
import subprocess
import glob

def main():
    if len(sys.argv) < 2:
        print("Usage: python submit_batch.py <submission_scripts_folder>")
        print("Example: python submit_batch.py batch_work/submission_scripts")
        sys.exit(1)
        
    scripts_dir = sys.argv[1]
    
    if not os.path.isdir(scripts_dir):
        print(f"Error: {scripts_dir} is not a directory")
        sys.exit(1)
        
    # Look for .sh files
    scripts = sorted(glob.glob(os.path.join(scripts_dir, "submit_*.sh")))
    
    if not scripts:
        print("No submit_*.sh scripts found in the folder.")
        sys.exit(1)
        
    print(f"Found {len(scripts)} submission scripts in {scripts_dir}")
    print("Starting batch submission...")
    
    success_count = 0
    fail_count = 0
    
    for script in scripts:
        script_name = os.path.basename(script)
        print(f"Submitting {script_name}...")
        
        # We assume 'bash' is available to run the .sh files
        # If on Windows without bash in path, this might fail.
        # But 'bsub' environment strongly implies Linux/Bash availability.
        cmd = ["bash", script]
        
        try:
            subprocess.check_call(cmd)
            success_count += 1
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit {script_name}. Error: {e}")
            fail_count += 1
        except FileNotFoundError:
            print("Error: 'bash' command not found. Ensure you are in a bash-compatible environment.")
            fail_count += 1
            break

    print("-" * 30)
    print(f"Submission Summary: {success_count} succeeded, {fail_count} failed.")

if __name__ == "__main__":
    main()
