import os
import sys
import subprocess
import platform

def main():
    if len(sys.argv) < 2:
        print("Usage: python submit_batch.py <submission_scripts_dir>")
        sys.exit(1)

    scripts_dir = sys.argv[1]
    if not os.path.isdir(scripts_dir):
        print(f"Error: {scripts_dir} is not a directory.")
        sys.exit(1)

    # Find scripts
    scripts = [f for f in os.listdir(scripts_dir) if f.startswith("submit_struct_") and f.endswith(".sh")]
    # Sort by number in filename if possible
    try:
        scripts.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
    except:
        scripts.sort()

    print(f"Found {len(scripts)} submission scripts in {scripts_dir}")

    is_windows = (os.name == 'nt')

    for script in scripts:
        script_path = os.path.join(scripts_dir, script)
        abs_path = os.path.abspath(script_path)
        
        print(f"Processing {script}...")
        
        if is_windows:
            print(f"  [Windows Warning] Cannot execute bash submission script directly.")
            print(f"  Script content: {open(abs_path, 'r').read().strip()}")
            print(f"  Skipping actual execution.")
        else:
            try:
                # Ensure executable
                st = os.stat(abs_path)
                os.chmod(abs_path, st.st_mode | 0o111)
                
                # Execute
                # We use bash to run the script
                subprocess.check_call(["bash", abs_path])
                print(f"  Submitted.")
            except subprocess.CalledProcessError as e:
                print(f"  Failed to submit {script}: {e}")
            except Exception as e:
                print(f"  Error: {e}")

if __name__ == "__main__":
    main()
