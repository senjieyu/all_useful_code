import os
import sys
import glob

def main():
    # Default batch folder name if not provided
    default_batch_dir = "batch_work"
    
    if len(sys.argv) > 1:
        batch_dir = sys.argv[1]
    else:
        batch_dir = default_batch_dir
        
    if not os.path.exists(batch_dir):
        print(f"Error: Directory '{batch_dir}' not found.")
        print("Usage: python collect.py [batch_work_folder]")
        sys.exit(1)

    # Create output directory
    output_dir = "out_xyz"
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory created: {os.path.abspath(output_dir)}")

    # Find all struct_* directories
    struct_dirs = sorted(glob.glob(os.path.join(batch_dir, "struct_*")))
    
    if not struct_dirs:
        print(f"No 'struct_*' directories found in '{batch_dir}'.")
        sys.exit(0)

    print(f"Found {len(struct_dirs)} structure directories.")

    count = 0
    for s_dir in struct_dirs:
        struct_name = os.path.basename(s_dir) # e.g., "struct_0"
        
        # Path to eos_xyz inside the struct folder
        eos_xyz_dir = os.path.join(s_dir, "eos_xyz")
        
        if not os.path.exists(eos_xyz_dir):
            print(f"Warning: 'eos_xyz' folder not found in {struct_name}, skipping.")
            continue
            
        # Find all .extxyz files (based on submit script) or .xyz files
        # We look for both just in case
        xyz_files = glob.glob(os.path.join(eos_xyz_dir, "*.extxyz")) + \
                    glob.glob(os.path.join(eos_xyz_dir, "*.xyz"))
        
        # Sort files to maintain order (e.g. by volume if named correctly)
        xyz_files.sort()
        
        if not xyz_files:
            print(f"Warning: No xyz files found in {struct_name}/eos_xyz, skipping.")
            continue

        # Target merged file path
        merged_filename = f"{struct_name}.xyz"
        merged_path = os.path.join(output_dir, merged_filename)
        
        # Merge content
        with open(merged_path, 'w', encoding='utf-8') as outfile:
            for f in xyz_files:
                try:
                    with open(f, 'r', encoding='utf-8') as infile:
                        outfile.write(infile.read())
                except Exception as e:
                    print(f"Error reading {f}: {e}")
        
        print(f"Generated: {merged_path} ({len(xyz_files)} frames)")
        count += 1

    print("-" * 30)
    print(f"Collection complete. {count} files merged into '{output_dir}'.")

if __name__ == "__main__":
    main()
