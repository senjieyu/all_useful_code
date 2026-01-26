import os
import sys
import shutil
import glob

def main():
    if len(sys.argv) < 3:
        print("Usage: python collect_elastic.py <work_dir> <final_out_dir>")
        sys.exit(1)
        
    work_dir = sys.argv[1]
    final_out_dir = sys.argv[2]
    
    if not os.path.exists(work_dir):
        print(f"Work directory {work_dir} not found.")
        sys.exit(1)
        
    os.makedirs(final_out_dir, exist_ok=True)
    
    struct_dirs = sorted(glob.glob(os.path.join(work_dir, "struct_*")))
    
    count = 0
    for s_dir in struct_dirs:
        s_name = os.path.basename(s_dir) # struct_0
        
        # Result file
        res_file = os.path.join(s_dir, "elastic_properties.txt")
        
        if os.path.exists(res_file):
            dest_name = f"{s_name}.txt"
            shutil.copy(res_file, os.path.join(final_out_dir, dest_name))
            count += 1
        else:
            print(f"Warning: No result found in {s_name}")
            
    print(f"Collected {count} result files to {final_out_dir}")

if __name__ == "__main__":
    main()
