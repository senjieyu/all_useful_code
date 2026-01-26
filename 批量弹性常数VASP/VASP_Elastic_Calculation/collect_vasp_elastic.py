import os
import sys
import shutil
import glob

def parse_vasp_elastic_outcar(outcar_path):
    """
    Parses OUTCAR for elastic moduli.
    Looks for "TOTAL ELASTIC MODULI (kBar)"
    Returns lines containing the matrix and moduli.
    """
    if not os.path.exists(outcar_path):
        return None
        
    elastic_data = []
    capture = False
    
    with open(outcar_path, 'r') as f:
        lines = f.readlines()
        
    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI (kBar)" in line:
            capture = True
            # Capture the next ~20 lines
            elastic_data.append("Found TOTAL ELASTIC MODULI (kBar) in OUTCAR:\n")
            # Skip header lines
            continue
            
        if capture:
            if "----------------------------------------------------------------" in line:
                # End of block usually looks like a separator
                # But VASP OUTCAR format varies.
                # Just capture until we see something else or fixed lines.
                pass
                
            elastic_data.append(line)
            
            # Stop condition (heuristic)
            if len(elastic_data) > 30: 
                capture = False
                
    if not elastic_data:
        return "Elastic moduli not found in OUTCAR (calculation might have failed or not finished)."
        
    return "".join(elastic_data)

def main():
    if len(sys.argv) < 3:
        print("Usage: python collect_vasp_elastic.py <work_dir> <final_out_dir>")
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
        outcar = os.path.join(s_dir, "OUTCAR")
        
        result_text = parse_vasp_elastic_outcar(outcar)
        
        if result_text:
            out_file = os.path.join(final_out_dir, f"{s_name}_elastic.txt")
            with open(out_file, 'w') as f:
                f.write(f"Results for {s_name}\n")
                f.write("-" * 30 + "\n")
                f.write(result_text)
            count += 1
        else:
            print(f"Warning: OUTCAR not found or invalid in {s_name}")
            
    print(f"Collected results for {count} structures to {final_out_dir}")

if __name__ == "__main__":
    main()
