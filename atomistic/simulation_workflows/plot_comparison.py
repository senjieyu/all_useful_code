import os
import matplotlib.pyplot as plt
import numpy as np

# Configuration
base_path = r'c:\Users\senji\Desktop\拉伸测试'
dirs = ['X', 'Y', 'Z']
structs = ['0', '1', '2', '3', '111']
output_dir = os.path.join(base_path, 'comparison_plots')

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def read_data(filepath):
    data = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip().startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    # strain, sxx, syy, szz
                    try:
                        data.append([float(x) for x in parts[:4]])
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return np.array([])
        
    return np.array(data)

for d in dirs:
    res_dir = f'results_{d}'
    # Column index for stress: X->1, Y->2, Z->3 (0 is strain)
    # The file format is: Strain Sxx Syy Szz
    stress_col = 1 if d == 'X' else (2 if d == 'Y' else 3)
    
    for s in structs:
        filename = f'{s}_stress_strain.dat'
        
        path_eam = os.path.join(base_path, 'EAM', res_dir, filename)
        path_sus2 = os.path.join(base_path, 'SUS2', res_dir, filename)
        
        data_eam = read_data(path_eam)
        data_sus2 = read_data(path_sus2)
        
        plt.figure(figsize=(10, 6))
        
        has_data = False
        if data_eam.size > 0:
            plt.plot(data_eam[:, 0], data_eam[:, stress_col], label='EAM', linewidth=1.5, alpha=0.8)
            has_data = True
        else:
            print(f"Warning: No data for EAM {d} {s}")
            
        if data_sus2.size > 0:
            plt.plot(data_sus2[:, 0], data_sus2[:, stress_col], label='SUS2', linewidth=1.5, alpha=0.8)
            has_data = True
        else:
            print(f"Warning: No data for SUS2 {d} {s}")
            
        if has_data:
            plt.title(f'Stress-Strain Comparison: Direction {d}, Structure {s}')
            plt.xlabel('Strain')
            plt.ylabel(f'Stress (Direction {d}) [GPa/Units]')
            plt.legend()
            plt.grid(True, linestyle='--', alpha=0.7)
            
            out_name = f'comparison_{d}_{s}.png'
            out_path = os.path.join(output_dir, out_name)
            plt.savefig(out_path, dpi=300)
            print(f"Saved {out_name}")
        else:
            print(f"Skipping plot for {d} {s} due to missing data")
            
        plt.close()

print("Processing complete.")
