import sys
import matplotlib.pyplot as plt
import numpy as np
import os

def analyze(data_file, plot_file, result_file, direction):
    # Load data
    try:
        # Skip the first line if it's a header (LAMMPS fix print usually has a header line)
        # But if screen no is used, it might just be the data.
        # fix print ... file ... adds a header line starting with # usually
        data = np.loadtxt(data_file, comments='#')
    except Exception as e:
        print(f"Error reading {data_file}: {e}")
        return

    if data.ndim == 1:
        # Only one data point
        data = data.reshape(1, -1)

    if data.shape[0] == 0:
        print("Empty data file")
        return

    strain = data[:, 0]
    # Select stress column: strain sigma_xx sigma_yy sigma_zz
    if direction == 'x':
        stress = data[:, 1]
    elif direction == 'y':
        stress = data[:, 2]
    elif direction == 'z':
        stress = data[:, 3]
    else:
        stress = data[:, 1]

    # Plot
    plt.figure()
    plt.plot(strain, stress, label=f'Tensile {direction.upper()}')
    plt.xlabel('Strain')
    plt.ylabel('Stress (GPa)')
    plt.title(f'Stress-Strain Curve ({direction.upper()})')
    plt.legend()
    plt.grid(True)
    plt.savefig(plot_file)
    plt.close()

    # Find yield point (Max Stress)
    max_stress_idx = np.argmax(stress)
    max_stress = stress[max_stress_idx]
    yield_strain = strain[max_stress_idx]
    
    # Constants from LAMMPS input
    print_freq = 200
    dump_freq = 250
    
    # Estimate timestep
    # fix print output usually starts at step 0 or step print_freq depending on run
    # Let's assume line 0 -> step 0 (if printed) or step 200.
    # Usually LAMMPS prints step 0. 
    current_step = max_stress_idx * print_freq 
    
    # Frame number for dump
    frame = int(round(current_step / dump_freq))

    # Append to result file
    with open(result_file, 'a') as f:
        f.write(f"Direction: {direction}\n")
        f.write(f"Max Stress: {max_stress:.4f} GPa\n")
        f.write(f"Yield Strain: {yield_strain:.4f}\n")
        f.write(f"Approximate Timestep: {current_step}\n")
        f.write(f"Yield Frame: {frame}\n")
        f.write("-" * 20 + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python analyze.py <data_file> <plot_file> <result_file> <direction>")
        sys.exit(1)
        
    analyze(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
