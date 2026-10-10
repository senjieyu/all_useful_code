import sys
import matplotlib.pyplot as plt
import numpy as np
import os

def write_lammps_data(filename, atoms, box_bounds, step):
    """
    Write a LAMMPS data file.
    atoms: list of [id, type, x, y, z, ...]
    box_bounds: [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
    """
    try:
        with open(filename, 'w') as f:
            f.write(f"LAMMPS data file via analyze.py at step {step}\n\n")
            f.write(f"{len(atoms)} atoms\n")
            f.write("1 atom types\n\n")
            
            f.write(f"{box_bounds[0][0]} {box_bounds[0][1]} xlo xhi\n")
            f.write(f"{box_bounds[1][0]} {box_bounds[1][1]} ylo yhi\n")
            f.write(f"{box_bounds[2][0]} {box_bounds[2][1]} zlo zhi\n\n")
            
            f.write("Atoms\n\n")
            for atom in atoms:
                # atom is [id, type, x, y, z, ...]
                # data file needs: id type x y z
                f.write(f"{int(atom[0])} {int(atom[1])} {atom[2]} {atom[3]} {atom[4]}\n")
    except Exception as e:
        print(f"Error writing data file {filename}: {e}")

def process_dump(dump_file, output_dir, file_prefix, yield_frame, limit_after_yield=50):
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError as e:
            pass # might exist now

    print(f"Extracting frames from {dump_file}...")
    print(f"Yield frame: {yield_frame}. Saving all before and {limit_after_yield} after.")
        
    try:
        with open(dump_file, 'r') as f:
            frame_idx = 0
            while True:
                line = f.readline()
                if not line:
                    break
                
                if "ITEM: TIMESTEP" in line:
                    try:
                        timestep_line = f.readline()
                        if not timestep_line: break
                        timestep = int(timestep_line.strip())
                        
                        f.readline() # ITEM: NUMBER OF ATOMS
                        natoms_line = f.readline()
                        if not natoms_line: break
                        natoms = int(natoms_line.strip())
                        
                        f.readline() # ITEM: BOX BOUNDS ...
                        
                        box_lines = [f.readline().split() for _ in range(3)]
                        box_bounds = [[float(v) for v in l[:2]] for l in box_lines]
                        
                        f.readline() # ITEM: ATOMS ...
                        
                        # Read atoms
                        atoms = []
                        for _ in range(natoms):
                            parts = f.readline().split()
                            # id type x y z ...
                            atom_data = [float(x) for x in parts] 
                            atoms.append(atom_data)
                        
                        # Decision to save
                        save_frame = False
                        if frame_idx < yield_frame:
                            save_frame = True
                        elif frame_idx < yield_frame + limit_after_yield:
                            save_frame = True
                        
                        if save_frame:
                            # Use zfill for correct sorting
                            out_name = os.path.join(output_dir, f"{file_prefix}_frame_{str(frame_idx).zfill(4)}.data")
                            write_lammps_data(out_name, atoms, box_bounds, timestep)
                        
                        frame_idx += 1
                        
                        if frame_idx >= yield_frame + limit_after_yield:
                            print(f"Reached limit of {limit_after_yield} frames after yield. Stopping extraction.")
                            break
                    except ValueError as e:
                        print(f"Error parsing frame {frame_idx}: {e}")
                        break
                        
    except Exception as e:
        print(f"Error processing dump file {dump_file}: {e}")

def analyze(data_file, plot_file, result_file, direction, dump_file=None, output_dir=None, job_name=None):
    # Load data
    try:
        data = np.loadtxt(data_file, comments='#')
    except Exception as e:
        print(f"Error reading {data_file}: {e}")
        return

    if data.ndim == 1:
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
    
    # Process dump if arguments are provided
    if dump_file and output_dir and job_name:
        process_dump(dump_file, output_dir, f"{job_name}_{direction}", frame, limit_after_yield=50)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python analyze.py <data_file> <plot_file> <result_file> <direction> [dump_file] [output_dir] [job_name]")
        sys.exit(1)
    
    dump_file = None
    output_dir = None
    job_name = None
    
    if len(sys.argv) >= 8:
        dump_file = sys.argv[5]
        output_dir = sys.argv[6]
        job_name = sys.argv[7]
        
    analyze(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], dump_file, output_dir, job_name)
