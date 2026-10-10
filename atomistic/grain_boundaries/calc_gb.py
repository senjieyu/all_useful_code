import os
import argparse
import numpy as np

def get_vasp_info(vasp_dir):
    poscar_path = os.path.join(vasp_dir, 'POSCAR')
    if not os.path.exists(poscar_path):
        return None, None, None
    
    with open(poscar_path, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        
    scale = float(lines[1])
    a = np.array([float(x) for x in lines[2].split()]) * scale
    b = np.array([float(x) for x in lines[3].split()]) * scale
    area = np.linalg.norm(np.cross(a, b))
    
    # Num atoms
    try:
        atom_counts = [int(x) for x in lines[5].split()]
    except ValueError:
        atom_counts = [int(x) for x in lines[6].split()]
    n_atoms = sum(atom_counts)
    
    energy = None
    # Parse OUTCAR to get the final energy without entropy
    outcar_path = os.path.join(vasp_dir, 'OUTCAR')
    if os.path.exists(outcar_path):
        with open(outcar_path, 'r') as f:
            for line in f:
                if 'energy  without entropy=' in line:
                    energy = float(line.split()[-1])
    
    # If OUTCAR not found, check vasp_scf.out
    if energy is None:
        out_path = os.path.join(vasp_dir, 'vasp_scf.out')
        if os.path.exists(out_path):
            with open(out_path, 'r') as f:
                for line in f:
                    if 'free  energy   TOTEN' in line:
                        parts = line.split('=')
                        energy_str = parts[1].replace('eV', '').strip()
                        energy = float(energy_str)
                        
    return energy, n_atoms, area

def get_lammps_info(lammps_dir):
    out_path = os.path.join(lammps_dir, 'lammps.out')
    if not os.path.exists(out_path):
        return None, None, None
    
    energy, n_atoms, area = None, None, None
    with open(out_path, 'r') as f:
        for line in f:
            if line.startswith('Final_Energy ='):
                energy = float(line.split('=')[1].strip())
            elif line.startswith('Num_Atoms ='):
                n_atoms = int(line.split('=')[1].strip())
            elif line.startswith('Area ='):
                area = float(line.split('=')[1].strip())
                
    return energy, n_atoms, area

def main():
    parser = argparse.ArgumentParser(description='Calculate Grain Boundary Energy')
    parser.add_argument('--e_eam', type=float, required=True, help='Bulk energy per atom for EAM')
    parser.add_argument('--e_sus2', type=float, required=True, help='Bulk energy per atom for SUS2')
    parser.add_argument('--e_vasp', type=float, required=True, help='Bulk energy per atom for VASP')
    args = parser.parse_args()

    methods = {
        'LAMMPS_EAM': args.e_eam,
        'LAMMPS_SUS2': args.e_sus2,
        'VASP': args.e_vasp
    }

    print(f"{'Method':<15} | {'Structure':<10} | {'E_GB(eV)':<12} | {'N_Atoms':<8} | {'Area(A^2)':<10} | {'Gamma_GB(J/m^2)':<15}")
    print("-" * 80)

    eV_A2_to_J_m2 = 16.0217662

    for method, e_bulk in methods.items():
        for i in range(10):
            calc_dir = os.path.join(method, str(i))
            if method == 'VASP':
                e_gb, n, area = get_vasp_info(calc_dir)
            else:
                e_gb, n, area = get_lammps_info(calc_dir)
                
            if e_gb is not None and n is not None and area is not None:
                gamma_eV = (e_gb - n * e_bulk) / (2 * area)
                gamma_J = gamma_eV * eV_A2_to_J_m2
                print(f"{method:<15} | {i:<10} | {e_gb:<12.4f} | {n:<8} | {area:<10.4f} | {gamma_J:<15.4f}")
            else:
                print(f"{method:<15} | {i:<10} | {'N/A':<12} | {'N/A':<8} | {'N/A':<10} | {'N/A':<15}")
        print("-" * 80)

if __name__ == "__main__":
    main()
