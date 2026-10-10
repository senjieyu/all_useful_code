import os
import shutil
import numpy as np
from ase.build import fcc111

def main():
    if os.path.exists('stru'):
        shutil.rmtree('stru')
    os.makedirs('stru')

    a = 3.615
    # Build a 1x2 surface cell with 24 layers (48 atoms)
    # 12 layers for bottom region, 12 layers for top region
    atoms = fcc111('Cu', size=(1, 2, 24), a=a, vacuum=0, orthogonal=True)
    
    # ASE standard fcc111 orthogonal:
    # X = [1 -1 0], Y = [1 1 -2], Z = [1 1 1]
    # We rotate it -90 degrees around Z to match LAMMPS:
    # X becomes [1 1 -2] (which is the <112> direction), Y becomes [-1 1 0].
    atoms.rotate(-90, 'z', center=(0,0,0))
    
    # Re-align cell matrix to be strictly diagonal
    cell = atoms.get_cell()
    lx = cell[1,1]
    ly = cell[0,0]
    
    # The Z spacing between layers is a / sqrt(3) ~ 2.087 A
    # Vacuum added to avoid interaction between periodic images in Z
    # One stacking fault plane (internal) -> SFE = ΔE / A
    lz_slab = 24 * (a / np.sqrt(3))
    vacuum = 15.0
    lz_box = lz_slab + vacuum
    
    # Center atoms in the new positive orthant
    pos = atoms.positions
    pos[:, 0] -= np.min(pos[:, 0])
    pos[:, 1] -= np.min(pos[:, 1])
    pos[:, 2] -= np.min(pos[:, 2])
    
    atoms.set_cell([lx, ly, lz_box])
    atoms.set_positions(pos)
    atoms.wrap()
    
    # Find the midpoint to separate region 1 and 2
    # Layers 0-11 (Z < 25) are bottom, Layers 12-23 (Z > 25) are top
    mid_z = 11.5 * (a / np.sqrt(3))
    
    displacements = np.arange(0.0, 2.1, 0.1)
    
    for x_disp in displacements:
        atoms_disp = atoms.copy()
        pos = atoms_disp.positions
        
        # Move the upper half (z > mid_z) by x_disp along X
        pos[pos[:, 2] > mid_z, 0] -= x_disp
        
        # Wrap only XY to keep atoms inside the box
        pos[:, 0] %= lx
        pos[:, 1] %= ly
        
        atoms_disp.set_positions(pos)
        
        filename = f"stru/Cu_GSFE_disp{x_disp:.1f}.vasp"
        with open(filename, 'w') as f:
            f.write(f"Cu GSFE disp {x_disp:.1f} (24 layers)\n")
            f.write("1.0\n")
            f.write(f"  {lx:.9f} 0.000000000 0.000000000\n")
            f.write(f"  0.000000000 {ly:.9f} 0.000000000\n")
            f.write(f"  0.000000000 0.000000000 {lz_box:.9f}\n")
            f.write("  Cu\n")
            f.write(f"  {len(atoms)}\n")
            f.write("Selective dynamics\n")
            f.write("Cartesian\n")
            for i in range(len(atoms)):
                # Fix XY, relax Z (F F T)
                f.write(f"  {pos[i,0]:.9f} {pos[i,1]:.9f} {pos[i,2]:.9f}  F F T\n")
                
    print(f"Generated {len(displacements)} structures in stru/ with 24 layers (48 atoms, no vacuum, shear model).")

if __name__ == '__main__':
    main()
