
import numpy as np
from ase.build import fcc111

atoms = fcc111('Cu', size=(1, 2, 24), a=3.615, vacuum=0, orthogonal=True)
cell = atoms.get_cell()
print("Original Cell:")
print(cell)
print("Positions first 5:")
print(atoms.positions[:5])

# Rotate -90 z
atoms.rotate(-90, 'z', center=(0,0,0))
print("\nRotated Positions first 5:")
print(atoms.positions[:5])

# Check what build_structures.py does
lx = cell[1,1]
ly = cell[0,0]
print(f"\nExtracted lx (cell[1,1]): {lx}")
print(f"Extracted ly (cell[0,0]): {ly}")

# Check if rotated atoms fit in [lx, ly] box
# Original X is along 1-10. Length?
# Original Y is along 11-2. Length?

# ASE fcc111 orthogonal=True:
# Vector 1 (X) is along <1 -1 0>
# Vector 2 (Y) is along <1 1 -2>
# Length of X (size 1) = a * sqrt(2) / 2? No.
# a_fcc = 3.615
# Nearest neighbor dist = a / sqrt(2) ~ 2.556
# Vector <1 -1 0> length = a * sqrt(2).
# But orthogonal cell usually uses minimal rectangular.
# Rectangular unit in (111):
# v1 = 1/2 [1, -1, 0] * a? No, that's length a/sqrt(2).
# v2 = 1/2 [1, 1, -2] * a? No.
# Let's see what ASE gives.
