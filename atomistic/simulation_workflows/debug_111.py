
import numpy as np
from ase.build import bulk, make_supercell
from ase.io import write

def is_orthorhombic(atoms, tol_deg=1e-3):
    ang = np.array(atoms.get_cell().angles(), dtype=float)
    print(f"Angles: {ang}")
    return np.max(np.abs(ang - 90.0)) < tol_deg

def build_oriented_orth_cell(symbol, a0, x, y, z, tol_deg=1e-3):
    base = bulk(symbol, "fcc", a=a0, cubic=True)
    x = np.array(x, dtype=int)
    y = np.array(y, dtype=int)
    z = np.array(z, dtype=int)
    
    print(f"Base cell:\n{base.cell}")
    
    if (np.dot(x, y) != 0) or (np.dot(x, z) != 0) or (np.dot(y, z) != 0):
        raise ValueError(f"Directions not orthogonal:\n x={x}, y={y}, z={z}")
    
    P = np.column_stack([x, y, z]).astype(int)
    print(f"Transformation Matrix P:\n{P}")
    
    cand1 = make_supercell(base, P)
    print("Candidate 1 cell:")
    print(cand1.cell)
    print(f"Lengths: {cand1.cell.lengths()}")
    print(f"Angles: {cand1.cell.angles()}")
    
    if is_orthorhombic(cand1, tol_deg=tol_deg):
        return cand1
        
    cand2 = make_supercell(base, P.T) # This logic in original code seems weird if P columns are new vectors
    print("Candidate 2 cell:")
    print(cand2.cell)
    
    if is_orthorhombic(cand2, tol_deg=tol_deg):
        return cand2
        
    print("Failed to find orthorhombic cell")
    return cand1 # Return anyway for inspection

A0 = 3.615
# x=[1-10], y=[11-2], z=[111]
atoms = build_oriented_orth_cell('Cu', A0, [1,-1,0], [1,1,-2], [1,1,1])
write("debug_111.vasp", atoms, format="vasp", vasp5=True, direct=True)
