from ase.build import fcc111
atoms = fcc111('Cu', size=(1, 2, 6), a=3.615, vacuum=10, orthogonal=True)
print(atoms.cell)
print(len(atoms))
