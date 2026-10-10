from ase.build import fcc111

atoms = fcc111('Cu', size=(1, 2, 24), a=3.615, vacuum=0, orthogonal=True)
print("ASE 24 layers thickness in Z:", atoms.cell[2,2])
