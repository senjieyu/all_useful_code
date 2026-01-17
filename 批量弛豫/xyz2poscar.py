from ase.io import iread,write
from ase import Atoms
import os
import sys

def mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

kk = sys.argv[1]
atoms = list(iread(kk))
dir = 'xyz2poscar'
mkdir(dir)

for index,atom in enumerate(atoms):
    file_name = index #atom.info['label']
    path = os.path.join(dir,str(file_name)+'.vasp')
    write(path, atom, format='vasp')


