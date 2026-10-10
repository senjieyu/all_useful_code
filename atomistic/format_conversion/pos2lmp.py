import os
from ase.data import atomic_numbers
from ase.io import read,write
from ase.data import atomic_masses
from ase.data import chemical_symbols
import sys

def pos2lmp(ele_,size,ele_model,POSCAR):
    # if os.path.exists('POSCAR'):
    #     POSCAR = 'POSCAR'
    # else:
    #     POSCAR = [f for f in os.listdir() if f.endswith('.vasp')][0]
    a = read(POSCAR)
    a = a.repeat(size)

    if ele_model == 1:
        ele = [atomic_numbers[i] for i in ele_]
        ele = sorted(ele)
        specorder = [chemical_symbols[i] for i in ele]

    elif ele_model == 2:
        specorder = ele_
    else:
        raise ValueError('ele_model must be 1 or 2!')

    write('data.in', a, format='lammps-data',masses=True,specorder=specorder,force_skew=True)

if __name__ == '__main__':
    ele_ = ['H', 'C', 'N', 'O', 'Cl', 'K', 'Fe']
    ele_model = 1
    size = (2,2,2)
    #POSCAR = 'LGPS_stable.vasp' #sys.argv[1]
    POSCAR = sys.argv[1]
    pos2lmp(ele_, size, ele_model,POSCAR)