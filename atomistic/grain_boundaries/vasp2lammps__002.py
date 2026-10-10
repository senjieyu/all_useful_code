import sys
from ase.io import read, write

try:
    input_vasp = sys.argv[1]
    output_data = sys.argv[2]
    atoms = read(input_vasp, format='vasp')
    write(output_data, atoms, format='lammps-data')
    print("Success")
except Exception as e:
    print(f"Error: {e}")
