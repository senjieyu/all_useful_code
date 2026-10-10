import os
import sys
from tqdm import tqdm
from ase.io import read, write

def find_files(root_dir, filename):
    file_list = []
    for root, dirs, files in os.walk(root_dir):
        if filename in files:
            file_list.append(os.path.join(root, filename))
    return file_list


def collect_poscar(file_list, out_file):
    for file in tqdm(file_list):
        try:
            atoms = read(file)
            atoms.info['path'] = os.path.dirname(file)
            write(out_file, atoms, format='extxyz', append=True)
        except Exception as e:
            print("skip:", file, e)


if __name__ == '__main__':

    root_dir = sys.argv[1]

    # 可以改成 POSCAR 或 CONTCAR
    filename = "POSCAR"

    out_file = "structures.xyz"

    file_list = find_files(root_dir, filename)

    collect_poscar(file_list, out_file)

