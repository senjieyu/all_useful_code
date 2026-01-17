import os.path
import sys

from ase.io import read,write
#from ase.geometry import orthorhombicize
import glob
import numpy as np
from ase import Atoms
import spglib
import ast
from pymatgen.io.ase import AseAtomsAdaptor

def find_primitive_cell(atoms: Atoms) -> Atoms:
    """
    使用 spglib 找到原胞。
    """
    # 获取 ASE Atoms 对象的晶格和原子信息
    lattice = atoms.get_cell().array
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()

    # 使用 spglib 标准化晶格并找到原胞
    cell = (lattice, positions, numbers)
    primitive_cell = spglib.standardize_cell(cell, to_primitive=True, no_idealize=False)

    # 将 spglib 的输出转换回 ASE Atoms 对象
    primitive_atoms = Atoms(
        numbers=primitive_cell[2],
        scaled_positions=primitive_cell[1],
        cell=primitive_cell[0],
        pbc=atoms.get_pbc()
    )
    return primitive_atoms

def super_atoms(input_files,output_dir):
    atoms = read(input_files)  # 替换为你的文件路径
    num = atoms.get_number_of_atoms()
    print(f'num:{num} cell:{atoms.cell.cellpar()}')
    size = input("Please enter size: ")

    atoms = atoms.repeat(ast.literal_eval(size))

    # if num < 95:
    #     a, b, c = atoms.cell.lengths()
    #     for index, t in enumerate([a,b,c]):
    #         if t<9:
    #             size[index] = 2
    #     # if size == [1,1,1]:
    #     #     print(os.path.basename(input_files),set(atoms.get_chemical_symbols()),atoms.get_chemical_formula())
    #     atoms = atoms.repeat(tuple(size))
    #     #print(atoms.get_cell_lengths_and_angles(),size,num,num_2)
    # else:
    #     #print(atoms.get_cell_lengths_and_angles(),num)
    #     pass

    b = atoms
    sorted_indices = sorted(range(len(b)), key=lambda i: b[i].symbol)
    sorted_atoms = b.__class__()
    sorted_atoms.cell = b.cell
    sorted_atoms.pbc = b.pbc
    for i in sorted_indices:
        sorted_atoms.append(b[i])
    write(os.path.join(output_dir,os.path.basename(input_files)), sorted_atoms, format='vasp')

    # else:
    #     kk = find_primitive_cell(atoms)
    #     a, b, c = kk.cell.lengths()
    #     print(a, b, c)

    # tt = atoms.get_cell_lengths_and_angles()
    # formula = atoms.get_chemical_formula()
    #
    # kk = find_primitive_cell(atoms)
    # print(formula,tt,num,kk.get_number_of_atoms(),kk.get_cell_lengths_and_angles())
    # # 获取正交后的晶格常数
    # a, b, c = atoms.cell.lengths()
    # print(f'正交后的晶格常数 a = {a}, b = {b}, c = {c}')
from pymatgen.core import Structure

def super_atoms_2(file,outpath):
    ori_structure = Structure.from_file(file)
    structure = Structure.from_file(file)

    # 定义目标晶格常数
    target_a = 9
    target_b = 9
    target_c = 9

    # 计算扩胞矩阵
    a, b, c = structure.lattice.abc
    scaling_factors = [max(1, int(np.ceil(target_a / a))),
                       max(1, int(np.ceil(target_b / b))),
                       max(1, int(np.ceil(target_c / c)))]
    supercell_matrix = np.diag(scaling_factors)

    # 扩胞
    supercell = structure.make_supercell(supercell_matrix)

    if max(ori_structure.lattice.abc)<100:
        # 打印扩胞后的结构信息
        print('----------------------------------------------------')
        print(f"Lattice: {ori_structure.lattice.abc}", )
        print(f"angles: {ori_structure.lattice.angles}")
        print("Supercell lattice parameters:", supercell.lattice.abc)
        print("Supercell matrix:", supercell_matrix)
        print("Number of atoms in supercell:", len(supercell.sites))
        out_file = os.path.join(outpath,os.path.basename(file))
        ase_atoms = AseAtomsAdaptor.get_atoms(supercell)
        write(out_file,ase_atoms,format='vasp')
        return True
    else:
        return False

if __name__ == '__main__':
    # files = glob.glob(os.path.join(r'B_vasp','*.vasp'))
    # out_dir = os.path.join(os.getcwd(),'super_B')
    # if not os.path.exists(out_dir):
    #     os.mkdir(out_dir)
    # for a in files:
    #     super_atoms(a,out_dir)
    #
    #
    # ff = glob.glob(os.path.join(out_dir,'*.vasp'))
    #
    # for a in ff:
    #     print(os.path.basename(a),read(a).get_cell_lengths_and_angles(),read(a).get_number_of_atoms())
        #print(os.path.basename(a),read(a).get_chemical_formula(),set(read(a).get_chemical_symbols()))


    input_dir = sys.argv[1]
    input_dir_files = glob.glob(os.path.join(input_dir,'*.vasp'))
    out_dir = sys.argv[2]

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    count = 0
    for a in input_dir_files:
        #if len(list(set(read(a).get_chemical_symbols()))) == 3:
        if super_atoms_2(a,out_dir):
            print(list(set(read(a).get_chemical_symbols())))
            count +=1
    print(count)
