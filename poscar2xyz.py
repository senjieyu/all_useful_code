from ase.io import read,iread,write
import os
import sys
import glob
from tqdm import tqdm
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator


def find_files_by_keyword(root_dir, keyword):
    """
    递归查找文件名中包含指定关键词的所有文件（忽略后缀）

    Args:
        root_dir (str): 要搜索的根目录路径
        keyword (str): 要匹配的关键词（大小写敏感）

    Returns:
        list: 匹配的文件完整路径列表
    """
    matched_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if keyword in filename:  # 检查文件名是否包含关键词
                full_path = os.path.join(dirpath, filename)
                matched_files.append(full_path)
    return matched_files


input_dir = sys.argv[1]

#dirs = glob.glob(os.path.join(input_dir,'*.vasp'))
#dirs = glob.glob(os.path.join(input_dir,'*.cif'))
root_dir = input_dir
keyword= 'POSCAR'
dirs = find_files_by_keyword(root_dir, keyword)

atom_list = []
for index,a in tqdm(enumerate(dirs)):
    atom = read(a)
    label = index
    label_name = os.path.basename(a)
    atom.info['label'] = label
    atom.info['label_name'] = label_name
    # 设置总能量为 0
    total_energy = 0.0

    # 设置所有原子的力为 0
    forces = np.zeros((len(atom), 3))  # 创建一个形状为 (N, 3) 的零矩阵，N 是原子数量

    # 设置应力为 0
    stress = np.zeros(6)  # 压力张量通常是一个 6 维向量

    # 使用 SinglePointCalculator 设置总能量、原子力和应力
    atom.calc = SinglePointCalculator(atom, energy=total_energy, forces=forces, stress=stress)

    # atom.info['energy'] = 0.0
    # forces = np.zeros((len(atom), 3))  # 创建一个形状为 (N, 3) 的零矩阵，N 是原子数量
    # atom.set_forces(forces)
    #
    # # 设置应力为 0
    virial = np.zeros(9)  # 压力张量通常是一个 6 维向量
    atom.info['virial'] = virial
    atom_list.append(atom)
write('out.xyz',atom_list,format='extxyz')


