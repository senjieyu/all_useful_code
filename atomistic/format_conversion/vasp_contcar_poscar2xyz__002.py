from ase.io import read, write
import os, sys
import numpy as np
from tqdm import tqdm
from ase.calculators.singlepoint import SinglePointCalculator

def find_poscar_like_files(root_dir):
    """递归找 POSCAR/CONTCAR/vasp 文件"""
    patterns = ("POSCAR", "CONTCAR")
    matched = []

    for dirpath, _, filenames in os.walk(root_dir):
        for fn in filenames:
            full = os.path.join(dirpath, fn)

            # 1) 文件名就是 POSCAR / CONTCAR（常见）
            if fn in patterns:
                matched.append(full)
                continue

            # 2) 以 .vasp 结尾（你这种）
            if fn.lower().endswith(".vasp"):
                matched.append(full)
                continue

    return sorted(matched)

def main():
    if len(sys.argv) < 3:
        print("用法: python poscar2xyz.py <输入文件夹> <输出xyz文件名>")
        sys.exit(1)

    input_dir = sys.argv[1]
    out_xyz = sys.argv[2]

    files = find_poscar_like_files(input_dir)
    print(f"找到 {len(files)} 个 VASP 结构文件")

    if len(files) == 0:
        print("❌ 没找到 POSCAR/CONTCAR 或 *.vasp 文件。检查目录和文件名。")
        sys.exit(1)

    atom_list = []
    for idx, fp in tqdm(list(enumerate(files)), total=len(files)):
        atoms = read(fp)  # ASE 自动识别 POSCAR/CONTCAR/vasp

        atoms.info["label"] = idx
        atoms.info["label_name"] = os.path.basename(fp)

        # 设置能量/力/应力为 0（extxyz 需要这些字段时好用）
        forces = np.zeros((len(atoms), 3))
        stress6 = np.zeros(6)
        atoms.calc = SinglePointCalculator(atoms, energy=0.0, forces=forces, stress=stress6)

        atoms.info["virial"] = np.zeros(9)  # 你原来的习惯保留

        atom_list.append(atoms)

    write(out_xyz, atom_list, format="extxyz")
    print(f"✅ 写出: {out_xyz}  (frames={len(atom_list)})")

if __name__ == "__main__":
    main()

