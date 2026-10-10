# batch_convert.py
import os
import shutil
from pathlib import Path
from ase.io import read, write
from ase.data import atomic_numbers, chemical_symbols


def convert_vasp_to_lammps(vasp_file, ele_list, size, ele_model=1, output_dir="lammps_data"):
    """单个文件转换函数"""
    atoms = read(vasp_file)
    atoms = atoms.repeat(size)

    if ele_model == 1:
        ele = [atomic_numbers[i] for i in ele_list]
        ele = sorted(ele)
        specorder = [chemical_symbols[i] for i in ele]
    elif ele_model == 2:
        specorder = ele_list
    else:
        raise ValueError('ele_model must be 1 or 2!')

    # 新版本：直接保存到输出文件夹根目录
    filename = Path(vasp_file).stem
    output_path = Path(output_dir)  # ← 直接使用输出文件夹
    output_path.mkdir(parents=True, exist_ok=True)

    output_file = output_path / f"{filename}.data"  # ← 直接保存到输出文件夹
    write(output_file, atoms, format='lammps-data',
          masses=True, specorder=specorder, force_skew=True)

    return str(output_file)


def batch_convert_folder(vasp_folder="poscar_folder", output_folder="lammps_data",
                         ele_list=None, size=(2, 2, 2), ele_model=1):
    """批量转换文件夹中的所有.vasp文件"""
    if ele_list is None:
        ele_list = ['H', 'C', 'N', 'O', 'Cl', 'K', 'Fe']

    vasp_dir = Path(vasp_folder)
    output_dir = Path(output_folder)
    output_dir.mkdir(exist_ok=True)

    vasp_files = list(vasp_dir.glob("*.vasp"))

    if not vasp_files:
        print(f"在 {vasp_folder} 中未找到.vasp文件")
        return

    print(f"找到 {len(vasp_files)} 个.vasp文件，开始转换...")

    success_count = 0
    for vasp_file in vasp_files:
        try:
            output_file = convert_vasp_to_lammps(
                vasp_file, ele_list, size, ele_model, output_folder
            )
            print(f"  ✓ {vasp_file.name} -> {Path(output_file).relative_to(output_dir)}")
            success_count += 1
        except Exception as e:
            print(f"  ✗ 转换 {vasp_file.name} 失败: {str(e)}")

    print(f"\n转换完成！成功: {success_count}/{len(vasp_files)}")
    print(f"输出目录: {output_dir.absolute()}")


# 使用示例
if __name__ == "__main__":
    # 配置参数
    ELEMENTS = ['Al', 'Ti', 'V', 'Nb', 'Hf']  # 修改为实际元素
    SIZE = (1, 1, 1)  # 扩胞倍数
    ELE_MODEL = 1

    # 运行批量转换
    batch_convert_folder(
        vasp_folder="poscar",  # 输入文件夹
        output_folder="lammps_data",  # 输出文件夹
        ele_list=ELEMENTS,
        size=SIZE,
        ele_model=ELE_MODEL
    )