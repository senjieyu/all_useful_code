# -*- coding: utf-8 -*-
# Task_Grid_Cutter.py
#
# 【均匀网格切胞任务模块】
# 功能说明：
# 1. 遍历指定目录（通常是 twin_complex）下的所有 VASP 结构文件。
# 2. 将每个结构文件按照设定的网格尺寸 (Grid Size) 进行均匀切割。
#    - 类似于切魔方：将大盒子切分为多个小盒子。
#    - 自动计算每个方向需要切多少份。
# 3. 筛选并保存符合原子数要求的子块。
#
# 主要用途：
# - 将大型复杂结构（如孪晶）分解为适合后续模拟的小型样本。

import os
import numpy as np
import tempfile
from math import ceil
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.atoms import Atoms

# 兼容不同版本 ASE 的 VASP writer
try:
    from ase.io.vasp import write_vasp as _write_vasp
except Exception:
    _write_vasp = None


def ensure_dir(d):
    """确保目录存在，若不存在则创建"""
    os.makedirs(d, exist_ok=True)


def write_vasp_compat(atoms: Atoms, path: str):
    """
    兼容旧版 ASE 的 VASP 文件写入方法
    """
    atoms.wrap() # 确保原子在盒子内

    if _write_vasp is None:
        from ase.io import write
        try:
            write(path, atoms, format="vasp", vasp5=True, sort=True)
        except TypeError:
            write(path, atoms, format="vasp")
        return

    try:
        _write_vasp(path, atoms, vasp5=True, sort=True, direct=True)
        return
    except TypeError:
        pass

    _write_vasp(path, atoms)


def read_vasp_safe(filepath, default_element="Cu"):
    """
    安全读取 VASP 文件，自动修复 'Type_1' 等未知元素符号
    """
    try:
        # 尝试直接读取
        return ase_read(filepath)
    except Exception as e:
        # 如果失败，尝试修复元素行
        # VASP 5 格式中，第 6 行 (index 5) 是元素符号
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            if len(lines) < 7:
                raise ValueError("File too short")
                
            # 检查第 6 行
            line5 = lines[5].strip()
            if "Type_" in line5 or not line5.isalpha(): # 简单启发式
                print(f"    [Fixing] Found invalid element line: '{line5}'. Replacing with '{default_element}'")
                # 假设所有原子都是同一种元素，或者简单替换
                # 如果有多个 Type_1 Type_2，这里只处理简单的单元素情况，或者全部替换为 default
                # 统计原子种类数
                parts = line5.split()
                new_line = " ".join([default_element] * len(parts)) + "\n"
                lines[5] = new_line
                
                with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vasp') as tmp:
                    tmp.writelines(lines)
                    tmp_path = tmp.name
                
                try:
                    atoms = ase_read(tmp_path)
                    return atoms
                finally:
                    try:
                        os.remove(tmp_path)
                    except:
                        pass
            
            # 如果不是元素行的问题，重新抛出异常
            raise e
            
        except Exception as fix_err:
            print(f"    [Read Failed] Could not fix or read file: {fix_err}")
            return None


def cut_single_file(source_file, output_base_dir, config):
    """
    对单个文件进行均匀网格切割
    """
    default_element = config.get("element", "Cu")
    atoms = read_vasp_safe(source_file, default_element)
    
    if atoms is None:
        return 0

    # 获取配置
    grid_size = np.array(config.get("grid_size", (10.0, 10.0, 10.0)))
    atom_limit = config.get("atom_limit", (10, 10000))
    perturb_amp = config.get("perturb_amp", 0.0)
    
    cell = atoms.get_cell()
    lengths = atoms.cell.lengths()
    
    # 计算切分份数 (向上取整)
    nx = max(1, int(ceil(lengths[0] / grid_size[0])))
    ny = max(1, int(ceil(lengths[1] / grid_size[1])))
    nz = max(1, int(ceil(lengths[2] / grid_size[2])))
    
    # 实际切块大小
    dx = cell[0] / nx
    dy = cell[1] / ny
    dz = cell[2] / nz
    
    print(f"  > Cutting {os.path.basename(source_file)}")
    print(f"    Box: {lengths}, Grid: {nx}x{ny}x{nz}, Sub-box size approx: {lengths/([nx,ny,nz])}")
    
    saved_count = 0
    base_name = os.path.splitext(os.path.basename(source_file))[0]
    
    # 为了处理跨边界原子，我们可以简单地将原子的分数坐标映射到 [0, 1] 区间
    # 然后根据分数坐标判断属于哪个 grid
    
    scaled = atoms.get_scaled_positions(wrap=True)
    
    for ix in range(nx):
        x_min, x_max = ix / nx, (ix + 1) / nx
        
        for iy in range(ny):
            y_min, y_max = iy / ny, (iy + 1) / ny
            
            for iz in range(nz):
                z_min, z_max = iz / nz, (iz + 1) / nz
                
                # 筛选原子
                mask = (
                    (scaled[:, 0] >= x_min) & (scaled[:, 0] < x_max) &
                    (scaled[:, 1] >= y_min) & (scaled[:, 1] < y_max) &
                    (scaled[:, 2] >= z_min) & (scaled[:, 2] < z_max)
                )
                
                sub_indices = np.where(mask)[0]
                n_atoms = len(sub_indices)
                
                if atom_limit[0] <= n_atoms <= atom_limit[1]:
                    # 构建子结构
                    # 1. 提取原子
                    sub_pos_scaled = scaled[sub_indices].copy()
                    
                    # 2. 将坐标平移回 [0, 1] 相对子盒子
                    # sub_x = (global_x - x_min) * nx
                    sub_pos_scaled[:, 0] = (sub_pos_scaled[:, 0] - x_min) * nx
                    sub_pos_scaled[:, 1] = (sub_pos_scaled[:, 1] - y_min) * ny
                    sub_pos_scaled[:, 2] = (sub_pos_scaled[:, 2] - z_min) * nz
                    
                    # 3. 定义新 Cell
                    new_cell = np.array([dx, dy, dz])
                    
                    # 4. 创建 Atoms 对象
                    # 注意：这里需要还原为笛卡尔坐标
                    sub_pos_cart = np.dot(sub_pos_scaled, new_cell)
                    
                    sub_atoms = Atoms(
                        symbols=np.array(atoms.get_chemical_symbols())[sub_indices],
                        positions=sub_pos_cart,
                        cell=new_cell,
                        pbc=[True, True, True]
                    )
                    
                    # 微扰
                    if perturb_amp > 0:
                        rng = np.random.default_rng()
                        disp = (rng.random((len(sub_atoms), 3)) - 0.5) * 2.0 * perturb_amp
                        sub_atoms.positions += disp
                        sub_atoms.wrap()
                    
                    # 保存
                    out_name = f"{base_name}_cut_{ix}_{iy}_{iz}_N{n_atoms}.vasp"
                    out_path = os.path.join(output_base_dir, out_name)
                    write_vasp_compat(sub_atoms, out_path)
                    saved_count += 1
                    
    print(f"    Saved {saved_count} blocks.")
    return saved_count


def run_from_config(config, out_root_dir):
    """
    主执行入口
    """
    if not config.get("enable", False):
        print("[GridCutter] Task disabled.")
        return

    source_dir = config.get("source_dir", "twin_complex")
    output_subdir = config.get("output_subdir", "20_grid_cuts")
    
    # 绝对路径处理
    # 假设脚本在项目根目录运行
    if not os.path.isabs(source_dir):
        source_dir = os.path.join(os.getcwd(), source_dir)
        
    out_dir = os.path.join(out_root_dir, output_subdir)
    ensure_dir(out_dir)
    
    print(f"\n>>> [GridCutter] Starting...")
    print(f"    Source: {source_dir}")
    print(f"    Output: {out_dir}")
    
    if not os.path.exists(source_dir):
        print(f"    [Error] Source directory does not exist!")
        return

    total_saved = 0
    
    # 递归遍历所有文件
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            # 简单判断是否为结构文件
            if file.endswith(".vasp") or "POSCAR" in file or file.isdigit(): 
                # file.isdigit() 是为了兼容那些名为 "1", "2" 的文件
                full_path = os.path.join(root, file)
                total_saved += cut_single_file(full_path, out_dir, config)
                
    print(f">>> [GridCutter] Done. Total blocks saved: {total_saved}\n")
