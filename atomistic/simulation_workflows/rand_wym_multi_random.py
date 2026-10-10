import numpy as np
from ase.io import read, write
from ase.build import make_supercell
from pathlib import Path
import sys

# ==========================================
#               配置区域 (Config)
# ==========================================
INPUT_DIR = "stru"                      # 存放 VASP 文件的输入目录
INPUT_PATTERN = "*.vasp"                # 匹配的文件名规则
OUTPUT_FILE = "out.xyz"                 # 输出的总文件路径
OUTPUT_FORMAT = "extxyz"                # 输出格式

SUPERCELL_MATRIX = [[1, 0, 0],          # 超胞扩建矩阵 (默认1x1x1不扩建)
                    [0, 1, 0], 
                    [0, 0, 1]]

# --- 随机扰动专属配置 ---
PERTURB_NUM = 5                         # 每个文件独立生成的随机扰动结构数量
MAX_DISPLACEMENT = 0.05                 # 原子的最大随机位移幅度 (单位: 埃, Angstrom)
RANDOM_SEED = 1234                      # 随机数种子
# ==========================================


def apply_random_perturbation(atoms, max_disp, rng):
    """
    对所有原子施加均匀的完全随机位置扰动，晶胞大小保持不变。
    """
    perturbed_atoms = atoms.copy()
    pos = perturbed_atoms.get_positions()
    
    # 对 X, Y, Z 三个方向分别生成完全独立的均匀随机分布
    displacements = rng.uniform(-max_disp, max_disp, size=pos.shape)
    
    # 将随机位移加到原子原本的坐标上
    perturbed_atoms.set_positions(pos + displacements)
    return perturbed_atoms


def generate_perturbed_structures(atoms, rng):
    """
    对单个输入结构仅进行随机位置扰动，不改变晶胞体积。
    """
    generated_structures = []
    
    # 建立超胞 (或保持原样)
    supercell_atoms = make_supercell(atoms, SUPERCELL_MATRIX)
    
    # 1. 记录无扰动的基准结构
    generated_structures.append(supercell_atoms.copy())
    
    # 2. 连续生成 PERTURB_NUM 个完全独立的随机扰动结构
    for _ in range(PERTURB_NUM):
        rand_atoms = apply_random_perturbation(supercell_atoms, MAX_DISPLACEMENT, rng)
        generated_structures.append(rand_atoms)
            
    return generated_structures


if __name__ == '__main__':
    input_path = Path(INPUT_DIR)
    
    if not input_path.exists() or not input_path.is_dir():
        print(f"错误: 找不到输入目录 '{INPUT_DIR}'。")
        sys.exit(1)

    vasp_files = list(input_path.glob(INPUT_PATTERN))
    if not vasp_files:
        print(f"警告: 在 '{INPUT_DIR}' 目录中没找到 '{INPUT_PATTERN}' 文件。")
        sys.exit(0)

    # --- 提前打印计算量预估 ---
    structures_per_file = 1 + PERTURB_NUM
    expected_total = len(vasp_files) * structures_per_file
    
    print("==============================")
    print("         预估计算量")
    print("==============================")
    print(f"输入文件: {len(vasp_files)} 个")
    print(f"单文件生成: 1个基准 + {PERTURB_NUM}次随机 = {structures_per_file} 个")
    print(f"预计总生成: {expected_total} 个结构")
    print("==============================\n")

    print(f"--- 开始处理 ---")
    total_structures = []
    stats = {}
    global_rng = np.random.RandomState(RANDOM_SEED)

    for vasp_file in vasp_files:
        try:
            original_atoms = read(vasp_file)
            perturbed_list = generate_perturbed_structures(original_atoms, global_rng)
            
            total_structures.extend(perturbed_list)
            stats[vasp_file.name] = len(perturbed_list)
        except Exception as e:
            print(f"处理文件 {vasp_file.name} 时出错: {e}")

    print("--- 正在写入文件 ---")
    write(OUTPUT_FILE, total_structures, format=OUTPUT_FORMAT)

    print("\n处理完成！最终生成的总结构数:", len(total_structures))

