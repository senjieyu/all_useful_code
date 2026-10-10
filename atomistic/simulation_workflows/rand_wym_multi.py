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
OUTPUT_FORMAT = "extxyz"                # 输出格式 (使用 extxyz 保存晶格和受力/能量等扩展信息)

SUPERCELL_MATRIX = [[1, 0, 0],          # 超胞扩建矩阵 (默认 1x1x1 即不扩建)
                    [0, 1, 0], 
                    [0, 0, 1]]
STRAINS = [1.05, 1.03, 1.0, 0.97, 0.95] # 晶胞缩放形变比例
RATTLE_STEPS = 4                        # 原子位置随机扰动(Rattle)的步数
RATTLE_BASE_STDEV = 0.02                # Rattle 的基础标准差增量
RANDOM_SEED = 1234                      # 随机数种子 (固定种子可确保每次运行生成的结构一致)
# ==========================================


def additional_rand(atom, std, target_elements=["Li", "Na"], seed=1234):
    """
    (可选) 对特定元素施加额外的随机高斯扰动。
    目前在主循环中未调用，可根据需要启用。
    """
    a = atom.arrays['positions']
    s = atom.get_chemical_symbols()
    rng = np.random.RandomState(seed)
    for i in range(len(s)):
        if s[i] in target_elements:
            a[i] = a[i] + rng.normal(scale=std, size=3)
    atom.set_positions(a)


def generate_perturbed_structures(atoms, rng):
    """
    核心逻辑：对单个输入结构进行缩放(Strain)和扰动(Rattle)，返回生成的结构列表。
    """
    generated_structures = []
    
    # 1. 建立超胞 (如果配置为1x1x1则不改变结构)
    supercell_atoms = make_supercell(atoms, SUPERCELL_MATRIX)
    spos = supercell_atoms.get_scaled_positions()
    cell = supercell_atoms.get_cell()
    
    # 2. 遍历形变比例
    for s in STRAINS:
        strained_atoms = supercell_atoms.copy()
        strained_atoms.set_cell(cell * s)
        strained_atoms.set_scaled_positions(spos)
        
        # 将基础形变结构加入列表
        generated_structures.append(strained_atoms)
        
        # 3. 对当前形变结构进行多次 Rattle 扰动
        for i in range(RATTLE_STEPS):
            rattled_atoms = strained_atoms.copy()
            stdev = (i + 1) * RATTLE_BASE_STDEV
            rattled_atoms.rattle(stdev=stdev, rng=rng)
            
            # 如果需要对特定元素额外扰动，可以在这里取消注释：
            # additional_rand(rattled_atoms, std=5 * (i + 1) * 0.01)
            
            generated_structures.append(rattled_atoms)
            
    return generated_structures


if __name__ == '__main__':
    input_path = Path(INPUT_DIR)
    
    # 检查输入目录是否存在
    if not input_path.exists() or not input_path.is_dir():
        print(f"错误: 找不到输入目录 '{INPUT_DIR}'。请确保该文件夹与脚本在同一目录下。")
        sys.exit(1)

    # 获取所有匹配的文件
    vasp_files = list(input_path.glob(INPUT_PATTERN))
    if not vasp_files:
        print(f"警告: 在 '{INPUT_DIR}' 目录中没有找到任何 '{INPUT_PATTERN}' 文件。")
        sys.exit(0)

    print(f"--- 开始处理 ---")
    print(f"输入目录: {INPUT_DIR} (找到 {len(vasp_files)} 个文件)")
    print(f"输出文件: {OUTPUT_FILE}\n")

    total_structures = []
    stats = {}  # 用于记录每个文件生成的结构数量
    
    # 统一使用带有随机种子的生成器，保证科学计算的可重复性
    global_rng = np.random.RandomState(RANDOM_SEED)

    # 遍历处理每个文件
    for vasp_file in vasp_files:
        try:
            # 读取原始结构
            original_atoms = read(vasp_file)
            
            # 生成扰动结构
            perturbed_list = generate_perturbed_structures(original_atoms, global_rng)
            
            # 合并到总列表并记录统计信息
            total_structures.extend(perturbed_list)
            stats[vasp_file.name] = len(perturbed_list)
            
        except Exception as e:
            print(f"读取或处理文件 {vasp_file.name} 时出错: {e}")

    # 将所有生成的结构一次性写入输出文件
    print("--- 正在写入文件 ---")
    write(OUTPUT_FILE, total_structures, format=OUTPUT_FORMAT)

    # 打印最终统计结果
    print("\n==============================")
    print("         处理结果统计")
    print("==============================")
    for filename, count in stats.items():
        print(f"文件: {filename:<20} -> 生成结构数: {count}")
    print("------------------------------")
    print(f"总计生成的结构数量: {len(total_structures)}")
    print("==============================")

