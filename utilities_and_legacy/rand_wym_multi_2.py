import numpy as np
from ase.io import read, write
from ase.build import make_supercell
from pathlib import Path
import argparse

# ================== 默认参数 ==================
SUPERCELL_MATRIX = [[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]]
STRAINS = [1.0]  # 如果不需要晶胞形变，只保留1.0
RATTLE_STEPS = 3
RATTLE_STD = 0.02
RANDOM_SEED = 1234
OUTPUT_FORMAT = 'xyz'
# ============================================

def generate_perturbed_structures(atoms, rattle_steps, rattle_std, rng):
    """对单个 ASE atoms 对象生成扰动结构列表"""
    structures = []

    # 建立超胞（默认1x1x1不改变）
    supercell_atoms = make_supercell(atoms, SUPERCELL_MATRIX)

    # 原始结构加入列表
    structures.append(supercell_atoms.copy())

    # Rattle 扰动
    for i in range(rattle_steps):
        rattled = supercell_atoms.copy()
        stdev = (i + 1) * rattle_std
        rattled.rattle(stdev=stdev, rng=rng)
        structures.append(rattled)

    return structures

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="读取 XYZ 文件中的所有结构并进行随机扰动")
    parser.add_argument("input", type=str, help="输入 XYZ 文件")
    parser.add_argument("output", type=str, help="输出文件路径")
    parser.add_argument("--rattle_steps", type=int, default=RATTLE_STEPS, help="Rattle 扰动次数")
    parser.add_argument("--rattle_std", type=float, default=RATTLE_STD, help="Rattle 标准差")
    parser.add_argument("--seed", type=int, default=RANDOM_SEED, help="随机种子")
    args = parser.parse_args()

    input_file = Path(args.input)
    if not input_file.exists():
        print(f"错误: 找不到输入文件 '{args.input}'")
        exit(1)

    rng = np.random.RandomState(args.seed)
    total_structures = []

    try:
        # 读取 XYZ 文件中所有结构
        atoms_list = read(str(input_file), index=":")  # index=":" 读取所有 frames
        print(f"读取到 {len(atoms_list)} 个结构")

        for atoms in atoms_list:
            perturbed = generate_perturbed_structures(atoms, args.rattle_steps, args.rattle_std, rng)
            total_structures.extend(perturbed)

    except Exception as e:
        print(f"读取或处理文件 {input_file.name} 出错: {e}")
        exit(1)

    # 写入输出文件
    write(args.output, total_structures, format=OUTPUT_FORMAT)
    print(f"生成 {len(total_structures)} 个扰动结构，保存到 {args.output}")

