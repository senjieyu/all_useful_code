# ============================================
# 用法（Usage）:
#   python xyz2poscar.py input.xyz
#
# 功能:
#   1) 读取 input.xyz（或其他 ASE 支持的轨迹文件，多帧结构）
#   2) 输出到文件夹 xyz2poscar/
#   3) 每个结构保存为: 化学式_0001.vasp, 化学式_0002.vasp ...
#   4) 生成 summary.txt：
#        - 先列出有哪些化学式
#        - 再统计每种化学式对应的结构数量
#
# 示例:
#   python xyz2poscar.py traj.xyz
# ============================================

from ase.io import iread, write
import os
import sys
from collections import Counter, defaultdict

def mkdir(d):
    if not os.path.exists(d):
        os.mkdir(d)

def main():
    if len(sys.argv) < 2:
        print("Usage: python xyz2poscar.py <input_file>")
        sys.exit(1)

    infile = sys.argv[1]
    atoms_list = list(iread(infile))

    outdir = "xyz2poscar"
    mkdir(outdir)

    # 统计：化学式 -> 结构数
    formula_counts = Counter()
    # 记录每个化学式已经写了多少个，用于编号
    formula_written = defaultdict(int)

    # 先算出所有化学式（保证 summary “先列出有哪些结构”）
    formulas = []
    for at in atoms_list:
        # 金属体系推荐 mode='metal'；不想要可改为 at.get_chemical_formula()
        f = at.get_chemical_formula(mode="metal")
        formulas.append(f)
        formula_counts[f] += 1

    # 写 POSCAR (.vasp)
    for at, f in zip(atoms_list, formulas):
        formula_written[f] += 1
        idx = formula_written[f]
        # 例：TiAl_0001.vasp
        filename = f"{f}_{idx:04d}.vasp"
        path = os.path.join(outdir, filename)
        write(path, at, format="vasp")

    # 写 summary.txt
    summary_path = os.path.join(outdir, "summary.txt")
    unique_formulas = sorted(formula_counts.keys())

    with open(summary_path, "w", encoding="utf-8") as fp:
        fp.write("=== 化学式列表（有哪些结构）===\n")
        fp.write(", ".join(unique_formulas) + "\n\n")

        fp.write("=== 各化学式结构数目统计 ===\n")
        for f in unique_formulas:
            fp.write(f"{f}\t{formula_counts[f]}\n")

        fp.write("\n=== 总结构数 ===\n")
        fp.write(str(len(atoms_list)) + "\n")

    print(f"Done. Wrote {len(atoms_list)} structures to '{outdir}/' and summary to '{summary_path}'.")

if __name__ == "__main__":
    main()
