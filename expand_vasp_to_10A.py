# ============================================
# 用法（Usage）:
#   python expand_vasp_to_10A.py 输入文件夹 输出文件夹
#
# 示例:
#   python expand_vasp_to_10A.py ./vasp_in ./vasp_out
#
# 说明:
#   - 递归扫描输入文件夹及子文件夹中的所有 .vasp 文件
#   - 若满足: (任一晶格边长 < 10Å) 且 (原子数 < 100)
#       则尝试扩胞，使扩胞后 a,b,c 均 >= 10Å
#       且扩胞后总原子数 <= 200（硬约束）
#       在原子数<=200前提下，三边尽可能接近10Å（取最小需要倍数）
#   - 输出:
#       1) 扩胞后的 .vasp 输出到新文件夹（保持子目录结构）
#       2) expansion_log.csv 记录扩胞前后变化/原因
# 依赖:
#   pip install ase
# ============================================

import os
import sys
import math
import csv
from ase.io import read, write

def mkdir(path: str):
    if not os.path.exists(path):
        os.makedirs(path)

def get_cell_lengths(atoms):
    # 返回 a,b,c（单位 Å）
    a, b, c = atoms.cell.lengths()
    return float(a), float(b), float(c)

def required_multiplier(L: float, target: float = 10.0) -> int:
    # 让 n*L >= target 的最小整数 n
    if L >= target:
        return 1
    return int(math.ceil(target / L))

def main():
    if len(sys.argv) < 3:
        print("Usage: python expand_vasp_to_10A.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.isdir(input_dir):
        print(f"Error: input_dir not found -> {input_dir}")
        sys.exit(1)

    mkdir(output_dir)

    log_rows = []
    total = 0
    expanded = 0
    skipped = 0
    failed = 0

    for root, _, files in os.walk(input_dir):
        for fname in files:
            if not fname.lower().endswith(".vasp"):
                continue

            total += 1
            in_path = os.path.join(root, fname)

            # 输出保持目录结构
            rel_dir = os.path.relpath(root, input_dir)
            out_subdir = os.path.join(output_dir, rel_dir)
            mkdir(out_subdir)
            out_path = os.path.join(out_subdir, fname)

            try:
                atoms = read(in_path, format="vasp")
            except Exception as e:
                failed += 1
                log_rows.append({
                    "file": os.path.relpath(in_path, input_dir),
                    "status": "FAIL_READ",
                    "reason": str(e),
                    "n1": "", "n2": "", "n3": "",
                    "a0": "", "b0": "", "c0": "",
                    "natoms0": "",
                    "a1": "", "b1": "", "c1": "",
                    "natoms1": ""
                })
                continue

            a0, b0, c0 = get_cell_lengths(atoms)
            natoms0 = len(atoms)

            # 默认：先把原结构也复制到输出目录？（这里选择：不扩胞的也输出原样，方便统一使用）
            # 如果你希望“只输出扩胞后的”，可以把下面 write(out_path, atoms...) 移到扩胞分支里。
            # write(out_path, atoms, format="vasp", vasp5=True, direct=True)

            # 判断是否需要扩胞：晶格参数有小于10 且 原子数<100
            need_expand = (min(a0, b0, c0) < 10.0) and (natoms0 < 100)

            if not need_expand:
                # 直接原样输出
                write(out_path, atoms, format="vasp", vasp5=True, direct=True)
                skipped += 1
                log_rows.append({
                    "file": os.path.relpath(in_path, input_dir),
                    "status": "SKIP",
                    "reason": "no_need_or_natoms>=100",
                    "n1": 1, "n2": 1, "n3": 1,
                    "a0": f"{a0:.6f}", "b0": f"{b0:.6f}", "c0": f"{c0:.6f}",
                    "natoms0": natoms0,
                    "a1": f"{a0:.6f}", "b1": f"{b0:.6f}", "c1": f"{c0:.6f}",
                    "natoms1": natoms0
                })
                continue

            # 计算每个方向所需的最小扩胞倍数（保证 >=10Å 且尽量接近10Å）
            n1 = required_multiplier(a0, 10.0)
            n2 = required_multiplier(b0, 10.0)
            n3 = required_multiplier(c0, 10.0)

            mult = n1 * n2 * n3
            natoms1 = natoms0 * mult

            # 原子数硬约束
            if natoms1 > 200:
                # 无解：因为要保证三边>=10，n1/n2/n3 都是最小需求；再减会不达标，再增只会更糟
                write(out_path, atoms, format="vasp", vasp5=True, direct=True)
                skipped += 1
                log_rows.append({
                    "file": os.path.relpath(in_path, input_dir),
                    "status": "SKIP",
                    "reason": f"would_exceed_200_atoms(min_mult={mult})",
                    "n1": n1, "n2": n2, "n3": n3,
                    "a0": f"{a0:.6f}", "b0": f"{b0:.6f}", "c0": f"{c0:.6f}",
                    "natoms0": natoms0,
                    "a1": "", "b1": "", "c1": "",
                    "natoms1": natoms1
                })
                continue

            # 执行扩胞（对角扩胞，不引入剪切）
            atoms_sc = atoms.repeat((n1, n2, n3))
            a1, b1, c1 = get_cell_lengths(atoms_sc)

            # 再保险检查
            if min(a1, b1, c1) < 10.0 - 1e-8:
                # 理论上不会发生，除非 cell lengths 异常
                write(out_path, atoms, format="vasp", vasp5=True, direct=True)
                skipped += 1
                log_rows.append({
                    "file": os.path.relpath(in_path, input_dir),
                    "status": "SKIP",
                    "reason": "post_check_failed(length<10)",
                    "n1": n1, "n2": n2, "n3": n3,
                    "a0": f"{a0:.6f}", "b0": f"{b0:.6f}", "c0": f"{c0:.6f}",
                    "natoms0": natoms0,
                    "a1": f"{a1:.6f}", "b1": f"{b1:.6f}", "c1": f"{c1:.6f}",
                    "natoms1": len(atoms_sc)
                })
                continue

            # 输出扩胞结构
            write(out_path, atoms_sc, format="vasp", vasp5=True, direct=True)
            expanded += 1
            log_rows.append({
                "file": os.path.relpath(in_path, input_dir),
                "status": "EXPANDED",
                "reason": "",
                "n1": n1, "n2": n2, "n3": n3,
                "a0": f"{a0:.6f}", "b0": f"{b0:.6f}", "c0": f"{c0:.6f}",
                "natoms0": natoms0,
                "a1": f"{a1:.6f}", "b1": f"{b1:.6f}", "c1": f"{c1:.6f}",
                "natoms1": len(atoms_sc)
            })

    # 写 CSV 日志
    log_path = os.path.join(output_dir, "expansion_log.csv")
    fieldnames = [
        "file", "status", "reason",
        "n1", "n2", "n3",
        "a0", "b0", "c0", "natoms0",
        "a1", "b1", "c1", "natoms1"
    ]
    with open(log_path, "w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(log_rows)

    print("========== Done ==========")
    print(f"Total .vasp files : {total}")
    print(f"Expanded          : {expanded}")
    print(f"Skipped           : {skipped}")
    print(f"Failed            : {failed}")
    print(f"Output dir        : {output_dir}")
    print(f"Log CSV           : {log_path}")

if __name__ == "__main__":
    main()
