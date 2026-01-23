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
#       则按“逐步试探”扩胞：
#         1) 找当前最短边，把该方向倍数设为 2（即原始边长*2）
#         2) 重新比较：
#              - 若仍是同一方向最短边，则该方向倍数设为 3、4、5...（原始边长*3,*4...）
#              - 若最短边换了方向，则对新的最短边先设为 2
#         3) 循环直到：
#              (a) min(a,b,c) >= 10Å  或
#              (b) 下一步扩胞会使总原子数 > 200（硬约束）
#   - 输出:
#       1) 输出最终结构到新文件夹（保持子目录结构）
#       2) expansion_log.csv 记录扩胞前后变化/原因
# 依赖:
#   pip install ase
# ============================================

import os
import sys
import csv
from ase.io import read, write

TARGET = 10.0
MAX_ATOMS = 200

def mkdir(path: str):
    if not os.path.exists(path):
        os.makedirs(path)

def get_cell_lengths(atoms):
    a, b, c = atoms.cell.lengths()
    return float(a), float(b), float(c)

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

            # 判断是否需要扩胞：晶格参数有小于10 且 原子数<100
            need_expand = (min(a0, b0, c0) < TARGET) and (natoms0 < 100)

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

            # ========= 新逻辑：逐步试探扩胞 =========
            # m[i] 表示该方向相对“原胞”的倍数（原始边长 * m[i]）
            m = [1, 1, 1]     # [n1, n2, n3]
            last_idx = None
            stop_reason = ""

            # 记录“当前最好（合法<=200）”的结构
            atoms_best = atoms.copy()

            while True:
                mult = m[0] * m[1] * m[2]
                nat = natoms0 * mult

                # 当前结构
                atoms_cur = atoms.repeat((m[0], m[1], m[2]))
                a, b, c = get_cell_lengths(atoms_cur)
                L = [a, b, c]
                minL = min(L)

                # 停止条件 1：达标
                if minL >= TARGET - 1e-8:
                    atoms_best = atoms_cur
                    stop_reason = "reached_min_length>=10"
                    break

                # 找当前最短边（并列时取索引最小者：a优先于b优先于c）
                idx = min(range(3), key=lambda i: L[i])

                # 生成下一步倍数（关键规则）
                m_next = m.copy()
                if idx == last_idx:
                    # 仍是同一条最短边：2->3->4...
                    m_next[idx] = m[idx] + 1
                else:
                    # 换最短边：先设为2；如果已>=2（可能由于反复切换/并列），则+1推进避免卡住
                    m_next[idx] = 2 if m[idx] < 2 else (m[idx] + 1)

                mult_next = m_next[0] * m_next[1] * m_next[2]
                nat_next = natoms0 * mult_next

                # 停止条件 2：下一步会超原子数上限 -> 停在当前
                if nat_next > MAX_ATOMS:
                    atoms_best = atoms_cur
                    stop_reason = f"hit_atoms_limit(next_mult={mult_next}_natoms={nat_next})"
                    break

                # 接受下一步
                m = m_next
                last_idx = idx
                atoms_best = atoms_cur  # 当前仍合法，更新“最好结构”

            # 最终结构信息
            n1, n2, n3 = m[0], m[1], m[2]
            a1, b1, c1 = get_cell_lengths(atoms_best)
            natoms1 = len(atoms_best)

            # 输出最终结构（可能因200原子限制未完全达标，但已是规则下的最优停止点）
            write(out_path, atoms_best, format="vasp", vasp5=True, direct=True)

            status = "EXPANDED_OK" if min(a1, b1, c1) >= TARGET - 1e-8 else "EXPANDED_STOP_ATOMS"
            expanded += 1
            log_rows.append({
                "file": os.path.relpath(in_path, input_dir),
                "status": status,
                "reason": stop_reason,
                "n1": n1, "n2": n2, "n3": n3,
                "a0": f"{a0:.6f}", "b0": f"{b0:.6f}", "c0": f"{c0:.6f}",
                "natoms0": natoms0,
                "a1": f"{a1:.6f}", "b1": f"{b1:.6f}", "c1": f"{c1:.6f}",
                "natoms1": natoms1
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
