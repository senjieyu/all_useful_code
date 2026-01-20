# ============================================
# 用法（Usage）:
#   python batch_cif2vasp.py 输入文件夹 输出文件夹
#
# 示例:
#   python batch_cif2vasp.py ./cif_data ./vasp_out
#
# 功能:
#   1) 递归扫描输入文件夹及其所有子文件夹
#   2) 找到所有 .cif 文件
#   3) 转换为 VASP POSCAR 格式（.vasp）
#   4) 输出到新的文件夹，并保持原有子目录结构
# ============================================

import os
import sys
from ase.io import read, write

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main():
    if len(sys.argv) < 3:
        print("Usage: python batch_cif2vasp.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.isdir(input_dir):
        print(f"Error: input_dir not found -> {input_dir}")
        sys.exit(1)

    count = 0
    fail = 0

    for root, dirs, files in os.walk(input_dir):
        for fname in files:
            if fname.lower().endswith(".cif"):
                cif_path = os.path.join(root, fname)

                # 保持原目录结构
                rel_path = os.path.relpath(root, input_dir)
                out_subdir = os.path.join(output_dir, rel_path)
                mkdir(out_subdir)

                # 输出文件名：同名改后缀
                base = os.path.splitext(fname)[0]
                out_path = os.path.join(out_subdir, base + ".vasp")

                try:
                    atoms = read(cif_path)
                    write(out_path, atoms, format="vasp")
                    count += 1
                    print(f"[OK] {cif_path}  ->  {out_path}")
                except Exception as e:
                    fail += 1
                    print(f"[FAIL] {cif_path}  ({e})")

    print("\n========== Done ==========")
    print(f"Converted: {count}")
    print(f"Failed   : {fail}")
    print(f"Output dir: {output_dir}")

if __name__ == "__main__":
    main()
