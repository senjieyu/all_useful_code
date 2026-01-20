#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
将 LAMMPS dump custom (单帧) 转成 VASP POSCAR (.vasp)

用法：
  python dump_to_vasp.py <dump文件夹> <输出文件夹>

示例：
  python dump_to_vasp.py analysis_results/accept_dumps vasp_accept
"""

import os
import sys
import glob
import numpy as np


# ===================== 用户可改参数 =====================
ELEMENT = "Cu"     # 元素名（当前是 Cu）
DIRECT = True      # True -> Direct 坐标；False -> Cartesian
# ======================================================


def read_dump(path):
    """读取 dump 文件，返回 box, positions"""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # ---- timestep 无关，直接顺序解析 ----
    natoms = None
    box = []
    header = None
    data_start = None

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("ITEM: NUMBER OF ATOMS"):
            natoms = int(lines[i + 1].strip())
            i += 2
            continue

        if line.startswith("ITEM: BOX BOUNDS"):
            box = []
            for j in range(3):
                lo, hi = map(float, lines[i + 1 + j].split()[:2])
                box.append((lo, hi))
            i += 4
            continue

        if line.startswith("ITEM: ATOMS"):
            header = line.split()[2:]
            data_start = i + 1
            break

        i += 1

    if natoms is None or header is None or not box:
        raise RuntimeError(f"{path}: dump 文件格式不完整")

    # ---- 读取原子数据 ----
    data = []
    for j in range(data_start, data_start + natoms):
        data.append(lines[j].split())

    arr = np.array(data, dtype=float)

    idx = {k: i for i, k in enumerate(header)}
    x = arr[:, idx["x"]]
    y = arr[:, idx["y"]]
    z = arr[:, idx["z"]]

    # ---- box 向量 ----
    lx = box[0][1] - box[0][0]
    ly = box[1][1] - box[1][0]
    lz = box[2][1] - box[2][0]

    cell = np.array([
        [lx, 0.0, 0.0],
        [0.0, ly, 0.0],
        [0.0, 0.0, lz],
    ])

    pos_cart = np.column_stack((x, y, z))

    return cell, pos_cart


def write_poscar(path, cell, pos_cart):
    """写 POSCAR"""
    if DIRECT:
        pos = pos_cart @ np.linalg.inv(cell)
        coord_type = "Direct"
    else:
        pos = pos_cart
        coord_type = "Cartesian"

    with open(path, "w") as f:
        f.write(f"{ELEMENT} converted from LAMMPS dump\n")
        f.write("1.0\n")
        for v in cell:
            f.write(f"{v[0]:16.10f} {v[1]:16.10f} {v[2]:16.10f}\n")
        f.write(f"{ELEMENT}\n")
        f.write(f"{len(pos)}\n")
        f.write(f"{coord_type}\n")
        for p in pos:
            f.write(f"{p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    dump_dir = sys.argv[1]
    out_dir = sys.argv[2]

    os.makedirs(out_dir, exist_ok=True)

    dump_files = sorted(glob.glob(os.path.join(dump_dir, "*.force.dump")))
    if not dump_files:
        sys.exit(f"❌ 在 {dump_dir} 中找不到 *.force.dump")

    for fp in dump_files:
        name = os.path.basename(fp).replace(".force.dump", "")
        out = os.path.join(out_dir, f"{name}.vasp")

        cell, pos = read_dump(fp)
        write_poscar(out, cell, pos)

        print(f"✔ {name}.vasp")

    print(f"\n✅ 共转换 {len(dump_files)} 个结构，输出在 {out_dir}/")


if __name__ == "__main__":
    main()

