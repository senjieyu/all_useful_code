#!/usr/bin/env python
import re
import os
import sys
import numpy as np
from ase.data import chemical_symbols


def parse_lattice(lattice_str):
    """把 Lattice=\"...\" 里的 9 个数字解析成 3x3 矩阵"""
    nums = [float(x) for x in lattice_str.strip().split()]
    if len(nums) != 9:
        raise ValueError(f"Lattice 不是 9 个数: {lattice_str}")
    return np.array(nums).reshape(3, 3)


def parse_n_floats(s, n):
    """解析包含 n 个浮点数的字符串"""
    nums = [float(x) for x in s.strip().split()]
    if len(nums) != n:
        raise ValueError(f"期望 {n} 个数，实际 {len(nums)}: {s}")
    return nums


def build_type_map(all_elements):
    """按 ASE 的 chemical_symbols 顺序给元素编号"""
    elements = sorted(all_elements, key=lambda x: chemical_symbols.index(x))
    return {e: i for i, e in enumerate(elements)}


def convert(fin, fout=None):
    if fout is None:
        root, _ = os.path.splitext(fin)
        fout = root + ".cfg"

    # ---------- 第一遍扫描：收集所有元素种类 ----------
    all_elements = set()

    with open(fin, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if line == "":
                continue

            # 这一行应该是原子数 nat
            try:
                nat = int(line)
            except ValueError:
                raise RuntimeError(f"这一行不是原子数: {line}")

            header = f.readline()
            if not header:
                raise RuntimeError("文件在 header 处截断")

            # 跳过 nat 行原子，只统计元素
            for _ in range(nat):
                atom_line = f.readline()
                if not atom_line:
                    raise RuntimeError("原子行数量不足")
                sp = atom_line.split()
                all_elements.add(sp[0])

    type_map = build_type_map(all_elements)
    print("元素到 type 映射:", type_map)

    # ---------- 第二遍：真正转换并写 cfg ----------
    with open(fin, "r") as f_in, open(fout, "w") as f_out:
        frame = 0
        while True:
            line = f_in.readline()
            if not line:
                break
            line = line.strip()
            if line == "":
                continue

            # 这一行：nat
            try:
                nat = int(line)
            except ValueError:
                raise RuntimeError(f"第 {frame} 个结构：这一行不是原子数: {line}")

            # header 行
            header = f_in.readline()
            if not header:
                raise RuntimeError("文件在 header 处截断")
            header = header.strip()

            # ---- Lattice ----
            lat_match = re.search(r'Lattice="([^"]+)"', header)
            if not lat_match:
                raise RuntimeError(f"第 {frame} 个结构：找不到 Lattice=\"...\"")
            lattice = parse_lattice(lat_match.group(1))
            volume = abs(np.linalg.det(lattice))  # 体积

            # ---- energy ----
            en_match = re.search(r'energy=([^\s]+)', header)
            if not en_match:
                raise RuntimeError(f"第 {frame} 个结构：找不到 energy=")
            energy = float(en_match.group(1))

            # ---- stress ----
            st_match = re.search(r'stress="([^"]+)"', header)
            if not st_match:
                raise RuntimeError(f"第 {frame} 个结构：找不到 stress=\"...\"")
            stress_mat = np.array(parse_n_floats(st_match.group(1), 9)).reshape(3, 3)

            # ⭐ 关键：virial = - stress * volume （单位 eV）
            virial_mat = -stress_mat * volume

            # ---- 读原子行 ----
            species = []
            pos = []
            forces = []
            for _ in range(nat):
                atom_line = f_in.readline()
                if not atom_line:
                    raise RuntimeError(f"第 {frame} 个结构：原子行不足")
                sp = atom_line.split()
                s = sp[0]
                x, y, z = map(float, sp[1:4])
                fx, fy, fz = map(float, sp[4:7])
                species.append(s)
                pos.append([x, y, z])
                forces.append([fx, fy, fz])

            pos = np.array(pos)
            forces = np.array(forces)

            # ---- 写 cfg ----
            f_out.write("BEGIN_CFG\n")
            f_out.write(" Size\n")
            f_out.write(f"  {nat:6d} \n")
            f_out.write(" Supercell \n")
            f_out.write("{:15.10f} {:15.10f} {:15.10f}\n".format(*lattice[0]))
            f_out.write("{:15.10f} {:15.10f} {:15.10f}\n".format(*lattice[1]))
            f_out.write("{:15.10f} {:15.10f} {:15.10f}\n".format(*lattice[2]))
            f_out.write(
                "AtomData:  id type       cartes_x      cartes_y      cartes_z     fx          fy          fz\n"
            )
            for i in range(nat):
                f_out.write(
                    " {:6d} {:6d} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(
                        i + 1,
                        type_map[species[i]],
                        pos[i, 0],
                        pos[i, 1],
                        pos[i, 2],
                        forces[i, 0],
                        forces[i, 1],
                        forces[i, 2],
                    )
                )
            f_out.write("Energy \n")
            f_out.write(f"\t{energy} \n")
            f_out.write("PlusStress:  xx          yy          zz          yz          xz          xy \n")
            # 写 virial 的 6 个分量 (xx, yy, zz, yz, xz, xy)
            f_out.write(
                f"\t{virial_mat[0,0]}  \t{virial_mat[1,1]}  \t{virial_mat[2,2]}  "
                f"\t{virial_mat[1,2]}  \t{virial_mat[0,2]}  \t{virial_mat[0,1]} \n"
            )
            f_out.write("END_CFG \n")

            frame += 1

    print(f"已生成: {fout}，共处理 {frame} 个结构")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python extxyz2cfg_direct.py data.extxyz [train.cfg]")
        sys.exit(1)

    fin = sys.argv[1]
    fout = sys.argv[2] if len(sys.argv) >= 3 else None
    convert(fin, fout)

