#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
use : python merge.py dft.xyz fccdft.xyz out.xyz
"""
import sys
from pathlib import Path

def count_xyz_frames(xyz_path: Path) -> int:
    """
    统计多帧 XYZ 文件中的结构(帧)数量。
    假设每一帧标准 XYZ 格式：
      第1行: 原子数 N
      第2行: 注释行
      接下来 N 行: 原子信息
    """
    count = 0
    with xyz_path.open("r", encoding="utf-8", errors="ignore") as f:
        while True:
            # 读第一行（原子数）
            line = f.readline()
            if not line:
                break  # EOF
            line = line.strip()
            if line == "":
                # 跳过空行
                continue

            try:
                n_atoms = int(line)
            except ValueError:
                raise ValueError(
                    f"[{xyz_path}] 解析失败：期望原子数整数，但读到：{line!r}"
                )

            # 读注释行
            comment = f.readline()
            if not comment:
                raise ValueError(f"[{xyz_path}] 文件不完整：缺少注释行")

            # 跳过 n_atoms 行的原子坐标
            for _ in range(n_atoms):
                atom_line = f.readline()
                if not atom_line:
                    raise ValueError(
                        f"[{xyz_path}] 文件不完整：原子坐标行不足 (需要 {n_atoms} 行)"
                    )

            count += 1

    return count


def merge_xyz(in1: Path, in2: Path, out: Path):
    """
    合并两个 xyz 文件到 out（顺序：in1 再 in2）
    """
    with out.open("w", encoding="utf-8") as fout:
        for p in (in1, in2):
            with p.open("r", encoding="utf-8", errors="ignore") as fin:
                for line in fin:
                    fout.write(line)


def main():
    if len(sys.argv) != 4:
        print("用法: python merge.py dft.xyz fccdft.xyz out.xyz")
        sys.exit(1)

    in1 = Path(sys.argv[1])
    in2 = Path(sys.argv[2])
    out = Path(sys.argv[3])

    if not in1.exists():
        print(f"找不到文件: {in1}")
        sys.exit(1)
    if not in2.exists():
        print(f"找不到文件: {in2}")
        sys.exit(1)

    # 统计结构数
    n1 = count_xyz_frames(in1)
    n2 = count_xyz_frames(in2)

    # 合并
    merge_xyz(in1, in2, out)

    # 合并后再统计一次（更保险）
    n_out = count_xyz_frames(out)

    print(f"{in1.name} 结构数: {n1}")
    print(f"{in2.name} 结构数: {n2}")
    print(f"{out.name} 合并后结构数: {n_out}")


if __name__ == "__main__":
    main()

