#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import shutil
import subprocess
from pathlib import Path
import sys


def run(cmd):
    cmd = [str(x) for x in cmd]
    print(">>", " ".join(cmd))
    subprocess.run(cmd, check=True)


def estimate_atoms_fcc(lx, ly, lz, a0):
    return 4.0 * lx * ly * lz / (a0 ** 3)


def generate_hex_nodes(nx_cols, ny_rows, d, lz, angle_pool, x_shift=0.0, y_shift=0.0):
    """
    Generate nodes on a triangular lattice so that Voronoi cells are hexagon-like.
    """
    dx = d
    dy = math.sqrt(3.0) * 0.5 * d

    nodes = []
    for row in range(ny_rows):
        y = row * dy + y_shift
        x_offset = 0.5 * dx if (row % 2 == 1) else 0.0

        for col in range(nx_cols):
            x = col * dx + x_offset + x_shift
            ang = angle_pool[(row + col) % len(angle_pool)]
            nodes.append((x, y, 0.5 * lz, ang))

    return nodes, dx, dy


def main():
    atomsk = shutil.which("atomsk")
    if atomsk is None:
        atomsk = r"D:\2025.916aftersoftware\atmosk\Atomsk\atomsk.EXE"

    if not Path(atomsk).exists():
        print("Error: Atomsk not found.")
        sys.exit(1)

    workdir = Path("atomsk_hex_poly_regular").resolve()
    workdir.mkdir(exist_ok=True)

    # -----------------------------
    # Basic parameters
    # -----------------------------
    a0 = 3.615  # Cu lattice constant (Angstrom)

    # 控制晶粒数：列数、行数
    nx_cols = 4
    ny_rows = 4

    # 晶粒中心距
    d = 20.0
    dx = d
    dy = math.sqrt(3.0) * 0.5 * d

    # 盒子尺寸严格匹配三角晶格
    # 这样周期复制更规整
    Lx = nx_cols * dx
    Ly = ny_rows * dy
    Lz = 18.0

    # 关键：整体平移，让边界切在更对称的位置
    # 比之前更规整
    x_shift = 0.25 * dx
    y_shift = 0.50 * dy

    angle_pool = [0, 30, 60, 90]

    est_n = estimate_atoms_fcc(Lx, Ly, Lz, a0)
    print(f"Estimated atom count before GB processing: {est_n:.0f}")

    seed_file = workdir / "seed.xsf"
    poly_file = workdir / "poly_hex_regular.txt"
    poly_cfg = workdir / "hex_poly_regular.cfg"
    wrapped_cfg = workdir / "hex_poly_regular_wrapped.cfg"

    # -----------------------------
    # 1) Create seed crystal
    # -----------------------------
    run([
        atomsk,
        "--create", "fcc", a0, "Cu",
        "orient", "[100]", "[010]", "[001]",
        seed_file
    ])

    # -----------------------------
    # 2) Generate nodes
    # -----------------------------
    nodes, dx, dy = generate_hex_nodes(
        nx_cols=nx_cols,
        ny_rows=ny_rows,
        d=d,
        lz=Lz,
        angle_pool=angle_pool,
        x_shift=x_shift,
        y_shift=y_shift
    )

    # -----------------------------
    # 3) Write polycrystal file
    # -----------------------------
    lines = [f"box {Lx:.6f} {Ly:.6f} {Lz:.6f}", ""]

    for x, y, z, ang in nodes:
        # 把节点包回盒子，避免 seed 本身跑出边界
        x = x % Lx
        y = y % Ly
        lines.append(f"node {x:.6f} {y:.6f} {z:.6f}   0   0   {ang}")

    lines.append("")
    poly_file.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote: {poly_file}")
    print(f"Number of nodes: {len(nodes)}")
    print(f"Box: {Lx:.3f} x {Ly:.3f} x {Lz:.3f} Angstrom")

    # -----------------------------
    # 4) Build polycrystal
    # -----------------------------
    run([
        atomsk,
        "--polycrystal",
        seed_file,
        poly_file,
        poly_cfg
    ])

    # -----------------------------
    # 5) Wrap back to main cell
    # -----------------------------
    run([
        atomsk,
        poly_cfg,
        "-wrap",
        wrapped_cfg
    ])

    # -----------------------------
    # 6) Export to LAMMPS
    # -----------------------------
    run([
        atomsk,
        wrapped_cfg,
        "lammps"
    ])

    print("\nDone.")
    print("Generated files:")
    for p in sorted(workdir.iterdir()):
        print(" -", p.name)

    print("\nTips:")
    print("1. In OVITO, use 'Wrap at periodic boundaries'")
    print("2. Then use 'Replicate' in x/y to verify periodic tiling")
    print("3. LAMMPS should use: boundary p p p")


if __name__ == "__main__":
    main()
