#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import shutil
import subprocess
from pathlib import Path
import sys
import argparse


def run(cmd):
    cmd = [str(x) for x in cmd]
    print(">>", " ".join(cmd))
    subprocess.run(cmd, check=True)


def estimate_atoms_fcc(lx, ly, lz, a0):
    """估算 FCC 体系原子数：N ≈ 4 * V / a0^3"""
    return 4.0 * lx * ly * lz / (a0 ** 3)


def generate_hex_nodes(nx_cols, ny_rows, d, lz, angle_pool, x_shift=0.0, y_shift=0.0):
    """
    生成三角晶格 seed 点。
    Voronoi 剖分后会得到近似六边形/蜂窝状晶粒。

    参数
    ----
    nx_cols, ny_rows : int
        seed 点列数、行数
    d : float
        相邻 seed 的水平间距（Å）
    lz : float
        盒子 z 长度（Å）
    angle_pool : list[float]
        晶粒取向角列表，例如 [0, 30, 60, 90]
    x_shift, y_shift : float
        seed 整体平移量（单位 Å）

    返回
    ----
    nodes : list[(x, y, z, ang)]
    dx, dy : float
        三角晶格在 x / y 方向的步长
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


def build_model(
    atomsk_path,
    workdir_name,
    a0,
    nx_cols,
    ny_rows,
    d,
    lz,
    angle_pool,
    x_shift_factor,
    y_shift_factor,
    seed_orient_x,
    seed_orient_y,
    seed_orient_z,
):
    """
    调用 Atomsk 生成蜂窝状多晶，并自动 wrap、导出 LAMMPS 文件。

    x_shift_factor / y_shift_factor 的意义
    ------------------------------------
    它们不是直接的长度，而是“相对于晶粒中心间距的比例系数”。

    先定义：
        dx = d
        dy = sqrt(3)/2 * d

    实际平移量为：
        x_shift = x_shift_factor * dx
        y_shift = y_shift_factor * dy

    例如：
        x_shift_factor = 0.25
        y_shift_factor = 0.50

    表示：
        seed 点整体沿 x 方向平移 0.25 * dx
        seed 点整体沿 y 方向平移 0.50 * dy

    为什么要调这个？
    ----------------
    因为主胞边界会“切到”最外圈六边形。
    如果 seed 的整体位置不合适，边界处看到的晶粒就会比较歪。
    通过调这两个系数，可以改变“边界切在什么位置”，让主胞截图更规整。

    怎么调？
    -------
    最常试这几组：
        (0.25, 0.50)
        (0.50, 0.50)
        (0.25, 0.25)
        (0.00, 0.00)

    经验：
    - 如果左右边界看着别扭，优先调 x_shift_factor
    - 如果上下边界看着别扭，优先调 y_shift_factor
    - 每次建议只改一个参数，方便比较效果
    """
    # 1) resolve atomsk path
    if atomsk_path == "auto":
        atomsk = shutil.which("atomsk")
    else:
        atomsk = atomsk_path

    if atomsk is None or not Path(atomsk).exists():
        print("Error: Atomsk not found.")
        print("Please check --atomsk or set atomsk_path correctly.")
        sys.exit(1)

    workdir = Path(workdir_name).resolve()
    workdir.mkdir(exist_ok=True)

    # 2) box from triangular lattice
    dx = d
    dy = math.sqrt(3.0) * 0.5 * d

    Lx = nx_cols * dx
    Ly = ny_rows * dy
    Lz = lz

    x_shift = x_shift_factor * dx
    y_shift = y_shift_factor * dy

    est_n = estimate_atoms_fcc(Lx, Ly, Lz, a0)
    print(f"Estimated atom count before GB processing: {est_n:.0f}")

    seed_file = workdir / "seed.xsf"
    poly_file = workdir / "poly_hex_regular.txt"
    poly_cfg = workdir / "hex_poly_regular.cfg"
    wrapped_cfg = workdir / "hex_poly_regular_wrapped.cfg"

    # 3) create seed crystal
    run([
        atomsk,
        "--create", "fcc", a0, "Cu",
        "orient", seed_orient_x, seed_orient_y, seed_orient_z,
        seed_file
    ])

    # 4) generate nodes
    nodes, dx, dy = generate_hex_nodes(
        nx_cols=nx_cols,
        ny_rows=ny_rows,
        d=d,
        lz=Lz,
        angle_pool=angle_pool,
        x_shift=x_shift,
        y_shift=y_shift
    )

    # 5) write polycrystal file
    lines = [f"box {Lx:.6f} {Ly:.6f} {Lz:.6f}", ""]
    for x, y, z, ang in nodes:
        x = x % Lx
        y = y % Ly
        lines.append(f"node {x:.6f} {y:.6f} {z:.6f}   0   0   {ang}")
    lines.append("")
    poly_file.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote: {poly_file}")
    print(f"Number of nodes: {len(nodes)}")
    print(f"Box: {Lx:.3f} x {Ly:.3f} x {Lz:.3f} Å")
    print(f"Shift factor: x = {x_shift_factor:.3f}, y = {y_shift_factor:.3f}")
    print(f"Actual shift: x = {x_shift:.3f} Å, y = {y_shift:.3f} Å")

    # 6) build polycrystal
    run([
        atomsk,
        "--polycrystal",
        seed_file,
        poly_file,
        poly_cfg
    ])

    # 7) wrap atoms back
    run([
        atomsk,
        poly_cfg,
        "-wrap",
        wrapped_cfg
    ])

    # 8) export to LAMMPS
    run([
        atomsk,
        wrapped_cfg,
        "lammps"
    ])

    print("\nDone.")
    print("Generated files:")
    for p in sorted(workdir.iterdir()):
        print(" -", p.name)

    print("\nRecommended usage:")
    print("1. OVITO: Wrap at periodic boundaries")
    print("2. OVITO: Replicate in x/y to check periodic tiling")
    print("3. LAMMPS: boundary p p p")
    print(f"Box: {Lx:.3f} x {Ly:.3f} x {Lz:.3f} Å")


def parse_angle_pool(text):
    """
    把命令行里的字符串解析成角度列表。
    例如：
        '0,30,60,90' -> [0.0, 30.0, 60.0, 90.0]
    """
    parts = [x.strip() for x in text.split(",") if x.strip()]
    if not parts:
        raise argparse.ArgumentTypeError("angle pool cannot be empty")
    try:
        return [float(x) for x in parts]
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"invalid angle list: {text!r}. Use format like 0,30,60,90"
        ) from e


def build_parser(defaults):
    description = """
Generate a hex-like Voronoi polycrystal Cu model with Atomsk.

About x_shift_factor and y_shift_factor:
  They control how much all seed points are shifted before Voronoi construction.

  First define:
      dx = d
      dy = sqrt(3)/2 * d

  Then the actual shift is:
      x_shift = x_shift_factor * dx
      y_shift = y_shift_factor * dy

Why tune them?
  Because the periodic box cuts through the outer hexagons.
  Different shifts change where the box boundary cuts the pattern.
  This can make boundary grains look more or less regular.

Typical values to try:
  --shiftx 0.25 --shifty 0.50
  --shiftx 0.50 --shifty 0.50
  --shiftx 0.25 --shifty 0.25
  --shiftx 0.00 --shifty 0.00
"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=description
    )

    parser.add_argument(
        "--atomsk",
        default=defaults["atomsk_path"],
        help=f"Path to Atomsk executable, or 'auto' to search in PATH. Default: {defaults['atomsk_path']}"
    )
    parser.add_argument(
        "--workdir",
        default=defaults["workdir_name"],
        help=f"Output directory. Default: {defaults['workdir_name']}"
    )
    parser.add_argument(
        "--a0",
        type=float,
        default=defaults["a0"],
        help=f"Lattice constant in Angstrom. Default: {defaults['a0']}"
    )
    parser.add_argument(
        "--nx",
        type=int,
        default=defaults["nx_cols"],
        help=f"Number of seed columns. Default: {defaults['nx_cols']}"
    )
    parser.add_argument(
        "--ny",
        type=int,
        default=defaults["ny_rows"],
        help=f"Number of seed rows. Default: {defaults['ny_rows']}"
    )
    parser.add_argument(
        "--d",
        type=float,
        default=defaults["d"],
        help=f"Horizontal spacing between neighboring seeds in Angstrom. Default: {defaults['d']}"
    )
    parser.add_argument(
        "--lz",
        type=float,
        default=defaults["lz"],
        help=f"Box thickness along z in Angstrom. Default: {defaults['lz']}"
    )
    parser.add_argument(
        "--angles",
        type=parse_angle_pool,
        default=defaults["angle_pool"],
        help=f"Comma-separated grain angle list. Example: 0,30,60,90  Default: {','.join(map(str, defaults['angle_pool']))}"
    )
    parser.add_argument(
        "--shiftx",
        type=float,
        default=defaults["x_shift_factor"],
        help=(
            "Seed shift factor along x.\n"
            "Actual x shift = shiftx * d.\n"
            "Used to tune how the boundary cuts hexagons.\n"
            f"Default: {defaults['x_shift_factor']}"
        )
    )
    parser.add_argument(
        "--shifty",
        type=float,
        default=defaults["y_shift_factor"],
        help=(
            "Seed shift factor along y.\n"
            "Actual y shift = shifty * (sqrt(3)/2*d).\n"
            "Used to tune how the boundary cuts hexagons.\n"
            f"Default: {defaults['y_shift_factor']}"
        )
    )
    parser.add_argument(
        "--ox",
        default=defaults["seed_orient_x"],
        help=f"Seed crystal orientation x. Must be orthogonal to oy and oz. Default: {defaults['seed_orient_x']}"
    )
    parser.add_argument(
        "--oy",
        default=defaults["seed_orient_y"],
        help=f"Seed crystal orientation y. Must be orthogonal to ox and oz. Default: {defaults['seed_orient_y']}"
    )
    parser.add_argument(
        "--oz",
        default=defaults["seed_orient_z"],
        help=f"Seed crystal orientation z. Must be orthogonal to ox and oy. Default: {defaults['seed_orient_z']}"
    )


    return parser


if __name__ == "__main__":
    # ============================================================
    # 默认参数区：直接改这里就行
    # 运行 python script.py 时，会使用这里的默认值
    # 如果命令行传了参数，就会覆盖这里的默认值
    # ============================================================

    # Atomsk 路径：
    # - 如果 atomsk 已加入 PATH，可以写 "auto"
    # - 如果没有，就写 atomsk.exe 的完整路径
    atomsk_path = r"D:\2025.916aftersoftware\atmosk\Atomsk\atomsk.EXE"
    # atomsk_path = "auto"

    # 输出目录名
    workdir_name = "atomsk_hex_poly_regular"

    # Cu 晶格常数（Å）
    a0 = 3.615

    # seed 点列数、行数
    # 越大 -> 晶粒数越多，盒子越大，总原子数越多
    nx_cols = 4
    ny_rows = 4

    # seed 点水平间距 d（Å）
    # 越小 -> 晶粒更小更密
    d = 40

    # z 方向厚度（Å）
    # 越大 -> 原子数越多
    # lz = 18.0
    nz_period = 7
    lz = nz_period * a0

    # 晶粒取向角列表
    # 会在所有晶粒中循环分配
    angle_pool = [0, 30, 60, 90]

    # ------------------------------------------------------------
    # 重点：边界规整度调节参数
    #
    # x_shift_factor / y_shift_factor 不是直接长度，
    # 而是“按晶粒中心间距的比例”来平移所有 seed 点。
    #
    # 定义：
    #   dx = d
    #   dy = sqrt(3)/2 * d
    #
    # 实际平移量：
    #   x_shift = x_shift_factor * dx
    #   y_shift = y_shift_factor * dy
    #
    # 例子：
    #   x_shift_factor = 0.25
    #   y_shift_factor = 0.50
    #
    # 表示：
    #   所有 seed 点整体向 x 方向平移 0.25*dx
    #   所有 seed 点整体向 y 方向平移 0.50*dy
    #
    # 为什么要调这个？
    #   主胞边界会截断最外圈六边形。
    #   改这个平移量，本质上是在改“边界切到蜂窝图案的哪里”。
    #   调得合适，边界看起来就更规整。
    #
    # 怎么调？
    #   常试这几组：
    #   1) 0.25, 0.50   -> 当前常用，通常比较均衡
    #   2) 0.50, 0.50   -> 往中心再推一点
    #   3) 0.25, 0.25   -> 上下边界切法会不同
    #   4) 0.00, 0.00   -> 不平移，作为对照
    #
    # 调参经验：
    #   - 左右边界不顺眼：先调 x_shift_factor
    #   - 上下边界不顺眼：先调 y_shift_factor
    #   - 每次只改一个，容易比较
    # ------------------------------------------------------------
    x_shift_factor = 0.25
    y_shift_factor = 0.50

    # seed 单晶取向，必须两两正交
    seed_orient_x = "[100]"
    seed_orient_y = "[010]"
    seed_orient_z = "[001]"

    defaults = {
        "atomsk_path": atomsk_path,
        "workdir_name": workdir_name,
        "a0": a0,
        "nx_cols": nx_cols,
        "ny_rows": ny_rows,
        "d": d,
        "lz": lz,
        "angle_pool": angle_pool,
        "x_shift_factor": x_shift_factor,
        "y_shift_factor": y_shift_factor,
        "seed_orient_x": seed_orient_x,
        "seed_orient_y": seed_orient_y,
        "seed_orient_z": seed_orient_z,
    }

    parser = build_parser(defaults)
    args = parser.parse_args()

    build_model(
        atomsk_path=args.atomsk,
        workdir_name=args.workdir,
        a0=args.a0,
        nx_cols=args.nx,
        ny_rows=args.ny,
        d=args.d,
        lz=args.lz,
        angle_pool=args.angles,
        x_shift_factor=args.shiftx,
        y_shift_factor=args.shifty,
        seed_orient_x=args.ox,
        seed_orient_y=args.oy,
        seed_orient_z=args.oz,
    )

