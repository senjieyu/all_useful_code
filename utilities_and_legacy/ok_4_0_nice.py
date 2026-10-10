#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import shutil
import subprocess
from pathlib import Path
import sys
import argparse


def run(cmd, cwd=None):
    """运行外部命令。"""
    cmd = [str(x) for x in cmd]
    print(">>", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def estimate_atoms_fcc(lx, ly, lz, a0):
    """
    估算 FCC 体系原子数：
        N ≈ 4 * V / a0^3
    """
    return 4.0 * lx * ly * lz / (a0 ** 3)


def parse_angle_pool(text):
    """
    解析角度列表。
    例如：
        '0,30,60,90' -> [0.0, 30.0, 60.0, 90.0]
    """
    parts = [x.strip() for x in text.split(",") if x.strip()]
    if not parts:
        raise argparse.ArgumentTypeError("角度列表不能为空")
    try:
        return [float(x) for x in parts]
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"角度列表格式错误：{text!r}，应写成 0,30,60,90"
        ) from e


def parse_float_list(text):
    """
    解析浮点数列表。
    例如：
        '0.98,1.00,1.02' -> [0.98, 1.0, 1.02]
    """
    parts = [x.strip() for x in text.split(",") if x.strip()]
    if not parts:
        raise argparse.ArgumentTypeError("浮点数列表不能为空")
    try:
        return [float(x) for x in parts]
    except ValueError as e:
        raise argparse.ArgumentTypeError(
            f"浮点数列表格式错误：{text!r}，应写成 0.98,1.00,1.02"
        ) from e


def generate_hex_nodes(nx_cols, ny_rows, d, lz, angle_pool, x_shift=0.0, y_shift=0.0):
    """
    生成三角晶格 seed 点。
    对这些点做 Voronoi 剖分后，会得到近似六边形/蜂窝状晶粒。

    参数
    ----
    nx_cols, ny_rows : int
        seed 点列数、行数
    d : float
        相邻 seed 的水平间距（Å）
    lz : float
        盒子 z 长度（Å）
    angle_pool : list[float]
        晶粒取向角列表，会循环分配
    x_shift, y_shift : float
        对所有 seed 点施加整体平移（单位 Å）

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


def rescale_nodes_about_box_center(nodes, lx, ly, scale):
    """
    以盒子中心为参考，对 seed 点在 x/y 平面内做整体缩放。

    含义
    ----
    scale > 1.0
        seed 点整体远离盒子中心
        晶粒之间更松、更分开

    scale < 1.0
        seed 点整体靠近盒子中心
        晶粒之间更紧、更挤

    scale = 1.0
        不变

    说明
    ----
    这不是修改名义上的 d，
    而是在 d 固定时，用来调节 seed 的“实际相对位置”。

    这是一个很实用的“晶界筛选参数”：
    - 太挤时可以稍微放松
    - 太稀时可以稍微收紧
    """
    cx = 0.5 * lx
    cy = 0.5 * ly

    new_nodes = []
    for x, y, z, ang in nodes:
        x_new = cx + scale * (x - cx)
        y_new = cy + scale * (y - cy)

        # 包回主胞
        x_new = x_new % lx
        y_new = y_new % ly

        new_nodes.append((x_new, y_new, z, ang))

    return new_nodes


def format_scale_tag(scale):
    """
    把 packing_scale 转成适合目录名的字符串。
    例如：
        0.98 -> pscale_0p98
        1.00 -> pscale_1p00
    """
    return f"pscale_{scale:.2f}".replace(".", "p")


def build_one_model(
    atomsk_path,
    root_workdir,
    a0,
    nx_cols,
    ny_rows,
    d,
    lz,
    angle_pool,
    x_shift_factor,
    y_shift_factor,
    packing_scale,
    seed_orient_x,
    seed_orient_y,
    seed_orient_z,
):
    """
    用单个 packing_scale 生成一个模型。
    每个 packing_scale 都写入独立子目录，避免文件覆盖。
    """
    # 解析 Atomsk 路径
    if atomsk_path == "auto":
        atomsk = shutil.which("atomsk")
    else:
        atomsk = atomsk_path

    if atomsk is None or not Path(atomsk).exists():
        print("错误：找不到 Atomsk。")
        print("请检查 atomsk_path，或用 --atomsk 指定。")
        sys.exit(1)

    # 三角晶格步长
    dx = d
    dy = math.sqrt(3.0) * 0.5 * d

    # 盒子尺寸
    lx = nx_cols * dx
    ly = ny_rows * dy
    lz_box = lz

    # 用于调节边界切割位置的整体平移
    x_shift = x_shift_factor * dx
    y_shift = y_shift_factor * dy

    # 当前 packing_scale 对应的输出目录
    subdir = root_workdir / format_scale_tag(packing_scale)
    subdir.mkdir(parents=True, exist_ok=True)

    est_n = estimate_atoms_fcc(lx, ly, lz_box, a0)

    print("\n" + "=" * 72)
    print(f"开始生成：packing_scale = {packing_scale:.4f}")
    print(f"输出目录：{subdir}")
    print(f"估算原子数：{est_n:.0f}")
    print(f"盒子尺寸：{lx:.3f} x {ly:.3f} x {lz_box:.3f} Å")
    print(f"shift 系数：x = {x_shift_factor:.3f}, y = {y_shift_factor:.3f}")
    print(f"实际 shift：x = {x_shift:.3f} Å, y = {y_shift:.3f} Å")
    print("=" * 72)

    # 为避免 Atomsk 路径歧义，这里统一使用短文件名 + cwd
    seed_name = "seed.xsf"
    poly_name = "poly_hex_regular.txt"
    cfg_name = "poly.cfg"
    wrapped_name = "wrapped.cfg"

    poly_file = subdir / poly_name

    # 1) 生成 seed 单晶
    run([
        atomsk,
        "--create", "fcc", a0, "Cu",
        "orient", seed_orient_x, seed_orient_y, seed_orient_z,
        seed_name
    ], cwd=subdir)

    # 2) 生成三角晶格 seed 点
    nodes, dx, dy = generate_hex_nodes(
        nx_cols=nx_cols,
        ny_rows=ny_rows,
        d=d,
        lz=lz_box,
        angle_pool=angle_pool,
        x_shift=x_shift,
        y_shift=y_shift
    )

    # 3) 用 packing_scale 调整 seed 相对位置
    nodes = rescale_nodes_about_box_center(nodes, lx, ly, packing_scale)

    # 4) 写入 polycrystal 参数文件
    lines = [f"box {lx:.6f} {ly:.6f} {lz_box:.6f}", ""]
    for x, y, z, ang in nodes:
        x = x % lx
        y = y % ly
        lines.append(f"node {x:.6f} {y:.6f} {z:.6f}   0   0   {ang}")
    lines.append("")

    poly_file.write_text("\n".join(lines), encoding="utf-8")
    print(f"已写入：{poly_file}")

    # 5) 用 Atomsk 构建多晶
    run([
        atomsk,
        "--polycrystal",
        seed_name,
        poly_name,
        cfg_name
    ], cwd=subdir)

    # 6) 把越界原子包回主胞
    run([
        atomsk,
        cfg_name,
        "-wrap",
        wrapped_name
    ], cwd=subdir)

    # 7) 导出为 LAMMPS data
    run([
        atomsk,
        wrapped_name,
        "lammps"
    ], cwd=subdir)

    print(f"完成：packing_scale = {packing_scale:.4f}")
    print("生成文件：")
    for p in sorted(subdir.iterdir()):
        print(" -", p.name)


def build_parser(defaults):
    """构建命令行参数解析器。"""
    description = """
用 Atomsk 生成六边形/蜂窝风格的 Voronoi 多晶 Cu 模型。

核心参数说明
------------
1) d
   名义上的 seed 间距，主要控制晶粒大小。

2) shiftx / shifty
   用来调“边界怎么切六边形”的。
   它们会在 Voronoi 前整体平移所有 seed 点：

       x_shift = shiftx * d
       y_shift = shifty * (sqrt(3)/2 * d)

3) pscales
   packing_scale 列表。
   会一次生成多组模型，用于筛选更合理的晶界。

   含义：
       > 1.0  -> 更松、更分开
       < 1.0  -> 更紧、更靠近
       = 1.0  -> 不变

   推荐先试：
       0.98,0.99,1.00,1.01,1.02

建议流程
--------
固定 d 不变，扫描 pscales：
    0.98, 0.99, 1.00, 1.01, 1.02

然后对每个模型统一做：
    1. delete_atoms overlap
    2. minimize
    3. 比较：
       - 删除了多少原子
       - 每原子势能
       - 晶粒内部是不是主要保持 FCC
"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=description
    )

    parser.add_argument(
        "--atomsk",
        default=defaults["atomsk_path"],
        help=f"Atomsk 可执行文件路径，或写 auto 从 PATH 中寻找。默认：{defaults['atomsk_path']}"
    )
    parser.add_argument(
        "--workdir",
        default=defaults["workdir_name"],
        help=f"输出根目录。默认：{defaults['workdir_name']}"
    )
    parser.add_argument(
        "--a0",
        type=float,
        default=defaults["a0"],
        help=f"晶格常数（Å）。默认：{defaults['a0']}"
    )
    parser.add_argument(
        "--nx",
        type=int,
        default=defaults["nx_cols"],
        help=f"seed 点列数。默认：{defaults['nx_cols']}"
    )
    parser.add_argument(
        "--ny",
        type=int,
        default=defaults["ny_rows"],
        help=f"seed 点行数。默认：{defaults['ny_rows']}"
    )
    parser.add_argument(
        "--d",
        type=float,
        default=defaults["d"],
        help=f"名义上的相邻 seed 水平间距（Å）。默认：{defaults['d']}"
    )
    parser.add_argument(
        "--lz",
        type=float,
        default=defaults["lz"],
        help=f"盒子 z 方向厚度（Å）。默认：{defaults['lz']}"
    )
    parser.add_argument(
        "--angles",
        type=parse_angle_pool,
        default=defaults["angle_pool"],
        help=f"晶粒取向角列表，例如 0,30,60,90。默认：{','.join(map(str, defaults['angle_pool']))}"
    )
    parser.add_argument(
        "--shiftx",
        type=float,
        default=defaults["x_shift_factor"],
        help=(
            "seed 在 x 方向的整体平移系数。\n"
            "实际 x_shift = shiftx * d\n"
            f"默认：{defaults['x_shift_factor']}"
        )
    )
    parser.add_argument(
        "--shifty",
        type=float,
        default=defaults["y_shift_factor"],
        help=(
            "seed 在 y 方向的整体平移系数。\n"
            "实际 y_shift = shifty * (sqrt(3)/2*d)\n"
            f"默认：{defaults['y_shift_factor']}"
        )
    )
    parser.add_argument(
        "--pscales",
        type=parse_float_list,
        default=defaults["packing_scales"],
        help=(
            "packing_scale 列表，会一次生成多组模型。\n"
            "例如：0.98,1.00,1.02\n"
            ">1 更松，<1 更紧\n"
            f"默认：{','.join(f'{x:.2f}' for x in defaults['packing_scales'])}"
        )
    )
    parser.add_argument(
        "--ox",
        default=defaults["seed_orient_x"],
        help=f"seed 单晶 x 方向晶向。默认：{defaults['seed_orient_x']}"
    )
    parser.add_argument(
        "--oy",
        default=defaults["seed_orient_y"],
        help=f"seed 单晶 y 方向晶向。默认：{defaults['seed_orient_y']}"
    )
    parser.add_argument(
        "--oz",
        default=defaults["seed_orient_z"],
        help=f"seed 单晶 z 方向晶向。默认：{defaults['seed_orient_z']}"
    )

    return parser


if __name__ == "__main__":
    # ============================================================
    # 默认参数区
    # 平时你主要改下面这些参数即可
    # 命令行参数如果提供，会覆盖这里的默认值
    # ============================================================

    # ------------------------------------------------------------
    # Atomsk 路径
    # 如果 atomsk 已加入 PATH，可以写 "auto"
    # 否则写 atomsk.exe 的完整路径
    # ------------------------------------------------------------
    atomsk_path = r"D:\2025.916aftersoftware\atmosk\Atomsk\atomsk.EXE"
    # atomsk_path = "auto"

    # ------------------------------------------------------------
    # 输出根目录
    # 每个 packing_scale 会在这个目录下生成一个子目录
    # ------------------------------------------------------------
    workdir_name = "atomsk_hex_poly_batch"

    # ------------------------------------------------------------
    # Cu 晶格常数（Å）
    # ------------------------------------------------------------
    a0 = 3.615

    # ------------------------------------------------------------
    # seed 点列数、行数
    # 越大 -> 晶粒数越多，体系尺寸越大
    # ------------------------------------------------------------
    nx_cols = 4
    ny_rows = 4

    # ------------------------------------------------------------
    # 名义上的 seed 间距 d（Å）
    # 越大 -> 晶粒整体越大
    # ------------------------------------------------------------
    d = 20.0

    # ------------------------------------------------------------
    # z 方向厚度（Å）
    # 对 [100]/[010]/[001] 取向，建议尽量取接近 a0 的整数倍
    # 例如：
    #   5*a0
    #   6*a0
    #   7*a0
    # ------------------------------------------------------------
    lz = 5 * a0

    # ------------------------------------------------------------
    # 晶粒取向角列表（单位：度）
    # 会循环分配给各个晶粒
    # ------------------------------------------------------------
    angle_pool = [0, 30, 60, 90]

    # ------------------------------------------------------------
    # 边界切割位置调节参数
    #
    # 它们不会直接改变晶粒大小，
    # 而是改变“主胞边界切到蜂窝图案的哪里”。
    #
    # 定义：
    #   dx = d
    #   dy = sqrt(3)/2 * d
    #
    # 实际平移：
    #   x_shift = x_shift_factor * dx
    #   y_shift = y_shift_factor * dy
    #
    # 常试值：
    #   (0.25, 0.50)
    #   (0.50, 0.50)
    #   (0.25, 0.25)
    #   (0.00, 0.00)
    # ------------------------------------------------------------
    x_shift_factor = 0.25
    y_shift_factor = 0.50

    # ------------------------------------------------------------
    # packing_scale 列表：一次生成一组模型
    #
    # 含义：
    #   > 1.0
    #       seed 点整体远离盒子中心
    #       晶粒之间更松、更分开
    #
    #   < 1.0
    #       seed 点整体靠近盒子中心
    #       晶粒之间更紧、更靠近
    #
    #   = 1.0
    #       不变
    #
    # 这个参数非常适合在 d 固定时筛选“更合理的晶界”。
    #
    # 推荐第一轮就扫这一组：
    #   [0.98, 0.99, 1.00, 1.01, 1.02]
    #
    # 如果你发现：
    # - 晶界太挤、重叠偏多 -> 试更大的值
    # - 晶界太稀、边界开口偏大 -> 试更小的值
    #
    # 最后再配合：
    #   delete_atoms overlap + minimize
    # 去比较哪一组最合理。
    # ------------------------------------------------------------
    packing_scales = [0.98, 0.99, 1.00, 1.01, 1.02]

    # ------------------------------------------------------------
    # seed 单晶取向
    # 这三个方向必须两两正交
    # 稳妥默认值：
    #   [100], [010], [001]
    # ------------------------------------------------------------
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
        "packing_scales": packing_scales,
        "seed_orient_x": seed_orient_x,
        "seed_orient_y": seed_orient_y,
        "seed_orient_z": seed_orient_z,
    }

    parser = build_parser(defaults)
    args = parser.parse_args()

    root_workdir = Path(args.workdir).resolve()
    root_workdir.mkdir(parents=True, exist_ok=True)

    print("将要生成的 packing_scale 列表：")
    for val in args.pscales:
        print(f" - {val:.4f}")

    for pscale in args.pscales:
        build_one_model(
            atomsk_path=args.atomsk,
            root_workdir=root_workdir,
            a0=args.a0,
            nx_cols=args.nx,
            ny_rows=args.ny,
            d=args.d,
            lz=args.lz,
            angle_pool=args.angles,
            x_shift_factor=args.shiftx,
            y_shift_factor=args.shifty,
            packing_scale=pscale,
            seed_orient_x=args.ox,
            seed_orient_y=args.oy,
            seed_orient_z=args.oz,
        )

    print("\n全部模型生成完成。")
    print("建议后续对每个子目录中的 wrapped.lmp 统一做：")
    print("1. delete_atoms overlap")
    print("2. minimize")
    print("3. 比较删除原子数、每原子势能、CNA/PTM 结果")
