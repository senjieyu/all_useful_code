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

    注意：
    当前脚本默认使用固定蜂窝角度分布：
        偶数行：0, 30, 0, 30, ...
        奇数行：90, 60, 90, 60, ...

    因此 --angles 主要保留为兼容参数。
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


def get_honeycomb_angle(row, col):
    """
    固定蜂窝角度分布。

    目标图案：

        0°    30°    0°    30°
          90°    60°    90°    60°
        0°    30°    0°    30°
          90°    60°    90°    60°

    规则：
        偶数行：
            偶数列 -> 0°
            奇数列 -> 30°

        奇数行：
            偶数列 -> 90°
            奇数列 -> 60°
    """
    if row % 2 == 0:
        if col % 2 == 0:
            return 0.0
        else:
            return 30.0
    else:
        if col % 2 == 0:
            return 90.0
        else:
            return 60.0


def generate_hex_nodes(nx_cols, ny_rows, d, lz, angle_pool=None, x_shift=0.0, y_shift=0.0):
    """
    生成三角晶格 seed 点。
    对这些点做 Voronoi 剖分后，会得到近似六边形/蜂窝状晶粒。

    当前版本使用固定角度分布：

        0°    30°    0°    30°
          90°    60°    90°    60°
        0°    30°    0°    30°
          90°    60°    90°    60°

    这样可以避免相邻晶粒出现相同取向。

    参数
    ----
    nx_cols, ny_rows : int
        seed 点列数、行数
    d : float
        相邻 seed 的水平间距，单位 Å
    lz : float
        盒子 z 长度，单位 Å
    angle_pool : list[float] or None
        保留参数，当前默认不使用
    x_shift, y_shift : float
        对所有 seed 点施加整体平移，单位 Å

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

            ang = get_honeycomb_angle(row, col)

            nodes.append((x, y, 0.5 * lz, ang))

    return nodes, dx, dy


def rescale_nodes_about_box_center(nodes, lx, ly, scale):
    """
    以盒子中心为参考，对 seed 点在 x/y 平面内做整体缩放。

    scale > 1.0:
        seed 点整体远离盒子中心，晶粒之间更松、更分开

    scale < 1.0:
        seed 点整体靠近盒子中心，晶粒之间更紧、更挤

    scale = 1.0:
        不变
    """
    cx = 0.5 * lx
    cy = 0.5 * ly

    new_nodes = []
    for x, y, z, ang in nodes:
        x_new = cx + scale * (x - cx)
        y_new = cy + scale * (y - cy)

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


def write_nanotwin_seed_xsf(
    output_path,
    a0,
    element,
    twin_spacing,
    twin_periods,
    nx_inplane,
    ny_inplane,
):
    """
    生成一个周期性纳米孪晶 FCC seed，并写成 XSF 文件。

    结构设定
    --------
    取向固定为：
        x ≈ [11-2]
        y ≈ [1-10]
        z ≈ [111]

    z 方向为 FCC {111} 层法向。
    通过改变 close-packed 层的 ABC 堆垛方向来形成孪晶：

        正常 FCC: A B C A B C ...
        孪晶反向: C B A C B A ...

    twin_spacing 控制相邻孪晶界之间的名义距离，单位 Å。

    注意
    ----
    由于 FCC {111} 层间距是：
        d111 = a0 / sqrt(3)

    实际孪晶间距会被自动取整到：
        layers_per_twin * d111

    所以实际值会接近但不一定完全等于输入值。
    """
    if twin_spacing <= 0:
        raise ValueError("twin_spacing 必须大于 0")

    if twin_periods < 2:
        raise ValueError("twin_periods 至少为 2")

    if nx_inplane < 1 or ny_inplane < 1:
        raise ValueError("nx_inplane 和 ny_inplane 必须为正整数")

    d111 = a0 / math.sqrt(3.0)
    a2d = a0 / math.sqrt(2.0)

    layers_per_twin = max(2, int(round(twin_spacing / d111)))
    actual_spacing = layers_per_twin * d111

    n_layers = layers_per_twin * twin_periods

    print("纳米孪晶 seed 参数：")
    print(f"  输入 twin_spacing     = {twin_spacing:.6f} Å")
    print(f"  FCC d111              = {d111:.6f} Å")
    print(f"  layers_per_twin       = {layers_per_twin}")
    print(f"  实际 twin spacing     = {actual_spacing:.6f} Å")
    print(f"  twin_periods          = {twin_periods}")
    print(f"  total {111} layers    = {n_layers}")
    print(f"  in-plane repeats      = {nx_inplane} x {ny_inplane}")

    # close-packed (111) 面内二维基矢
    # v1 和 v2 夹角 60°
    v1 = (a2d, 0.0, 0.0)
    v2 = (0.5 * a2d, 0.5 * math.sqrt(3.0) * a2d, 0.0)
    v3 = (0.0, 0.0, n_layers * d111)

    # A/B/C 堆垛在二维原胞内的偏移
    stack_offsets = [
        (0.0, 0.0, 0.0),
        ((v1[0] + v2[0]) / 3.0, (v1[1] + v2[1]) / 3.0, 0.0),
        (2.0 * (v1[0] + v2[0]) / 3.0, 2.0 * (v1[1] + v2[1]) / 3.0, 0.0),
    ]

    atoms = []

    stack_index = 0
    direction = 1

    for layer in range(n_layers):
        ox, oy, _ = stack_offsets[stack_index % 3]
        z = layer * d111

        for i in range(nx_inplane):
            for j in range(ny_inplane):
                x = i * v1[0] + j * v2[0] + ox
                y = i * v1[1] + j * v2[1] + oy
                atoms.append((element, x, y, z))

        # 每隔 layers_per_twin 层反转 ABC 堆垛方向，形成孪晶界
        if (layer + 1) % layers_per_twin == 0:
            direction *= -1

        stack_index += direction

    cell_a = (nx_inplane * v1[0], nx_inplane * v1[1], nx_inplane * v1[2])
    cell_b = (ny_inplane * v2[0], ny_inplane * v2[1], ny_inplane * v2[2])
    cell_c = v3

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("CRYSTAL\n")
        f.write("PRIMVEC\n")
        f.write(f"{cell_a[0]:.12f} {cell_a[1]:.12f} {cell_a[2]:.12f}\n")
        f.write(f"{cell_b[0]:.12f} {cell_b[1]:.12f} {cell_b[2]:.12f}\n")
        f.write(f"{cell_c[0]:.12f} {cell_c[1]:.12f} {cell_c[2]:.12f}\n")
        f.write("PRIMCOORD\n")
        f.write(f"{len(atoms)} 1\n")
        for sym, x, y, z in atoms:
            f.write(f"{sym} {x:.12f} {y:.12f} {z:.12f}\n")

    print(f"已写入纳米孪晶 seed：{output_path}")


def build_seed(
    atomsk,
    subdir,
    a0,
    element,
    seed_name,
    use_twin,
    twin_spacing,
    twin_periods,
    twin_nx,
    twin_ny,
    seed_orient_x,
    seed_orient_y,
    seed_orient_z,
):
    """
    生成 seed.xsf。

    use_twin = False:
        用 Atomsk 生成纯 FCC seed。

    use_twin = True:
        用 Python 直接生成周期性纳米孪晶 seed。
    """
    seed_path = subdir / seed_name

    if use_twin:
        write_nanotwin_seed_xsf(
            output_path=seed_path,
            a0=a0,
            element=element,
            twin_spacing=twin_spacing,
            twin_periods=twin_periods,
            nx_inplane=twin_nx,
            ny_inplane=twin_ny,
        )
    else:
        run([
            atomsk,
            "--create", "fcc", a0, element,
            "orient", seed_orient_x, seed_orient_y, seed_orient_z,
            seed_name
        ], cwd=subdir)


def write_poly_file(poly_file, lx, ly, lz_box, nodes):
    """
    写入 Atomsk --polycrystal 所需的 polycrystal 参数文件。
    """
    lines = [f"box {lx:.6f} {ly:.6f} {lz_box:.6f}", ""]

    lines.append("# Honeycomb orientation pattern:")
    lines.append("# even rows: 0 30 0 30 ...")
    lines.append("# odd  rows: 90 60 90 60 ...")
    lines.append("")

    for x, y, z, ang in nodes:
        x = x % lx
        y = y % ly

        # Atomsk node:
        # node x y z  rotX rotY rotZ
        #
        # 这里绕 z 方向旋转晶粒。
        lines.append(f"node {x:.6f} {y:.6f} {z:.6f}   0   0   {ang:.6f}")

    lines.append("")
    poly_file.write_text("\n".join(lines), encoding="utf-8")
    print(f"已写入：{poly_file}")


def print_angle_map(nx_cols, ny_rows):
    """
    打印角度分布，方便检查是否符合目标图案。
    """
    print("\n蜂窝角度分布检查：")
    for row in range(ny_rows):
        indent = "    " if row % 2 == 1 else ""
        vals = []
        for col in range(nx_cols):
            vals.append(f"{get_honeycomb_angle(row, col):.0f}°")
        print(indent + "    ".join(vals))
    print("")


def check_even_grid(nx_cols, ny_rows):
    """
    提醒用户最好使用偶数行列。
    """
    if nx_cols % 2 != 0 or ny_rows % 2 != 0:
        print("\n警告：")
        print("  当前 nx 或 ny 不是偶数。")
        print("  如果使用周期边界，奇数行/列可能导致盒子边界相接处出现相同角度相邻。")
        print("  建议使用 nx = 4, 6, 8 ...；ny = 4, 6, 8 ...\n")


def build_one_model(
    atomsk_path,
    root_workdir,
    element,
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
    use_twin,
    twin_spacing,
    twin_periods,
    twin_nx,
    twin_ny,
):
    """
    用单个 packing_scale 生成一个模型。
    每个 packing_scale 都写入独立子目录，避免文件覆盖。
    """
    if atomsk_path == "auto":
        atomsk = shutil.which("atomsk")
    else:
        atomsk = atomsk_path

    if atomsk is None or not Path(atomsk).exists():
        print("错误：找不到 Atomsk。")
        print("请检查 atomsk_path，或用 --atomsk 指定。")
        sys.exit(1)

    check_even_grid(nx_cols, ny_rows)

    dx = d
    dy = math.sqrt(3.0) * 0.5 * d

    lx = nx_cols * dx
    ly = ny_rows * dy
    lz_box = lz

    x_shift = x_shift_factor * dx
    y_shift = y_shift_factor * dy

    model_tag = format_scale_tag(packing_scale)
    if use_twin:
        twin_tag = f"twin_{twin_spacing:.2f}A".replace(".", "p")
        subdir = root_workdir / twin_tag / model_tag
    else:
        subdir = root_workdir / model_tag

    subdir.mkdir(parents=True, exist_ok=True)

    est_n = estimate_atoms_fcc(lx, ly, lz_box, a0)

    print("\n" + "=" * 72)
    print(f"开始生成：packing_scale = {packing_scale:.4f}")
    print(f"输出目录：{subdir}")
    print(f"元素：{element}")
    print(f"估算原子数：{est_n:.0f}")
    print(f"盒子尺寸：{lx:.3f} x {ly:.3f} x {lz_box:.3f} Å")
    print(f"seed 网格：nx = {nx_cols}, ny = {ny_rows}")
    print(f"shift 系数：x = {x_shift_factor:.3f}, y = {y_shift_factor:.3f}")
    print(f"实际 shift：x = {x_shift:.3f} Å, y = {y_shift:.3f} Å")
    print(f"是否使用纳米孪晶 seed：{use_twin}")
    if use_twin:
        print(f"目标孪晶间距：{twin_spacing:.3f} Å")
    print("=" * 72)

    print_angle_map(nx_cols, ny_rows)

    seed_name = "seed.xsf"
    poly_name = "poly_hex_regular.txt"
    cfg_name = "poly.cfg"
    wrapped_name = "wrapped.cfg"

    poly_file = subdir / poly_name

    # 1) 生成 seed
    build_seed(
        atomsk=atomsk,
        subdir=subdir,
        a0=a0,
        element=element,
        seed_name=seed_name,
        use_twin=use_twin,
        twin_spacing=twin_spacing,
        twin_periods=twin_periods,
        twin_nx=twin_nx,
        twin_ny=twin_ny,
        seed_orient_x=seed_orient_x,
        seed_orient_y=seed_orient_y,
        seed_orient_z=seed_orient_z,
    )

    # 2) 生成三角晶格 Voronoi seed 点
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

    # 4) 写入 Atomsk polycrystal 参数文件
    write_poly_file(
        poly_file=poly_file,
        lx=lx,
        ly=ly,
        lz_box=lz_box,
        nodes=nodes,
    )

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
用 Atomsk 生成六边形/蜂窝风格的 Voronoi 多晶模型。

当前角度分布固定为：

    0°    30°    0°    30°
      90°    60°    90°    60°
    0°    30°    0°    30°
      90°    60°    90°    60°

这样可以避免相邻晶粒出现 0° 和 0° 相邻的问题。

新增功能
--------
可以用 --twin 开启纳米孪晶 seed。

开启后：
    每个 Voronoi 晶粒内部不再是纯 FCC，
    而是周期性 {111} 孪晶层片结构。

关键参数
--------
1) --d
   名义上的 seed 间距，主要控制晶粒大小。

2) --twin-spacing
   纳米孪晶层片间距，单位 Å。

3) --pscales
   packing_scale 列表。
   会一次生成多组模型，用于筛选更合理的晶界。

建议
----
为了避免周期边界处出现同角度相邻，建议：

    --nx 4 --ny 4
    --nx 6 --ny 4
    --nx 6 --ny 6
    --nx 8 --ny 6

尽量不要一开始使用奇数行列，例如 5 x 5。
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
        "--element",
        default=defaults["element"],
        help=f"元素符号。默认：{defaults['element']}"
    )
    parser.add_argument(
        "--a0",
        type=float,
        default=defaults["a0"],
        help=f"晶格常数，单位 Å。默认：{defaults['a0']}"
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
        help=f"名义上的相邻 seed 水平间距，单位 Å。默认：{defaults['d']}"
    )
    parser.add_argument(
        "--lz",
        type=float,
        default=defaults["lz"],
        help=f"盒子 z 方向厚度，单位 Å。默认：{defaults['lz']}"
    )
    parser.add_argument(
        "--angles",
        type=parse_angle_pool,
        default=defaults["angle_pool"],
        help=(
            "保留兼容参数。当前默认使用固定蜂窝角度分布，"
            "不是直接按这个列表循环分配。"
        )
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

    # 纯 FCC seed 的取向参数
    parser.add_argument(
        "--ox",
        default=defaults["seed_orient_x"],
        help=f"纯 FCC seed 单晶 x 方向晶向。默认：{defaults['seed_orient_x']}"
    )
    parser.add_argument(
        "--oy",
        default=defaults["seed_orient_y"],
        help=f"纯 FCC seed 单晶 y 方向晶向。默认：{defaults['seed_orient_y']}"
    )
    parser.add_argument(
        "--oz",
        default=defaults["seed_orient_z"],
        help=f"纯 FCC seed 单晶 z 方向晶向。默认：{defaults['seed_orient_z']}"
    )

    # 纳米孪晶参数
    parser.add_argument(
        "--twin",
        action="store_true",
        default=defaults["use_twin"],
        help="开启纳米孪晶 seed。开启后晶粒内部为周期性孪晶结构。"
    )
    parser.add_argument(
        "--twin-spacing",
        type=float,
        default=defaults["twin_spacing"],
        help=f"目标孪晶间距，单位 Å。默认：{defaults['twin_spacing']}"
    )
    parser.add_argument(
        "--twin-periods",
        type=int,
        default=defaults["twin_periods"],
        help=f"seed 内包含的孪晶周期数。默认：{defaults['twin_periods']}"
    )
    parser.add_argument(
        "--twin-nx",
        type=int,
        default=defaults["twin_nx"],
        help=f"孪晶 seed 面内 x 重复数。默认：{defaults['twin_nx']}"
    )
    parser.add_argument(
        "--twin-ny",
        type=int,
        default=defaults["twin_ny"],
        help=f"孪晶 seed 面内 y 重复数。默认：{defaults['twin_ny']}"
    )

    return parser


if __name__ == "__main__":
    # ============================================================
    # 默认参数区
    # 平时主要改这里
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
    # ------------------------------------------------------------
    workdir_name = "atomsk_hex_nanotwin_poly_honeycomb_angle"

    # ------------------------------------------------------------
    # 元素和晶格常数
    # Cu: a0 ≈ 3.615 Å
    # Ni: a0 ≈ 3.52 Å
    # Al: a0 ≈ 4.05 Å
    # ------------------------------------------------------------
    element = "Cu"
    a0 = 3.615

    # ------------------------------------------------------------
    # seed 点列数、行数
    #
    # 强烈建议用偶数：
    #   4 x 4
    #   6 x 4
    #   6 x 6
    #   8 x 6
    #
    # 这样周期边界处也不容易出现相同角度相邻。
    # ------------------------------------------------------------
    nx_cols = 2
    ny_rows = 2

    # ------------------------------------------------------------
    # 名义上的 seed 间距 d，单位 Å
    # 越大 -> 晶粒整体越大
    # ------------------------------------------------------------
    d = 40.0

    # ------------------------------------------------------------
    # z 方向厚度，单位 Å
    #
    # 如果使用孪晶 seed，z 方向是孪晶层片法向方向。
    # lz 最好取为 twin spacing 的整数倍附近。
    # 例如 twin_spacing = 8 Å，可以用 lz = 32, 40, 48 Å。
    # ------------------------------------------------------------
    lz = 40.0

    # ------------------------------------------------------------
    # 保留参数。
    # 当前角度分布由 get_honeycomb_angle(row, col) 固定控制。
    # ------------------------------------------------------------
    angle_pool = [0, 30, 60, 90]

    # ------------------------------------------------------------
    # 边界切割位置调节参数
    # ------------------------------------------------------------
    x_shift_factor = 0.25
    y_shift_factor = 0.50

    # ------------------------------------------------------------
    # packing_scale 列表
    # 推荐先固定为 [1.00]，等角度分布和孪晶结构确认后再扫。
    # ------------------------------------------------------------
    packing_scales = [1.00]
    # packing_scales = [0.98, 0.99, 1.00, 1.01, 1.02]

    # ------------------------------------------------------------
    # 纯 FCC seed 单晶取向
    # 只有在不使用 --twin 时生效。
    # ------------------------------------------------------------
    seed_orient_x = "[100]"
    seed_orient_y = "[010]"
    seed_orient_z = "[001]"

    # ------------------------------------------------------------
    # 纳米孪晶 seed 参数
    #
    # use_twin = True:
    #     默认生成孪晶晶粒。
    #
    # use_twin = False:
    #     默认生成纯 FCC 晶粒，但仍可通过命令行 --twin 开启。
    #
    # twin_spacing:
    #     目标孪晶间距，单位 Å。
    #
    # twin_periods:
    #     seed 内孪晶周期数。
    #
    # twin_nx / twin_ny:
    #     孪晶 seed 在 {111} 面内的重复数。
    # ------------------------------------------------------------
    use_twin = True
    twin_spacing = 8.0
    twin_periods = 8
    twin_nx = 8
    twin_ny = 8

    defaults = {
        "atomsk_path": atomsk_path,
        "workdir_name": workdir_name,
        "element": element,
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
        "use_twin": use_twin,
        "twin_spacing": twin_spacing,
        "twin_periods": twin_periods,
        "twin_nx": twin_nx,
        "twin_ny": twin_ny,
    }

    parser = build_parser(defaults)
    args = parser.parse_args()

    root_workdir = Path(args.workdir).resolve()
    root_workdir.mkdir(parents=True, exist_ok=True)

    print("将要生成的 packing_scale 列表：")
    for val in args.pscales:
        print(f" - {val:.4f}")

    print("\n模型类型：")
    if args.twin:
        print(" - 纳米孪晶晶粒")
        print(f" - 目标孪晶间距：{args.twin_spacing:.4f} Å")
    else:
        print(" - 纯 FCC 晶粒")

    print("\n角度分布模式：")
    print(" - 固定蜂窝分布")
    print(" - 偶数行：0, 30, 0, 30, ...")
    print(" - 奇数行：90, 60, 90, 60, ...")

    for pscale in args.pscales:
        build_one_model(
            atomsk_path=args.atomsk,
            root_workdir=root_workdir,
            element=args.element,
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
            use_twin=args.twin,
            twin_spacing=args.twin_spacing,
            twin_periods=args.twin_periods,
            twin_nx=args.twin_nx,
            twin_ny=args.twin_ny,
        )

    print("\n全部模型生成完成。")
    print("建议后续对每个子目录中的 wrapped.lmp 统一做：")
    print("1. delete_atoms overlap")
    print("2. minimize")
    print("3. 比较删除原子数、每原子势能、CNA/PTM 结果")
