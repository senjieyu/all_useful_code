#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_eos_multi_xyz.py

====================================================
用法 (Usage)
====================================================

功能：
1) 自动读取指定文件夹中所有 *.xyz 文件（允许混合 extxyz 和普通 xyz）
2) 每个 xyz 文件代表一种方法（图例 label = 文件名去掉 .xyz）
3) 从每一帧的注释行读取能量：
      - energy=... / Energy=... / E=... / dft_energy=...
4) 读取体积（总 V）优先级：
      A) extxyz：注释行包含 Lattice="9 numbers" -> 用晶格行列式算体积
      B) 普通 xyz：注释行包含 V=... / vol=... / volume=... -> 直接读取体积
5) 计算并绘制 E/atom - V/atom
6) 对每种方法进行 BM3 (Birch-Murnaghan 3rd) 拟合，输出参数到表格：
      V0, E0, B0（并附带 Bp、B0(GPa)）

----------------------------------------------------
运行示例
----------------------------------------------------
(1) 画图 + 导出拟合表格
python plot_eos_multi_xyz.py ./xyz_data

(2) 指定输出图片名和表格名
python plot_eos_multi_xyz.py ./xyz_data --out eos.png --table eos_fit.csv

(3) 不画散点，只画拟合曲线
python plot_eos_multi_xyz.py ./xyz_data --no-points

(4) 只画散点，不拟合
python plot_eos_multi_xyz.py ./xyz_data --only-points

----------------------------------------------------
输出文件
----------------------------------------------------
默认输出：
- eos_multi.png
- eos_fit_table.csv

单位：
- V0：Å^3/atom
- E0：eV/atom
- B0：eV/Å^3（并换算为 GPa）
  1 eV/Å^3 ≈ 160.21766208 GPa

----------------------------------------------------
重要说明（普通 xyz 的限制）
----------------------------------------------------
如果某个 xyz 既没有 Lattice="..."，又没有在注释行提供 V=/volume=...，
则无法从该文件得到体积，会被脚本自动跳过并打印原因。

====================================================
"""

import re
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ===================== BM3 EOS =====================
def birch_murnaghan_E(V, E0, V0, B0, Bp):
    eta = (V0 / V) ** (2.0 / 3.0)
    x = eta - 1.0
    return E0 + (9.0 * V0 * B0 / 16.0) * ((x**3) * Bp + (x**2) * (6.0 - 4.0 * eta))


def fit_bm3(V, E):
    try:
        from scipy.optimize import curve_fit
    except ImportError as exc:
        raise ImportError("缺少 scipy：请先安装 `pip install scipy`") from exc

    i0 = int(np.argmin(E))
    p0 = [float(E[i0]), float(V[i0]), 0.5, 4.0]  # E0, V0, B0, Bp 初值
    popt, _ = curve_fit(birch_murnaghan_E, V, E, p0=p0, maxfev=20000)
    return popt  # E0, V0, B0, Bp


# ===================== Parsers =====================
def parse_energy(comment_line: str) -> float | None:
    """
    从注释行解析总能量 Etot
    支持：energy= / Energy= / E= / dft_energy=
    """
    m = re.search(r'(?:(?:dft_)?energy|Energy|E)\s*=\s*([-\d\.eE+]+)', comment_line)
    return None if not m else float(m.group(1))


def parse_lattice_volume(comment_line: str) -> float | None:
    """
    extxyz: 从 Lattice="9 numbers" 计算晶胞体积
    """
    m = re.search(r'Lattice="([^"]+)"', comment_line)
    if not m:
        return None
    vals = m.group(1).split()
    if len(vals) != 9:
        return None
    try:
        A = np.array([float(x) for x in vals], dtype=float).reshape(3, 3)
    except ValueError:
        return None
    return float(abs(np.linalg.det(A)))


def parse_volume_from_comment(comment_line: str) -> float | None:
    """
    普通 xyz: 尝试从注释行解析体积（总 V）
    支持示例：
      V=123.45
      vol=123.45
      volume=123.45
      Volume: 123.45
    """
    m = re.search(r'(?i)\b(?:v|vol|volume)\b\s*[:=]\s*([-\d\.eE+]+)', comment_line)
    return None if not m else float(m.group(1))


# ===================== Reader (supports extxyz + xyz) =====================
def read_xyz_or_extxyz_VE(path: Path):
    """
    读取一个文件（可能是 extxyz 或普通 xyz，多帧）:
      - Etot: 从注释行读 energy=...
      - Vtot: 优先 Lattice="..."；否则尝试 V=/volume=
    返回：V/atom, E/atom（按 V 排序）
    """
    V_list, E_list = [], []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue

            # 第一行必须是原子数
            try:
                n = int(line)
            except ValueError:
                raise ValueError(f"第一行不是原子数（可能不是 xyz 文件）。读到：{line!r}")

            comment = f.readline()
            if not comment:
                raise ValueError("缺少注释行")

            Etot = parse_energy(comment)
            if Etot is None:
                raise ValueError("注释行未找到能量字段 energy=...")

            # 体积：先尝试 extxyz 的 Lattice
            Vtot = parse_lattice_volume(comment)
            # 如果没有 Lattice，再尝试注释行 V=/volume=
            if Vtot is None:
                Vtot = parse_volume_from_comment(comment)

            if Vtot is None:
                raise ValueError("未找到体积信息：既没有 Lattice=... 也没有 V=/volume=...")

            # skip atom lines
            for _ in range(n):
                atom_line = f.readline()
                if not atom_line:
                    raise ValueError("缺少原子坐标行")

            V_list.append(Vtot / n)
            E_list.append(Etot / n)

    V = np.array(V_list, dtype=float)
    E = np.array(E_list, dtype=float)

    mask = np.isfinite(V) & np.isfinite(E)
    V, E = V[mask], E[mask]

    if len(V) < 4:
        raise ValueError(f"体积点太少（{len(V)}），BM3 拟合建议至少 4 个点")

    idx = np.argsort(V)
    return V[idx], E[idx]


# ===================== Main =====================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", help="包含多个 xyz/extxyz 文件的文件夹")
    ap.add_argument("--out", default="eos_multi.png", help="输出图片名")
    ap.add_argument("--table", default="eos_fit_table.csv", help="输出拟合参数表格 csv")
    ap.add_argument("--no-points", action="store_true", help="不画散点，只画拟合曲线")
    ap.add_argument("--only-points", action="store_true", help="只画散点，不拟合")
    args = ap.parse_args()

    folder = Path(args.folder)
    xyz_files = sorted(folder.glob("*.xyz"))
    if not xyz_files:
        raise SystemExit(f"未找到 .xyz 文件：{folder}")

    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">", "h", "H"]
    plt.figure()

    fit_rows = []
    EV_A3_to_GPa = 160.21766208

    valid_curve_count = 0

    for i, fp in enumerate(xyz_files):
        method = fp.stem

        try:
            V, E = read_xyz_or_extxyz_VE(fp)
        except Exception as e:
            print(f"[跳过] {fp.name}: {e}")
            continue

        mk = markers[valid_curve_count % len(markers)]
        valid_curve_count += 1

        # 画点
        if not args.no_points:
            plt.scatter(V, E, marker=mk, label=f"{method} (data)")

        # 拟合并画曲线 + 记录表格
        if not args.only_points:
            try:
                E0, V0, B0, Bp = fit_bm3(V, E)
            except Exception as e:
                print(f"[跳过拟合] {fp.name}: BM3 拟合失败：{e}")
                continue

            V_fit = np.linspace(V.min(), V.max(), 400)
            E_fit = birch_murnaghan_E(V_fit, E0, V0, B0, Bp)
            plt.plot(V_fit, E_fit, label=f"{method} (BM3)")

            fit_rows.append([
                method,
                V0, E0, B0, B0 * EV_A3_to_GPa,
                Bp
            ])

    if valid_curve_count == 0:
        raise SystemExit("没有任何有效的 xyz/extxyz 文件可绘图（都被跳过了）。")

    # 出图
    plt.xlabel(r"Volume/Atom ($\AA^3$)")
    plt.ylabel("Energy/Atom (eV)")
    plt.title("Equation of State (BM3 Fit)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.show()
    print(f"[OK] Saved figure: {args.out}")

    # 输出表格
    if fit_rows:
        fit_rows = sorted(fit_rows, key=lambda x: x[0])
        header = "method,V0(Angstrom^3/atom),E0(eV/atom),B0(eV/Angstrom^3),B0(GPa),Bp\n"
        with open(args.table, "w", encoding="utf-8") as f:
            f.write(header)
            for row in fit_rows:
                f.write(",".join(f"{v}" for v in row) + "\n")
        print(f"[OK] Saved table: {args.table}")
    else:
        print("[WARN] 没有任何方法成功拟合，因此未输出表格（可能你使用了 --only-points 或拟合都失败）。")


if __name__ == "__main__":
    main()
