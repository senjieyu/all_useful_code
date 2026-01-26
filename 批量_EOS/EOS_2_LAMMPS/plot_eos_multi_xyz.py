#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_eos_multi_xyz.py

====================================================
用法 (Usage)
====================================================

功能：
1) 自动读取指定文件夹中所有 *.xyz (extxyz) 文件
2) 每个 xyz 文件代表一种方法
3) 从 extxyz 注释行读取：
      - Lattice="..."  (3x3 晶胞，用于计算体积)
      - energy=...     (DFT 总能量)
4) 计算并绘制 E/atom - V/atom
5) 对每种方法进行 BM3 (Birch-Murnaghan 3rd) 拟合，并输出参数到表格：
      V0, E0, B0 (以及可选 Bp)

----------------------------------------------------
运行示例
----------------------------------------------------

(1) 基本用法：画图 + 导出拟合参数表格
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
- eos_multi.png     (图像)
- eos_fit_table.csv (拟合参数表格)

表格列含义：
method, V0(Å^3/atom), E0(eV/atom), B0(eV/Å^3), B0(GPa), Bp

其中：
1 eV/Å^3 ≈ 160.21766208 GPa

----------------------------------------------------
extxyz 格式要求（关键）
----------------------------------------------------
每一帧 xyz 的第二行必须包含：
    Lattice="a1 a2 a3 b1 b2 b3 c1 c2 c3"
    energy=-123.456

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
    # 需要 scipy
    try:
        from scipy.optimize import curve_fit
    except ImportError as exc:
        raise ImportError("缺少 scipy：请先安装 `pip install scipy`") from exc

    i0 = int(np.argmin(E))
    p0 = [float(E[i0]), float(V[i0]), 0.5, 4.0]  # 初值：E0, V0, B0, Bp
    popt, _ = curve_fit(birch_murnaghan_E, V, E, p0=p0, maxfev=20000)
    return popt  # E0, V0, B0, Bp


# ===================== extxyz reader =====================
def parse_lattice(comment_line: str) -> float:
    m = re.search(r'Lattice="([^"]+)"', comment_line)
    if not m:
        raise ValueError("未找到 Lattice=\"...\"")
    vals = [float(x) for x in m.group(1).split()]
    if len(vals) != 9:
        raise ValueError("Lattice 不是 9 个数（3x3 晶格矩阵）")
    A = np.array(vals, dtype=float).reshape(3, 3)
    return float(abs(np.linalg.det(A)))


def parse_energy(comment_line: str) -> float:
    m = re.search(r'(?:(?:dft_)?energy|Energy|E)\s*=\s*([-\d\.eE+]+)', comment_line)
    if not m:
        raise ValueError("未找到 energy=...（请确认 extxyz 注释行包含能量）")
    return float(m.group(1))


def read_extxyz_VE(path: Path):
    V_list, E_list = [], []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue

            n = int(line)
            comment = f.readline()
            if not comment:
                break

            Vtot = parse_lattice(comment)
            Etot = parse_energy(comment)

            # skip atoms
            for _ in range(n):
                atom_line = f.readline()
                if not atom_line:
                    raise ValueError(f"{path.name}: 文件不完整，缺原子行")

            V_list.append(Vtot / n)
            E_list.append(Etot / n)

    V = np.array(V_list, dtype=float)
    E = np.array(E_list, dtype=float)

    mask = np.isfinite(V) & np.isfinite(E)
    V, E = V[mask], E[mask]

    if len(V) < 4:
        raise ValueError(f"{path.name}: 点数太少（{len(V)}），BM3 拟合建议至少 4 个点")

    idx = np.argsort(V)
    return V[idx], E[idx]


# ===================== main =====================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", help="包含多个 extxyz 的文件夹")
    ap.add_argument("--out", default="eos_multi.png", help="输出图片名")
    ap.add_argument("--table", default="eos_fit_table.csv", help="输出拟合参数表格 csv")
    ap.add_argument("--no-points", action="store_true", help="不画散点，只画拟合曲线")
    ap.add_argument("--only-points", action="store_true", help="只画散点，不拟合")
    args = ap.parse_args()

    folder = Path(args.folder)
    xyz_files = sorted(folder.glob("*.xyz"))
    if not xyz_files:
        raise SystemExit(f"未找到 xyz 文件：{folder}")

    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">", "h", "H"]
    plt.figure()

    fit_rows = []
    EV_A3_to_GPa = 160.21766208

    for i, fp in enumerate(xyz_files):
        method = fp.stem
        V, E = read_extxyz_VE(fp)
        mk = markers[i % len(markers)]

        # 画点
        if not args.no_points:
            plt.scatter(V, E, marker=mk, label=f"{method} (data)")

        # 拟合并画曲线 + 输出表格
        if not args.only_points:
            E0, V0, B0, Bp = fit_bm3(V, E)

            V_fit = np.linspace(V.min(), V.max(), 400)
            E_fit = birch_murnaghan_E(V_fit, E0, V0, B0, Bp)
            plt.plot(V_fit, E_fit, label=f"{method} (BM3)")

            fit_rows.append([
                method,
                V0, E0, B0, B0 * EV_A3_to_GPa,
                Bp
            ])

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
        fit_rows = sorted(fit_rows, key=lambda x: x[0])  # 按方法名排序
        header = "method,V0(Angstrom^3/atom),E0(eV/atom),B0(eV/Angstrom^3),B0(GPa),Bp\n"
        with open(args.table, "w", encoding="utf-8") as f:
            f.write(header)
            for row in fit_rows:
                f.write(",".join(f"{v}" for v in row) + "\n")
        print(f"[OK] Saved table: {args.table}")


if __name__ == "__main__":
    main()

