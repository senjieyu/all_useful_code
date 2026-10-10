#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
用法：
  python analyze_and_filter.py <dump文件夹> <功能编号...>

示例：
  python analyze_and_filter.py relaxed 1 2 3 4 5

功能编号说明：
  1 = 统计每个结构的力（Fmax / F95 / Frms / N），输出CSV到新文件夹
  2 = 画 Fmax 直方图（输出到新文件夹）
  3 = 画 F95 直方图（输出到新文件夹）
  4 = 画 F95 vs Fmax 散点图（输出到新文件夹）
  5 = 筛选（按分位数自动定阈值）+ 输出 accept/reject 名单 + 分文件夹

筛选策略（功能5）：
  - 剔除 F95 位于 top 3% 的结构（默认 q95=0.97）
  - 或剔除 Fmax 位于 top 1% 的结构（默认 qmax=0.99）

如果你想调严格程度：
  在脚本最下面 main() 里修改 Q_F95 / Q_FMAX 的值即可。
"""

import sys
import os
import glob
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ========== 这里可以改筛选严格程度 ==========
Q_F95 = 0.97   # F95 超过这个分位数 -> reject（top 3%）
Q_FMAX = 0.99  # Fmax 超过这个分位数 -> reject（top 1%）
# ===========================================


def read_dump_forces(path: str) -> np.ndarray:
    """读取 LAMMPS dump custom（单帧）里的力模长 |F|"""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header = None
    start = None
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            header = line.split()[2:]
            start = i + 1
            break
    if header is None:
        raise RuntimeError(f"{path}: 找不到 ITEM: ATOMS")

    data = []
    for j in range(start, len(lines)):
        if lines[j].startswith("ITEM:"):
            break
        s = lines[j].strip()
        if s:
            data.append(s.split())

    arr = np.array(data, dtype=float)
    idx = {h: i for i, h in enumerate(header)}

    # 优先用 v_fmag / fmag
    for k in ("v_fmag", "fmag"):
        if k in idx:
            return arr[:, idx[k]]

    # 否则由 fx fy fz 计算
    fx = arr[:, idx["fx"]]
    fy = arr[:, idx["fy"]]
    fz = arr[:, idx["fz"]]
    return np.sqrt(fx * fx + fy * fy + fz * fz)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    dump_dir = sys.argv[1]
    tasks = set(sys.argv[2:])

    dump_files = sorted(glob.glob(os.path.join(dump_dir, "*.force.dump")))
    if not dump_files:
        sys.exit(f"❌ 在目录 {dump_dir} 中没找到 *.force.dump")

    # 输出到新文件夹
    outdir = "analysis_results"
    os.makedirs(outdir, exist_ok=True)

    csv_path = os.path.join(outdir, "force_stats.csv")

    # ---------------- 功能 1：统计 ----------------
    if "1" in tasks:
        print("【1】统计每个结构的力数据 -> analysis_results/force_stats.csv")

        rows = []
        for fp in dump_files:
            fmag = read_dump_forces(fp)
            fmag = fmag[np.isfinite(fmag)]

            rows.append({
                "structure": os.path.basename(fp).replace(".force.dump", ""),
                "dump_path": fp,
                "N": int(len(fmag)),
                "Fmax": float(np.max(fmag)),
                "F95": float(np.percentile(fmag, 95)),
                "Frms": float(np.sqrt(np.mean(fmag ** 2))),
                "Fmean": float(np.mean(fmag)),
            })

        df = pd.DataFrame(rows).sort_values("Fmax", ascending=False).reset_index(drop=True)
        df.to_csv(csv_path, index=False)
        print(f"  -> 已输出 {csv_path}")

    # 如果没选 1，但后续要用 df，则从 CSV 读取
    df = pd.read_csv(csv_path) if os.path.exists(csv_path) else None
    if df is None:
        sys.exit("❌ 找不到 analysis_results/force_stats.csv，请先运行功能 1")

    # ---------------- 功能 2：画 Fmax 直方图 ----------------
    if "2" in tasks:
        print("【2】画 Fmax 直方图 -> analysis_results/hist_Fmax.png")
        plt.figure()
        plt.hist(df["Fmax"], bins=30)
        plt.xlabel("Fmax (eV/Å)")
        plt.ylabel("Count")
        plt.title("Distribution of Fmax")
        plt.savefig(os.path.join(outdir, "hist_Fmax.png"), dpi=200, bbox_inches="tight")
        plt.close()

    # ---------------- 功能 3：画 F95 直方图 ----------------
    if "3" in tasks:
        print("【3】画 F95 直方图 -> analysis_results/hist_F95.png")
        plt.figure()
        plt.hist(df["F95"], bins=30)
        plt.xlabel("F95 (eV/Å)")
        plt.ylabel("Count")
        plt.title("Distribution of F95")
        plt.savefig(os.path.join(outdir, "hist_F95.png"), dpi=200, bbox_inches="tight")
        plt.close()

    # ---------------- 功能 4：画散点图 ----------------
    if "4" in tasks:
        print("【4】画 F95 vs Fmax 散点图 -> analysis_results/scatter_F95_vs_Fmax.png")
        plt.figure()
        plt.scatter(df["F95"], df["Fmax"], s=15)
        plt.xlabel("F95 (eV/Å)")
        plt.ylabel("Fmax (eV/Å)")
        plt.title("F95 vs Fmax")
        plt.savefig(os.path.join(outdir, "scatter_F95_vs_Fmax.png"), dpi=200, bbox_inches="tight")
        plt.close()

    # ---------------- 功能 5：筛选 + 分文件夹 ----------------
    if "5" in tasks:
        print("【5】按分位数自动筛选 + 输出名单 + 分文件夹")

        thr_f95 = float(df["F95"].quantile(Q_F95))
        thr_fmax = float(df["Fmax"].quantile(Q_FMAX))

        df["reject"] = (df["F95"] > thr_f95) | (df["Fmax"] > thr_fmax)
        out_csv2 = os.path.join(outdir, "force_stats_with_reject.csv")
        df.to_csv(out_csv2, index=False)

        acc = df.loc[~df["reject"], "structure"].tolist()
        rej = df.loc[df["reject"], "structure"].tolist()

        with open(os.path.join(outdir, "accept_list.txt"), "w") as f:
            f.write("\n".join(acc) + "\n")
        with open(os.path.join(outdir, "reject_list.txt"), "w") as f:
            f.write("\n".join(rej) + "\n")

        print(f"  阈值：F95 > {thr_f95:.6g} (q={Q_F95}) 或 Fmax > {thr_fmax:.6g} (q={Q_FMAX})")
        print(f"  Accept={len(acc)}  Reject={len(rej)}")
        print(f"  -> 已输出 {out_csv2}, accept_list.txt, reject_list.txt")

        # 把 dump 文件分到新文件夹（这里用 copy，保险；你想 move 我也可以改）
        acc_dir = os.path.join(outdir, "accept_dumps")
        rej_dir = os.path.join(outdir, "reject_dumps")
        os.makedirs(acc_dir, exist_ok=True)
        os.makedirs(rej_dir, exist_ok=True)

        # 建立结构名 -> dump 路径映射
        name_to_path = {row["structure"]: row["dump_path"] for _, row in df.iterrows()}

        for name in acc:
            shutil.copy2(name_to_path[name], os.path.join(acc_dir, f"{name}.force.dump"))
        for name in rej:
            shutil.copy2(name_to_path[name], os.path.join(rej_dir, f"{name}.force.dump"))

        print(f"  -> dump 已复制到 {acc_dir}/ 和 {rej_dir}/")

    print("✅ 完成。所有输出都在 analysis_results/ 里。")


if __name__ == "__main__":
    main()

