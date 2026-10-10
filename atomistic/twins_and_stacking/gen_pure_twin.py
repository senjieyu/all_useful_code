# -*- coding: utf-8 -*-
"""
=== 纯孪晶铜结构生成器 ===
生成无晶界的纯孪晶单晶结构 (Matrix-Twin 交替堆垛)
孪晶面垂直 Y 轴 ([111] 方向)

用法: python gen_pure_twin.py
修改下方 CONFIG 区域的参数即可
"""
import os, subprocess, shutil, tempfile, sys

# ============================================================
#                    CONFIG - 修改这里
# ============================================================

# --- 材料 ---
ELEMENT = "Cu"
LATTICE_A = 3.615                          # 晶格常数 (Angstrom)
ORIENT = ["[11-2]", "[111]", "[-110]"]     # 晶向: X, Y(堆垛方向), Z

# --- 盒子尺寸 (通过晶胞重复数控制) ---
DUP_X = 8       # X 方向晶胞重复数 (实际长度 = DUP_X * ax)
DUP_Z = 4       # Z 方向晶胞重复数 (实际长度 = DUP_Z * az)
                 # Y 方向由孪晶间距和周期数自动决定

# --- 孪晶参数 ---
T_START = 1     # 起始孪晶间距 (晶胞数)
T_END   = 11    # 结束孪晶间距 (晶胞数)
T_STEP  = 1     # 间距步长

# --- 周期数控制 ---
# "auto": 自动计算使原子数不超过 MAX_ATOMS
# 整数: 所有结构使用固定周期数
N_PERIODS = "auto"
MAX_ATOMS = 5000    # 仅 auto 模式生效
MAX_PERIODS = 8     # auto 模式下的周期数上限

# --- 输出 ---
OUT_DIR = r"D:\铜正式实验\5000以内\pure_twin_series"

# ============================================================
#                    以下不需要修改
# ============================================================

ATOMSK = shutil.which("atomsk") or shutil.which("atomsk.exe")
if not ATOMSK:
    print("[ERROR] atomsk not found in PATH"); sys.exit(1)

CELL_Y = LATTICE_A * (3 ** 0.5)  # [111] 方向单胞高度

def calc_periods(t):
    if N_PERIODS != "auto":
        return int(N_PERIODS)
    atoms_per_period = DUP_X * 2 * t * DUP_Z * 6
    n = max(1, MAX_ATOMS // atoms_per_period)
    return min(n, MAX_PERIODS)

def run(cmd, cwd=None):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError("atomsk rc=%d\n%s" % (r.returncode, r.stderr[-200:] if r.stderr else ""))

def read_atoms(path):
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and p[1] == "atoms":
                return int(p[0])
    return -1

def build(t, out_dir):
    np_ = calc_periods(t)
    os.makedirs(out_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        o1, o2, o3 = ORIENT
        mat = os.path.join(tmp, "mat.xsf")
        run([ATOMSK, "--create", "fcc", str(LATTICE_A), ELEMENT,
             "orient", o1, o2, o3, mat, "-overwrite"], cwd=tmp)

        twn = os.path.join(tmp, "twn.xsf")
        run([ATOMSK, mat, "-mirror", "0", "Y", "-wrap", twn, "-overwrite"], cwd=tmp)

        layers = []
        for i in range(np_ * 2):
            src = twn if i % 2 else mat
            fn = "L%03d.xsf" % i
            run([ATOMSK, src, "-duplicate", str(DUP_X), str(t), str(DUP_Z),
                 os.path.join(tmp, fn), "-overwrite"], cwd=tmp)
            layers.append(fn)

        merged = os.path.join(tmp, "merged.xsf")
        run([ATOMSK, "--merge", "Y", str(len(layers))] + layers +
            [merged, "-overwrite"], cwd=tmp)

        out = os.path.join(out_dir, "twin_t%02d.lmp" % t)
        run([ATOMSK, merged, out, "-overwrite"], cwd=tmp)

    return out, read_atoms(out), np_

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    spacings = list(range(T_START, T_END + 1, T_STEP))

    print("=" * 60)
    print("  Pure Twin %s  |  DUP: %dx(t)x%d  |  t=%d..%d"
          % (ELEMENT, DUP_X, DUP_Z, T_START, T_END))
    print("=" * 60)

    results = []
    for t in spacings:
        try:
            path, na, np_ = build(t, OUT_DIR)
            results.append((t, na, np_, "OK"))
            print("  t=%-3d  spacing=%5.1f A  periods=%d  atoms=%-5d  OK" % (t, t*CELL_Y, np_, na))
        except Exception as e:
            results.append((t, 0, 0, "FAIL"))
            print("  t=%-3d  FAIL: %s" % (t, e))

    print("\n  Output: %s" % OUT_DIR)
    print("  Total: %d structures" % sum(1 for r in results if r[3] == "OK"))

if __name__ == "__main__":
    main()
