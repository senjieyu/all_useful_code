# -*- coding: utf-8 -*-
"""
=== 纳米孪晶铜 + 晶界 结构生成器 ===
生成含晶界的纳米孪晶多晶结构
柱状晶粒沿 Y 方向贯穿 -> 晶界垂直于孪晶面
孪晶面垂直 Y 轴 ([111] 方向)

用法: python gen_nanotwin_gb.py
修改下方 CONFIG 区域的参数即可
"""
import os, subprocess, shutil, tempfile, random, sys

# ============================================================
#                    CONFIG - 修改这里
# ============================================================

# --- 材料 ---
ELEMENT = "Cu"
LATTICE_A = 3.615
ORIENT = ["[11-2]", "[111]", "[-110]"]

# --- 盒子截面 (X, Z 固定, Y 由孪晶参数决定) ---
BOX_X = 50.0    # X 方向盒子长度 (Angstrom)
BOX_Z = 20.0    # Z 方向盒子长度 (Angstrom)
MAX_Y = 80.0    # Y 方向最大长度 (Angstrom), 超过则跳过

# --- 晶粒 ---
N_GRAINS = 2            # 柱状晶粒数量 (晶界数 = N_GRAINS - 1)
RANDOM_SEED = None      # 随机种子, None=随机, 整数=可复现

# --- 种子尺寸 ---
SEED_DUP_X = 5          # 种子 X 方向重复数 (需覆盖 BOX_X)
SEED_DUP_Z = 3          # 种子 Z 方向重复数 (需覆盖 BOX_Z)

# --- 孪晶参数 ---
T_START = 1             # 起始孪晶间距 (晶胞数)
T_END   = 12            # 结束孪晶间距
T_STEP  = 1             # 间距步长

# --- 层数控制 ---
# 每个间距尝试的层数列表 (从多到少, 取第一个满足约束的)
# 层数 = Matrix + Twin 交替层总数, 2层=1个周期, 4层=2个周期
LAYER_OPTIONS = [4, 2]

# --- 原子数限制 ---
MAX_ATOMS = 5000

# --- 输出 ---
OUT_DIR = r"D:\铜正式实验\5000以内\nanotwin_with_GB"

# ============================================================
#                    以下不需要修改
# ============================================================

ATOMSK = shutil.which("atomsk") or shutil.which("atomsk.exe")
if not ATOMSK:
    print("[ERROR] atomsk not found in PATH"); sys.exit(1)

CELL_Y = LATTICE_A * (3 ** 0.5)

# 估算密度: FCC 原子数 / 盒子体积
# 精确值取决于晶向, 近似: BOX_X * BOX_Z * CELL_Y 体积内有 6 个原子(基元)
# density_per_A_Y ≈ BOX_X * BOX_Z / (ax * az) * 6 / CELL_Y * CELL_Y = BOX_X*BOX_Z*6/(ax*az)
# ax = LATTICE_A * sqrt(6), az = LATTICE_A * sqrt(2) for [11-2]/[-110]
_ax = LATTICE_A * (6 ** 0.5)
_az = LATTICE_A * (2 ** 0.5)
DENSITY_PER_Y = BOX_X * BOX_Z * 6 / (_ax * _az * CELL_Y)  # atoms per Angstrom of Y

def run(cmd, cwd=None):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError("atomsk rc=%d" % r.returncode)

def read_atoms(path):
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and p[1] == "atoms":
                return int(p[0])
    return -1

def build(t, nl, box_y, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)

    with tempfile.TemporaryDirectory() as tmp:
        o1, o2, o3 = ORIENT
        mat = os.path.join(tmp, "mat.xsf")
        run([ATOMSK, "--create", "fcc", str(LATTICE_A), ELEMENT,
             "orient", o1, o2, o3, mat, "-overwrite"], cwd=tmp)
        twn = os.path.join(tmp, "twn.xsf")
        run([ATOMSK, mat, "-mirror", "0", "Y", "-wrap", twn, "-overwrite"], cwd=tmp)

        layers = []
        for i in range(nl):
            src = twn if i % 2 else mat
            fn = "L%03d.xsf" % i
            run([ATOMSK, src, "-duplicate", str(SEED_DUP_X), str(t), str(SEED_DUP_Z),
                 os.path.join(tmp, fn), "-overwrite"], cwd=tmp)
            layers.append(fn)

        seed = os.path.join(tmp, "seed.xsf")
        run([ATOMSK, "--merge", "Y", str(len(layers))] + layers +
            [seed, "-overwrite"], cwd=tmp)

        poly_txt = os.path.join(tmp, "poly.txt")
        with open(poly_txt, "w") as f:
            f.write("box %s %s %s\n" % (BOX_X, box_y, BOX_Z))
            for _ in range(N_GRAINS):
                x = random.uniform(0, BOX_X)
                z = random.uniform(0, BOX_Z)
                beta = random.uniform(0, 360)
                f.write("node %.4f %.4f %.4f 0.00 %.2f 0.00\n" % (x, box_y/2, z, beta))

        out_lmp = os.path.join(out_dir, "nt_gb_t%02d_L%d.lmp" % (t, nl))
        run([ATOMSK, "--polycrystal", seed, poly_txt, out_lmp, "-wrap", "-overwrite"], cwd=tmp)

    return out_lmp, read_atoms(out_lmp)

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 65)
    print("  Nanotwin %s + GB  |  Box: %.0f x Y x %.0f  |  Grains: %d"
          % (ELEMENT, BOX_X, BOX_Z, N_GRAINS))
    print("=" * 65)

    results = []
    for t in range(T_START, T_END + 1, T_STEP):
        built = False
        for nl in sorted(LAYER_OPTIONS, reverse=True):  # 优先多层
            by = nl * t * CELL_Y
            if by > MAX_Y:
                continue
            est = DENSITY_PER_Y * by
            if est > MAX_ATOMS * 1.1:  # 10% 余量
                continue
            try:
                path, na = build(t, nl, by, OUT_DIR)
                if na > MAX_ATOMS:
                    os.remove(path)
                    continue
                results.append((t, nl, by, na, "OK"))
                print("  t=%-2d  layers=%d  Y=%5.1f A  atoms=%-5d  OK" % (t, nl, by, na))
                built = True
                break
            except Exception as e:
                results.append((t, nl, by, 0, "FAIL"))
                print("  t=%-2d  layers=%d  FAIL: %s" % (t, nl, e))
                built = True
                break
        if not built:
            print("  t=%-2d  SKIP (all options exceed limits)" % t)

    print("\n  Output: %s" % OUT_DIR)
    ok = sum(1 for r in results if r[4] == "OK")
    print("  Total: %d structures (%d OK)" % (len(results), ok))

if __name__ == "__main__":
    main()
