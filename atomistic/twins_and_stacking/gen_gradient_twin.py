# -*- coding: utf-8 -*-
"""
=== 梯度孪晶铜结构生成器 ===
沿 Y 轴堆叠不同孪晶间距的区域, 形成梯度结构
可生成 2/3/4/N 梯度, 每个区域指定间距和层数
孪晶面垂直 Y 轴, 柱状晶粒 (晶界垂直孪晶面)

用法: python gen_gradient_twin.py
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

# --- 盒子截面 ---
BOX_X = 50.0    # X (Angstrom)
BOX_Z = 20.0    # Z (Angstrom)
                # Y 由各区域层总高自动计算

# --- 晶粒 ---
N_GRAINS = 2            # 柱状晶粒数
RANDOM_SEED = 42        # 固定种子保证各梯度结构的晶界位置一致

# --- 种子覆盖尺寸 ---
SEED_DUP_X = 5
SEED_DUP_Z = 3

# --- 梯度方案定义 ---
# 每个方案: (输出文件名, [(t1, n_layers1), (t2, n_layers2), ...])
#   t        = 该区域孪晶间距 (晶胞数, 1晶胞 ≈ 6.26 A)
#   n_layers = 该区域 Matrix+Twin 交替堆叠的层数
#              偶数层 = 完整周期, 奇数层最后一层无对应孪晶
#              建议: 间距小的区域多层(4~6), 间距大的少层(2)
#
# 原子数估算: ~85 * Y(A)  (对50x20截面)
# Y(A) = sum(t_i * n_i) * 6.26

GRADIENT_PLANS = [
    # === 2梯度 ===
    ("GNT_2grad_t1_t2",    [(1, 4), (2, 2)]),      # 密(4L) -> 疏(2L)
    ("GNT_2grad_t1_t3",    [(1, 4), (3, 2)]),      # 密(4L) -> 更疏(2L)

    # === 3梯度 ===
    ("GNT_3grad_t1_t2_t3", [(1, 4), (2, 2), (3, 2)]),   # 密->中->疏
    ("GNT_3grad_t1_t2_t4", [(1, 3), (2, 2), (4, 2)]),   # 密->中->最疏 (t4只1层也可写(4,1))

    # === 4梯度 ===
    ("GNT_4grad_t1_t2_t3_t4", [(1, 4), (2, 2), (3, 2), (4, 2)]),

    # === 自定义示例 ===
    # ("my_custom",  [(1,6), (3,4), (5,2), (8,2)]),  # 4梯度大跨度
    # ("symmetric",  [(1,2), (3,2), (3,2), (1,2)]),   # 对称夹心
]

# --- 输出 ---
OUT_DIR = r"D:\铜正式实验\5000以内\gradient_twin"

# ============================================================
#                    以下不需要修改
# ============================================================

ATOMSK = shutil.which("atomsk") or shutil.which("atomsk.exe")
if not ATOMSK:
    print("[ERROR] atomsk not found in PATH"); sys.exit(1)

CELL_Y = LATTICE_A * (3 ** 0.5)

def run(cmd, cwd=None):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        raise RuntimeError("atomsk rc=%d\n%s" % (r.returncode, r.stderr[-200:] if r.stderr else ""))

def read_atoms(path):
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and p[1] == "atoms":
                return int(p[0])
    return -1

def build_gradient(zone_spec, name, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    random.seed(RANDOM_SEED)

    with tempfile.TemporaryDirectory() as tmp:
        o1, o2, o3 = ORIENT
        mat = os.path.join(tmp, "mat.xsf")
        run([ATOMSK, "--create", "fcc", str(LATTICE_A), ELEMENT,
             "orient", o1, o2, o3, mat, "-overwrite"], cwd=tmp)
        twn = os.path.join(tmp, "twn.xsf")
        run([ATOMSK, mat, "-mirror", "0", "Y", "-wrap", twn, "-overwrite"], cwd=tmp)

        # 逐层构建: zone1 的所有层 + zone2 的所有层 + ...
        all_layers = []
        idx = 0
        total_y_cells = 0
        for t, nl in zone_spec:
            for i in range(nl):
                src = twn if (idx % 2 == 1) else mat
                fn = "layer_%03d.xsf" % idx
                run([ATOMSK, src, "-duplicate", str(SEED_DUP_X), str(t), str(SEED_DUP_Z),
                     os.path.join(tmp, fn), "-overwrite"], cwd=tmp)
                all_layers.append(fn)
                total_y_cells += t
                idx += 1

        # 合并种子
        seed = os.path.join(tmp, "seed.xsf")
        run([ATOMSK, "--merge", "Y", str(len(all_layers))] + all_layers +
            [seed, "-overwrite"], cwd=tmp)

        # polycrystal -> 柱状晶粒
        box_y = total_y_cells * CELL_Y
        poly_txt = os.path.join(tmp, "poly.txt")
        with open(poly_txt, "w") as f:
            f.write("box %s %s %s\n" % (BOX_X, box_y, BOX_Z))
            for _ in range(N_GRAINS):
                x = random.uniform(0, BOX_X)
                z = random.uniform(0, BOX_Z)
                beta = random.uniform(0, 360)
                f.write("node %.4f %.4f %.4f 0.00 %.2f 0.00\n" % (x, box_y/2, z, beta))

        out_lmp = os.path.join(out_dir, name + ".lmp")
        run([ATOMSK, "--polycrystal", seed, poly_txt, out_lmp, "-wrap", "-overwrite"], cwd=tmp)

    return out_lmp, read_atoms(out_lmp), box_y

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 70)
    print("  Gradient Nanotwin %s  |  Box: %.0f x Y x %.0f  |  Grains: %d"
          % (ELEMENT, BOX_X, BOX_Z, N_GRAINS))
    print("=" * 70)

    results = []
    for name, spec in GRADIENT_PLANS:
        desc = " -> ".join("t=%d(%dL)" % (t, nl) for t, nl in spec)
        total_cells = sum(t * nl for t, nl in spec)
        est_y = total_cells * CELL_Y
        est_atoms = int(85 * est_y)

        print("\n  [%s]" % name)
        print("    %s" % desc)
        print("    Y ~ %.0f A  |  est ~ %d atoms" % (est_y, est_atoms))

        try:
            path, na, by = build_gradient(spec, name, OUT_DIR)
            results.append((name, desc, na, by, "OK"))
            print("    => %d atoms, Y = %.1f A  OK" % (na, by))
        except Exception as e:
            results.append((name, desc, 0, 0, "FAIL"))
            print("    => FAIL: %s" % e)

    print("\n" + "=" * 70)
    print("  %-35s %-7s %-7s %s" % ("Structure", "Atoms", "Y(A)", "Gradient"))
    print("-" * 70)
    for name, desc, na, by, st in results:
        print("  %-35s %-7d %-7.1f %s" % (name, na, by, st))
    print("\n  Output: %s" % OUT_DIR)

if __name__ == "__main__":
    main()
