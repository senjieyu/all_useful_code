#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===========================================================
转成正交超胞.py
-----------------------------------------------------------
自动将任意 VASP 超胞（任意角度）转换为“尽量正交”的整数超胞（角度接近 90°）
使用 ase.build.make_supercell（整数变换），保证周期严格正确。

【输入】POSCAR / CONTCAR / *.vasp 或包含它们的文件夹
【输出】<输入名>_ortho/ 里仅输出 *.vasp

【推荐（层状/堆垛/孪晶体系）】
  python 转成正交超胞.py POSCAR --keep-c-axis

【更强搜索】
  python 转成正交超胞.py POSCAR --keep-c-axis --max-index 6

===========================================================
"""

import argparse
import numpy as np
from pathlib import Path
from itertools import combinations
from ase.io import read, write
from ase.build import make_supercell


def is_vasp_like(path: Path) -> bool:
    name = path.name.lower()
    return (path.suffix.lower() == ".vasp") or (name in ("poscar", "contcar"))


def collect_inputs(input_path: Path):
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        return [p for p in sorted(input_path.iterdir()) if p.is_file() and is_vasp_like(p)]
    raise FileNotFoundError(f"Input path not found: {input_path}")


def make_out_dir(input_path: Path, out_dir: str | None):
    if out_dir is not None:
        outp = Path(out_dir).expanduser().resolve()
    else:
        outp = (input_path.parent / f"{input_path.name}_ortho").resolve()
    outp.mkdir(parents=True, exist_ok=True)
    return outp


def out_name_for(infile: Path) -> str:
    low = infile.name.lower()
    if low in ("poscar", "contcar"):
        return infile.name + ".vasp"
    if infile.suffix.lower() == ".vasp":
        return infile.stem + "_ortho.vasp"
    return infile.name + "_ortho.vasp"


def angle_deg(u, v) -> float:
    uu = np.linalg.norm(u)
    vv = np.linalg.norm(v)
    if uu < 1e-12 or vv < 1e-12:
        return 180.0
    c = np.dot(u, v) / (uu * vv)
    c = np.clip(c, -1.0, 1.0)
    return float(np.degrees(np.arccos(c)))


def find_integer_orthorhombic_P(cell, max_index=4, tol_deg=0.5, keep_c_axis=False,
                                pair_gate_deg=10.0, max_vecs=220, max_pairs=6000):
    """
    找整数超胞矩阵 P，使新晶胞尽量接近 90°。
    评分：key = (err_sum, detP)，err_sum=|a12-90|+|a13-90|+|a23-90|

    为避免暴力 C(N,3)，使用“先找近似正交对，再找第三向量”的剪枝策略。
    - pair_gate_deg: 允许作为“近似正交对”的阈值（越小越快，但可能错过解）
    - max_vecs:      候选向量上限（越小越快）
    - max_pairs:     候选对上限（越小越快）
    """
    a1, a2, a3 = cell[0], cell[1], cell[2]

    rng = range(-max_index, max_index + 1)

    # 生成候选整数向量
    vecs = []
    for i in rng:
        for j in rng:
            for k in rng:
                if i == 0 and j == 0 and k == 0:
                    continue
                # keep_c_axis 时，第三向量固定为 (0,0,1)，这里仍然允许其它向量含 k 分量
                v = i * a1 + j * a2 + k * a3
                lv = np.linalg.norm(v)
                if lv > 1e-6:
                    vecs.append(((i, j, k), v, lv))

    # 优先短向量（更可能得到小 detP）
    vecs.sort(key=lambda x: x[2])
    vecs = vecs[:min(len(vecs), max_vecs)]

    def eval_triplet(c1, v1, c2, v2, c3, v3):
        P = np.array([c1, c2, c3], dtype=int)
        detP = int(round(np.linalg.det(P)))
        if detP <= 0:
            return None
        a12 = angle_deg(v1, v2)
        a13 = angle_deg(v1, v3)
        a23 = angle_deg(v2, v3)
        err = abs(a12 - 90.0) + abs(a13 - 90.0) + abs(a23 - 90.0)
        return err, (a12, a13, a23), detP, P

    best = None  # (key, err, angs, detP, P)

    if keep_c_axis:
        # 固定第三向量为 a3
        c3 = (0, 0, 1)
        v3 = a3
        for (c1, v1, _), (c2, v2, _) in combinations(vecs, 2):
            res = eval_triplet(c1, v1, c2, v2, c3, v3)
            if res is None:
                continue
            err, angs, detP, P = res
            key = (err, detP)
            if best is None or key < best[0]:
                best = (key, err, angs, detP, P)

    else:
        # 先建“近似正交对”
        pairs = []
        for (c1, v1, _,), (c2, v2, _,) in combinations(vecs, 2):
            a12 = angle_deg(v1, v2)
            if abs(a12 - 90.0) <= pair_gate_deg:
                pairs.append((c1, v1, c2, v2, abs(a12 - 90.0)))

        # 只保留最接近 90° 的 pairs
        pairs.sort(key=lambda x: x[4])
        pairs = pairs[:min(len(pairs), max_pairs)]

        # 对每个 pair 找第三向量，同时近似正交
        for c1, v1, c2, v2, _ in pairs:
            for (c3, v3, _) in vecs:
                # 快速角度门限（再算 det 之前先筛）
                if abs(angle_deg(v1, v3) - 90.0) > pair_gate_deg:
                    continue
                if abs(angle_deg(v2, v3) - 90.0) > pair_gate_deg:
                    continue
                res = eval_triplet(c1, v1, c2, v2, c3, v3)
                if res is None:
                    continue
                err, angs, detP, P = res
                key = (err, detP)
                if best is None or key < best[0]:
                    best = (key, err, angs, detP, P)

    if best is None:
        raise RuntimeError("Failed to find (near) orthorhombic supercell. Try larger --max-index or --keep-c-axis.")

    _, err, angs, detP, P = best
    _ = tol_deg  # 仅作参考，不做硬阈值
    return err, angs, detP, P


def convert_one(infile: Path, out_dir: Path, fmt_hint=None, max_index=4, tol_deg=0.5, keep_c_axis=False):
    atoms = read(str(infile), format=fmt_hint) if fmt_hint else read(str(infile))

    err, angs, detP, P = find_integer_orthorhombic_P(
        atoms.cell.array, max_index=max_index, tol_deg=tol_deg, keep_c_axis=keep_c_axis
    )

    atoms_o = make_supercell(atoms, P)
    atoms_o.wrap(eps=1e-12)

    out_file = out_dir / out_name_for(infile)
    write(str(out_file), atoms_o, format="vasp")

    return {
        "infile": infile,
        "out_file": out_file,
        "P": P,
        "err": err,
        "angs": angs,
        "detP": detP,
        "atoms_in": len(atoms),
        "atoms_out": len(atoms_o),
    }


def main():
    ap = argparse.ArgumentParser(description="Convert any VASP cell to (near) orthorhombic via integer supercell search.")
    ap.add_argument("input", help="Input file (POSCAR/CONTCAR/*.vasp) OR directory containing them.")
    ap.add_argument("--out-dir", default=None, help="Output directory. Default: <input>_ortho")
    ap.add_argument("--format", default=None, help="ASE format hint (e.g. vasp). Usually not needed.")
    ap.add_argument("--max-index", type=int, default=4, help="Search coefficient range [-N,N]. Default 4.")
    ap.add_argument("--tol-deg", type=float, default=0.5, help="Angle tolerance reference in degrees. Default 0.5.")
    ap.add_argument("--keep-c-axis", action="store_true", help="Keep c-axis unchanged (use original a3). Recommended for layered/stacked systems.")

    args = ap.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    files = collect_inputs(input_path)
    if len(files) == 0:
        raise RuntimeError(f"No POSCAR/CONTCAR/*.vasp found in: {input_path}")

    out_dir = make_out_dir(input_path, args.out_dir)

    print(f"[INFO] Input      : {input_path}")
    print(f"[INFO] Found      : {len(files)} file(s)")
    print(f"[INFO] Output dir : {out_dir}")
    print(f"[INFO] max_index  : {args.max_index} | tol_deg(ref): {args.tol_deg} | keep_c_axis: {args.keep_c_axis}")
    print("")

    for f in files:
        info = convert_one(
            f, out_dir,
            fmt_hint=args.format,
            max_index=args.max_index,
            tol_deg=args.tol_deg,
            keep_c_axis=args.keep_c_axis
        )
        a12, a13, a23 = info["angs"]
        print(f"[OK] {info['infile'].name} -> {info['out_file'].name}")
        print(f"     det(P)={info['detP']}  atoms {info['atoms_in']} -> {info['atoms_out']}")
        print(f"     angles(deg)=({a12:.4f}, {a13:.4f}, {a23:.4f})  err_sum={info['err']:.4f}")
        print(f"     P=\n{info['P']}\n")

    print("Done.")


if __name__ == "__main__":
    main()
