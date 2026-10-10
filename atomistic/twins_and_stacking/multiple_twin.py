import os
import numpy as np
from ase import Atoms
from ase.io import write


def build_cu111_from_stacking(stacking: str, nx: int = 4, ny: int = 4, a: float = 3.615):
    stacking = stacking.strip().upper()
    if len(stacking) == 0 or any(ch not in "ABC" for ch in stacking):
        raise ValueError("stacking must be a non-empty string containing only A/B/C")

    nlayer = len(stacking)

    a_hex = a / np.sqrt(2.0)
    d111 = a / np.sqrt(3.0)

    a1 = np.array([a_hex, 0.0, 0.0])
    a2 = np.array([0.5 * a_hex, 0.5 * np.sqrt(3.0) * a_hex, 0.0])
    a3 = np.array([0.0, 0.0, nlayer * d111])  # periodic along z, no vacuum

    shift = (a1 + a2) / 3.0
    shifts = {"A": 0.0 * shift, "B": 1.0 * shift, "C": 2.0 * shift}

    positions, symbols = [], []
    for k, letter in enumerate(stacking):
        lateral = shifts[letter]
        z = k * d111
        for i in range(nx):
            for j in range(ny):
                r = i * a1 + j * a2 + lateral + np.array([0.0, 0.0, z])
                positions.append(r)
                symbols.append("Cu")

    atoms = Atoms(
        symbols=symbols,
        positions=np.array(positions, dtype=float),
        cell=np.vstack([nx * a1, ny * a2, a3]),
        pbc=(True, True, True),
    )
    atoms.wrap(eps=1e-12)
    return atoms


def _grow_layers(n, next_map, first_letter):
    """Generate n letters starting from first_letter using next_map."""
    out = [first_letter]
    cur = first_letter
    for _ in range(n - 1):
        cur = next_map[cur]
        out.append(cur)
    return out


def make_fcc111_multictb_stacking(
    d1: int, d2: int, d3: int,
    d_middle: int,
    d4: int, d5: int, d6: int,
    start: str = "A",
    build_from: str = "bottom",   # "bottom" or "top"
):
    """
    生成你图里那种：上3条CTB + 中间隔层 + 下3条CTB 的(111)堆垛序列（A/B/C）。
    每条CTB都是“单层共享面”（不会重复插层）。

    解释（按从上往下的段长）：
      seg_top0 = d1
      seg_top1 = d2
      seg_top2 = d3
      seg_mid  = d_middle
      seg_bot0 = d4
      seg_bot1 = d5
      seg_bot2 = d6

    总段数 = 7 段；段与段之间有 6 个 CTB（每过一个 CTB 翻转一次堆垛方向）。

    注意：由于CTB为共享层，
      总层数 = sum(segments) - number_of_CTBs = (d1+d2+d3+d_middle+d4+d5+d6) - 6
      （因为每个界面共享一层）
    """
    vals = [d1, d2, d3, d_middle, d4, d5, d6]
    names = ["d1", "d2", "d3", "d_middle", "d4", "d5", "d6"]
    for name, v in zip(names, vals):
        if not isinstance(v, int) or v <= 0:
            raise ValueError(f"{name} must be a positive integer (layer count).")

    start = start.strip().upper()
    if start not in "ABC":
        raise ValueError("start must be 'A'/'B'/'C'.")

    # matrix / twin mapping
    fwd = {"A": "B", "B": "C", "C": "A"}  # matrix: ABCABC...
    bwd = {"A": "C", "C": "B", "B": "A"}  # twin:   ACBACB...

    # 你图是“上面一组 + middle + 下面一组”
    segments_topdown = [d1, d2, d3, d_middle, d4, d5, d6]

    # 我们实际建模通常是 bottom->top 叠层更直观
    if build_from.lower() == "bottom":
        segments = list(reversed(segments_topdown))  # bottom->top: d6,d5,d4,d_middle,d3,d2,d1
    elif build_from.lower() == "top":
        segments = segments_topdown[:]               # top->bottom
    else:
        raise ValueError("build_from must be 'bottom' or 'top'.")

    seq = []
    sense = "fwd"  # start with matrix by default（你也可以改成 bwd）
    next_map = fwd if sense == "fwd" else bwd

    # 第一段：直接 grow segments[0]
    seq += _grow_layers(segments[0], next_map, start)

    # 后续段：每遇到一个 CTB 翻转一次 sense，且首层共享（单层CTB）
    for seg_len in segments[1:]:
        # flip stacking sense at each CTB
        sense = "bwd" if sense == "fwd" else "fwd"
        next_map = fwd if sense == "fwd" else bwd

        shared = seq[-1]  # CTB 单层共享面
        seg_letters = _grow_layers(seg_len, next_map, shared)
        seq += seg_letters[1:]  # skip first because shared

    return "".join(seq)


def build_cu111_multictb(
    d1, d2, d3, d_middle, d4, d5, d6,
    nx=4, ny=4, a=3.615,
    start="A",
    build_from="bottom",
):
    stacking = make_fcc111_multictb_stacking(
        d1=d1, d2=d2, d3=d3,
        d_middle=d_middle,
        d4=d4, d5=d5, d6=d6,
        start=start,
        build_from=build_from,
    )
    atoms = build_cu111_from_stacking(stacking=stacking, nx=nx, ny=ny, a=a)
    return atoms, stacking


if __name__ == "__main__":
    # ======================
    # 你图里的可调参数（层数）
    # 顶部三条孪晶界间距：d1,d2,d3 (中间的间隔为边界的间隔为d-1，d-2)
    d1, d2, d3 = 5, 4, 3

    # 两组之间的距离（层数）
    d_middle = 10

    # 底部三条孪晶界间距：d4,d5,d6 (中间的间隔为边界的间隔为d-1，d-2)
    d4, d5, d6 = 5, 4, 5
    # ======================

    # 面内尺寸
    nx, ny = 5, 5
    a = 3.615

    # 从 bottom->top 构建（推荐）
    atoms, stacking = build_cu111_multictb(
        d1, d2, d3, d_middle, d4, d5, d6,
        nx=nx, ny=ny, a=a,
        start="A",
        build_from="bottom",
    )

    # 输出到新文件夹，只保留 .vasp
    out_dir = "Cu111_multiCTB_"
    os.makedirs(out_dir, exist_ok=True)

    vasp_path = os.path.join(out_dir, f"d1{d1}_d2{d2}_d3{d3}_mid{d_middle}_d4{d4}_d5{d5}_d6{d6}POSCAR.vasp")
    write(vasp_path, atoms, format="vasp")

    # 打印信息
    total_expected = (d1 + d2 + d3 + d_middle + d4 + d5 + d6) - 6
    print("Done:")
    print(f"  d1 d2 d3 = {d1} {d2} {d3}")
    print(f"  d_middle = {d_middle}")
    print(f"  d4 d5 d6 = {d4} {d5} {d6}")
    print(f"  Stacking = {stacking}")
    print(f"  Layers   = {len(stacking)} (expected {total_expected})")
    print(f"  Atoms    = {len(atoms)}")
    print(f"  Output   = {vasp_path}")
