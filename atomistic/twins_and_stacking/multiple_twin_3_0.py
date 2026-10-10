import os
import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import make_supercell


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


def make_fcc111_multictb_stacking_lists(
    d_up: list,
    d_middle: int,
    d_down: list,
    start: str = "A",
    build_from: str = "bottom",  # build stacking bottom->top or top->bottom (just for convenience)
):
    """
    你的新定义（列表输入）：

    - d_up  : 上部的段长度列表（从上往下数） 例如 [du1, du2, du3, ...]
    - d_down: 下部的段长度列表（从上往下数） 例如 [dd1, dd2, dd3, ...]
    - d_middle: 上下两组之间的中间段长度

    结构从“上到下”的段序列为：
        segments_topdown = d_up + [d_middle] + d_down

    段与段之间都是 CTB（孪晶界），每过一个 CTB 翻转一次堆垛方向。
    CTB 使用“共享单层”实现：相邻两段的连接层只出现一次。

    重要：如果 segments 有 N 段，则 CTB 数 = N-1
          总层数 = sum(segments) - (N-1)
    """
    if not isinstance(d_middle, int) or d_middle <= 0:
        raise ValueError("d_middle must be a positive integer (layer count).")

    if not isinstance(d_up, (list, tuple)) or len(d_up) == 0:
        raise ValueError("d_up must be a non-empty list of positive integers.")
    if not isinstance(d_down, (list, tuple)) or len(d_down) == 0:
        raise ValueError("d_down must be a non-empty list of positive integers.")

    for i, v in enumerate(d_up, 1):
        if not isinstance(v, int) or v <= 0:
            raise ValueError(f"d_up[{i}] must be a positive integer.")
    for i, v in enumerate(d_down, 1):
        if not isinstance(v, int) or v <= 0:
            raise ValueError(f"d_down[{i}] must be a positive integer.")

    start = start.strip().upper()
    if start not in "ABC":
        raise ValueError("start must be 'A'/'B'/'C'.")

    fwd = {"A": "B", "B": "C", "C": "A"}  # matrix: ABCABC...
    bwd = {"A": "C", "C": "B", "B": "A"}  # twin:   ACBACB...

    segments_topdown = list(d_up) + [d_middle] + list(d_down)

    if build_from.lower() == "bottom":
        segments = list(reversed(segments_topdown))  # bottom->top
    elif build_from.lower() == "top":
        segments = segments_topdown[:]               # top->bottom
    else:
        raise ValueError("build_from must be 'bottom' or 'top'.")

    seq = []
    sense = "fwd"
    next_map = fwd

    seq += _grow_layers(segments[0], next_map, start)

    for seg_len in segments[1:]:
        # flip at each CTB
        sense = "bwd" if sense == "fwd" else "fwd"
        next_map = fwd if sense == "fwd" else bwd

        shared = seq[-1]  # CTB shared single layer
        seg_letters = _grow_layers(seg_len, next_map, shared)
        seq += seg_letters[1:]  # skip shared

    return "".join(seq), segments_topdown


def build_cu111_multictb_lists(
    d_up, d_middle, d_down,
    nx=4, ny=4, a=3.615,
    start="A",
    build_from="bottom",
):
    stacking, segments_topdown = make_fcc111_multictb_stacking_lists(
        d_up=d_up, d_middle=d_middle, d_down=d_down,
        start=start, build_from=build_from
    )
    atoms = build_cu111_from_stacking(stacking=stacking, nx=nx, ny=ny, a=a)
    return atoms, stacking, segments_topdown


def to_orthorhombic_cell_fixed_111(atoms: Atoms):
    """
    对 fcc(111) 面内六方(60°)晶格的稳定正交化超胞：
      A1 =  a1 + a2
      A2 = -a1 + a2
      A3 =  a3
    """
    P = np.array([[1, 1, 0],
                  [-1, 1, 0],
                  [0,  0, 1]], dtype=int)
    ortho = make_supercell(atoms, P)
    ortho.wrap(eps=1e-12)
    return ortho, P


if __name__ == "__main__":
    # ======================
    # 新输入方式：列表（有几个数字就是几段）
    # 上部（从上往下的段长度）
    d_up = [2, 4, 4]     #FCC层数为d1-1 d2-2 d3-2 ..     # 3段 => 2条CTB（在上部内部）
    # 中间段
    d_middle = 5
    # 下部（从上往下的段长度）
    d_down = [4, 4, 4, 2]     # #FCC层数为d1-2 d2-2...dn-1 4段 => 3条CTB（在下部内部）
    # ======================

    nx, ny = 5, 5
    a = 3.615

    rect_cell = True  # 是否转正交晶胞
    start = "A"

    atoms, stacking, segments_topdown = build_cu111_multictb_lists(
        d_up=d_up, d_middle=d_middle, d_down=d_down,
        nx=nx, ny=ny, a=a,
        start=start,
        build_from="bottom",
    )

    P_used = None
    if rect_cell:
        atoms, P_used = to_orthorhombic_cell_fixed_111(atoms)

    # 输出
    out_dir = "Cu111_multiCTB_list"
    os.makedirs(out_dir, exist_ok=True)

    # 文件名带上列表信息（用-连接）
    up_tag = "-".join(map(str, d_up))
    down_tag = "-".join(map(str, d_down))
    tag = "_rect" if rect_cell else "_hex"

    vasp_path = os.path.join(out_dir, f"up[{up_tag}]_mid{d_middle}_down[{down_tag}]{tag}_POSCAR.vasp")
    write(vasp_path, atoms, format="vasp")

    # 统计
    Nseg = len(d_up) + 1 + len(d_down)
    Nctb = Nseg - 1
    total_expected = sum(d_up) + d_middle + sum(d_down) - Nctb

    print("Done:")
    print(f"  d_up     = {d_up}  (segments={len(d_up)}, CTB inside up={len(d_up)-1})")
    print(f"  d_middle = {d_middle}")
    print(f"  d_down   = {d_down}  (segments={len(d_down)}, CTB inside down={len(d_down)-1})")
    print(f"  segments_topdown = {segments_topdown}")
    print(f"  CTB total = {Nctb}")
    print(f"  Layers   = {len(stacking)} (expected {total_expected})")
    print(f"  Atoms    = {len(atoms)}")
    if P_used is not None:
        print(f"  Rect supercell P:\n{P_used}")
    print(f"  Output   = {vasp_path}")
