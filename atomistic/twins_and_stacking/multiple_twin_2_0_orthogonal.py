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
    build_from: str = "bottom",
):
    vals = [d1, d2, d3, d_middle, d4, d5, d6]
    names = ["d1", "d2", "d3", "d_middle", "d4", "d5", "d6"]
    for name, v in zip(names, vals):
        if not isinstance(v, int) or v <= 0:
            raise ValueError(f"{name} must be a positive integer (layer count).")

    start = start.strip().upper()
    if start not in "ABC":
        raise ValueError("start must be 'A'/'B'/'C'.")

    fwd = {"A": "B", "B": "C", "C": "A"}  # matrix: ABCABC...
    bwd = {"A": "C", "C": "B", "B": "A"}  # twin:   ACBACB...

    segments_topdown = [d1, d2, d3, d_middle, d4, d5, d6]

    if build_from.lower() == "bottom":
        segments = list(reversed(segments_topdown))
    elif build_from.lower() == "top":
        segments = segments_topdown[:]
    else:
        raise ValueError("build_from must be 'bottom' or 'top'.")

    seq = []
    sense = "fwd"
    next_map = fwd

    seq += _grow_layers(segments[0], next_map, start)

    for seg_len in segments[1:]:
        sense = "bwd" if sense == "fwd" else "fwd"
        next_map = fwd if sense == "fwd" else bwd

        shared = seq[-1]  # CTB 单层共享面
        seg_letters = _grow_layers(seg_len, next_map, shared)
        seq += seg_letters[1:]  # skip shared

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


# =============================
# 新增：转正交晶胞接口（物理严格）
# =============================
def to_orthorhombic_cell(atoms: Atoms, keep_c_axis: bool = True):
    """
    将当前晶胞转换为正交晶胞（三边90°），采用整数超胞变换 make_supercell。
    对于 fcc(111) 这种六方(60°)面内晶格，正交晶胞必然是“超胞”。

    keep_c_axis=True:
        - 固定 a3 不变（更适合(111)堆垛/孪晶体系）
        - 仅对面内 (a1,a2) 做正交化超胞：
              A1 =  a1 + a2
              A2 = -a1 + a2
              A3 =  a3
        这个变换得到面内90°且与原晶格严格周期等价。

    keep_c_axis=False:
        - 这里仍然给同样的 P（最稳），你也可以扩展成自动搜索版本。
    """
    # 对(111)面内六方最经典、最稳定的正交超胞变换
    P = np.array([[1, 1, 0],
                  [-1, 1, 0],
                  [0,  0, 1]], dtype=int)

    # keep_c_axis=True 时 P 已经保证 A3=a3
    ortho = make_supercell(atoms, P)
    ortho.wrap(eps=1e-12)
    return ortho, P


if __name__ == "__main__":
    # ======================
    # 你图里的可调参数（层数）
    d1, d2, d3 = 5, 4, 3
    d_middle = 10
    d4, d5, d6 = 5, 4, 5
    # ======================

    nx, ny = 5, 5
    a = 3.615

    # ===== 新增开关：是否转正交晶胞 =====
    rect_cell = True          # True: 输出正交超胞；False: 输出原六方晶胞
    keep_c_axis = True        # 推荐 True（保持(111)堆垛方向为z）

    atoms, stacking = build_cu111_multictb(
        d1, d2, d3, d_middle, d4, d5, d6,
        nx=nx, ny=ny, a=a,
        start="A",
        build_from="bottom",
    )

    P_used = None
    if rect_cell:
        atoms, P_used = to_orthorhombic_cell(atoms, keep_c_axis=keep_c_axis)

    # 输出到新文件夹，只保留 .vasp
    out_dir = "Cu111_multiCTB_orthogonal"
    os.makedirs(out_dir, exist_ok=True)

    tag = "_rect" if rect_cell else "_hex"
    vasp_path = os.path.join(
        out_dir,
        f"d1_{d1}_d2_{d2}_d3_{d3}_mid_{d_middle}_d4_{d4}_d5_{d5}_d6_{d6}{tag}_POSCAR.vasp"
    )
    write(vasp_path, atoms, format="vasp")

    total_expected = (d1 + d2 + d3 + d_middle + d4 + d5 + d6) - 6

    print("Done:")
    print(f"  rect_cell = {rect_cell} | keep_c_axis = {keep_c_axis}")
    print(f"  d1 d2 d3 = {d1} {d2} {d3}")
    print(f"  d_middle = {d_middle}")
    print(f"  d4 d5 d6 = {d4} {d5} {d6}")
    print(f"  Stacking = {stacking}")
    print(f"  Layers   = {len(stacking)} (expected {total_expected})")
    print(f"  Atoms    = {len(atoms)}")
    if P_used is not None:
        print(f"  Supercell transform P:\n{P_used}")
    print(f"  Output   = {vasp_path}")
