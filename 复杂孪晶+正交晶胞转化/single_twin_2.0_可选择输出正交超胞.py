import os
import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import make_supercell


def build_cu111_from_stacking(stacking: str, nx: int = 4, ny: int = 4, a: float = 3.615):
    """
    Build periodic Cu close-packed (111) structure with arbitrary stacking sequence.
    - stacking: e.g. "ABCABC", "ACB", ...
    - nx, ny: in-plane repeats of the (111) primitive net
    - No vacuum: periodic in x, y, z
    """
    stacking = stacking.strip().upper()
    if len(stacking) == 0 or any(ch not in "ABC" for ch in stacking):
        raise ValueError("stacking must be a non-empty string containing only A/B/C")

    nlayer = len(stacking)

    a_hex = a / np.sqrt(2.0)
    d111 = a / np.sqrt(3.0)

    # (111) plane primitive vectors (60 deg)
    a1 = np.array([a_hex, 0.0, 0.0])
    a2 = np.array([0.5 * a_hex, 0.5 * np.sqrt(3.0) * a_hex, 0.0])

    # z periodic
    a3 = np.array([0.0, 0.0, nlayer * d111])

    # A/B/C lateral shifts
    shift = (a1 + a2) / 3.0
    shifts = {"A": 0.0 * shift, "B": 1.0 * shift, "C": 2.0 * shift}

    positions, symbols = [], []

    # bottom -> top
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


def make_fcc111_twin_stacking_single_plane(d1: int, d2: int, d3: int, start: str = "A"):
    """
    bottom matrix (d3) | twin (d2) | top matrix (d1)
    CTB is a single atomic layer (shared plane), no duplicated layer.
    Total layers = d1 + d2 + d3 - 2
    """
    for name, d in [("d1", d1), ("d2", d2), ("d3", d3)]:
        if not isinstance(d, int) or d <= 0:
            raise ValueError(f"{name} must be a positive integer (layer count).")
    if d2 < 2:
        raise ValueError("d2 must be >= 2 (twin lamella must span two CTBs).")

    start = start.upper()
    if start not in "ABC":
        raise ValueError("start must be 'A'/'B'/'C'.")

    fwd = {"A": "B", "B": "C", "C": "A"}  # matrix
    bwd = {"A": "C", "C": "B", "B": "A"}  # twin (reversed)

    def grow(n, next_map, first_letter):
        out = [first_letter]
        cur = first_letter
        for _ in range(n - 1):
            cur = next_map[cur]
            out.append(cur)
        return out

    seq = []

    bottom = grow(d3, fwd, start)
    seq += bottom

    # twin: first layer shared with last of bottom
    twin_first = seq[-1]
    twin = grow(d2, bwd, twin_first)
    seq += twin[1:]

    # top matrix: first layer shared with last of twin
    top_first = seq[-1]
    top = grow(d1, fwd, top_first)
    seq += top[1:]

    return "".join(seq)


def build_cu111_twin_single_plane(d1, d2, d3, nx=4, ny=4, a=3.615, start="A"):
    stacking = make_fcc111_twin_stacking_single_plane(d1=d1, d2=d2, d3=d3, start=start)
    atoms = build_cu111_from_stacking(stacking=stacking, nx=nx, ny=ny, a=a)
    return atoms, stacking


def make_rectangular_cell_via_supercell(atoms: Atoms):
    """
    将(111)面内60°六方晶格 -> 面内90°矩形晶胞（正交超胞）
    关键：必须用 make_supercell 同时变换 cell 和原子，保证周期严格一致。

    使用整数变换矩阵 P，使新基矢：
      A1 =  a1 + a2
      A2 = -a1 + a2
      A3 =  a3
    其中 A1 ⟂ A2（面内90°）。
    """
    P = np.array([
        [ 1,  1, 0],
        [-1,  1, 0],
        [ 0,  0, 1],
    ], dtype=int)

    sup = make_supercell(atoms, P)
    sup.wrap(eps=1e-12)
    return sup


if __name__ == "__main__":
    # ===== 参数区 =====
    d1 = 7
    d2 = 12
    d3 = 5

    nx, ny = 5, 5
    a = 3.615
    start = "A"

    # 是否生成矩形边界条件（面内90°，三边90°）
    rect_cell = True

    # 输出到新文件夹，只保留 .vasp
    out_dir = "正交Cu111_single"
    os.makedirs(out_dir, exist_ok=True)

    atoms, stacking = build_cu111_twin_single_plane(
        d1=d1, d2=d2, d3=d3,
        nx=nx, ny=ny,
        a=a, start=start
    )

    if rect_cell:
        atoms = make_rectangular_cell_via_supercell(atoms)

    vasp_path = os.path.join(out_dir, f"POSCAR_Cu111_twin_singleplane_d1_{d1}_d2_{d2}_d3_{d3}.vasp")
    write(vasp_path, atoms, format="vasp")

    print("Done:")
    print(f"  rect_cell = {rect_cell}")
    print(f"  d1/d2/d3  = {d1}/{d2}/{d3}")
    print(f"  nx/ny     = {nx}/{ny}")
    print(f"  Stacking  = {stacking}")
    print(f"  Layers    = {len(stacking)} (should be d1+d2+d3-2 = {d1+d2+d3-2})")
    print(f"  Atoms     = {len(atoms)}")
    print(f"  Output    = {vasp_path}")
