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
    a3 = np.array([0.0, 0.0, nlayer * d111])  # periodic in z, no vacuum

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


def make_fcc111_twin_stacking_single_plane(
    d1: int, d2: int, d3: int,
    start: str = "A",
):
    """
    生成 fcc(111) 三段结构：bottom matrix (d3) | twin (d2) | top matrix (d1)
    且两条孪晶面(CTB)都是“单层原子面”：界面层采用“共享层”，不重复插层。

    典型界面堆垛：...ABC|CBA... （C 只出现一次）

    约定（重要）：
    - d3 = 底部母相“总层数”
    - d2 = 中间孪晶区“总层数”
    - d1 = 顶部母相“总层数”
    - 由于两条界面层分别被相邻两段共享，最终总层数 = d1 + d2 + d3 - 2
      （共享了 2 个界面单层）
    """
    for name, d in [("d1", d1), ("d2", d2), ("d3", d3)]:
        if not isinstance(d, int) or d <= 0:
            raise ValueError(f"{name} must be a positive integer (layer count).")
    if d2 < 2:
        raise ValueError("d2 must be >= 2 for a meaningful twin lamella with two CTBs.")

    start = start.upper()
    if start not in "ABC":
        raise ValueError("start must be 'A'/'B'/'C'.")

    fwd = {"A": "B", "B": "C", "C": "A"}  # matrix: A->B->C->A
    bwd = {"A": "C", "C": "B", "B": "A"}  # twin:   A->C->B->A

    def grow(n, next_map, first_letter):
        out = [first_letter]
        cur = first_letter
        for _ in range(n - 1):
            cur = next_map[cur]
            out.append(cur)
        return out

    seq = []

    # 1) bottom matrix: d3 layers
    bottom = grow(d3, fwd, start)
    seq += bottom

    # 2) twin lamella: d2 layers, first layer is SHARED with last of bottom
    # i.e. boundary is single plane: ...X | X ...
    twin_first = seq[-1]
    twin = grow(d2, bwd, twin_first)
    seq += twin[1:]  # skip first because it's shared

    # 3) top matrix: d1 layers, first layer is SHARED with last of twin
    top_first = seq[-1]
    top = grow(d1, fwd, top_first)
    seq += top[1:]  # skip first because it's shared

    return "".join(seq)


def build_cu111_twin_single_plane(d1, d2, d3, nx=4, ny=4, a=3.615, start="A"):
    stacking = make_fcc111_twin_stacking_single_plane(d1=d1, d2=d2, d3=d3, start=start)
    atoms = build_cu111_from_stacking(stacking=stacking, nx=nx, ny=ny, a=a)
    return atoms, stacking


if __name__ == "__main__":
    # 你调这三个：上/中/下 三段层数（注意：最终总层数 = d1 + d2 + d3 - 2）
    d1 = 4   # top matrix layers
    d2 = 10  # twin lamella layers (>=2)
    d3 = 5   # bottom matrix layers

    nx, ny = 5, 5
    a = 3.615

    atoms, stacking = build_cu111_twin_single_plane(d1=d1, d2=d2, d3=d3, nx=nx, ny=ny, a=a, start="A")

    # write(f"Cu111_twin_singleplane_d1{d1}_d2{d2}_d3{d3}.xyz", atoms)
    write(f"POSCAR_Cu111_twin_singleplane_d1{d1}_d2{d2}_d3{d3}", atoms, format="vasp")

    print("Done:")
    print(f"  d1/d2/d3 = {d1}/{d2}/{d3}")
    print(f"  Stacking = {stacking}")
    print(f"  Layers   = {len(stacking)}  (should be d1+d2+d3-2 = {d1+d2+d3-2})")
    print(f"  Atoms    = {len(atoms)}")
