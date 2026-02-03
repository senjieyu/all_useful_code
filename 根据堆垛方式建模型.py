import numpy as np
from ase import Atoms
from ase.io import write


def build_cu111_from_stacking(
    stacking: str,
    nx: int = 4,
    ny: int = 4,
    a: float = 3.615,   # Cu fcc lattice constant (Å)
):
    """
    Build periodic Cu close-packed (111) structure with arbitrary stacking sequence.
    - stacking: e.g. "ACBABC", "BACBA" ... (length = number of layers)
    - nx, ny: in-plane repeats of the (111) primitive net
    - No vacuum: periodic in x, y, z
    Returns: ASE Atoms
    """

    stacking = stacking.strip().upper()
    if len(stacking) == 0 or any(ch not in "ABC" for ch in stacking):
        raise ValueError("stacking must be a non-empty string containing only A/B/C, e.g. 'ACBABC' or 'BACBA'")

    nlayer = len(stacking)

    # (111) in-plane hex lattice constant and interlayer spacing
    a_hex = a / np.sqrt(2.0)
    d111 = a / np.sqrt(3.0)

    # 2D primitive vectors of (111) hex net
    a1 = np.array([a_hex, 0.0, 0.0])
    a2 = np.array([0.5 * a_hex, 0.5 * np.sqrt(3.0) * a_hex, 0.0])

    # Stack direction (normal to plane)
    a3 = np.array([0.0, 0.0, nlayer * d111])  # periodic along z, no vacuum

    # A/B/C lateral shifts
    shift = (a1 + a2) / 3.0
    shifts = {
        "A": 0.0 * shift,
        "B": 1.0 * shift,
        "C": 2.0 * shift,
    }

    positions = []
    symbols = []

    # Build atoms layer by layer: bottom -> top
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
        pbc=(True, True, True),  # no vacuum: periodic in z
    )

    # Wrap into cell
    atoms.wrap(eps=1e-12)
    return atoms


if __name__ == "__main__":
    # 你只改这一个：想要什么堆垛就写什么
    stacking = "ABCACBACBC"   # 比如 "ACBABC", "BACBA", "ABCABC", "ABAB" (注意只能用A/B/C)

    # 面内尺寸
    nx, ny = 3, 3

    atoms = build_cu111_from_stacking(stacking=stacking, nx=nx, ny=ny, a=3.615)

    # 输出
    write(f"Cu_111_{stacking}.xyz", atoms)
    write(f"POSCAR_Cu_111_{stacking}", atoms, format="vasp")

    print("Done:")
    print(f"  Layers: {len(stacking)}  |  Stacking: {stacking}  |  Atoms: {len(atoms)}")
    print(f"  Wrote: Cu_111_{stacking}.xyz, POSCAR_Cu_111_{stacking}")
