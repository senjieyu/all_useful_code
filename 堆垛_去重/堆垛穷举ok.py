import os
import numpy as np
from ase import Atoms
from ase.io import write


def build_cu111_from_stacking(
    stacking: str,
    nx: int = 4,
    ny: int = 4,
    a: float = 3.615,   # Cu fcc lattice constant (Å)
):
    stacking = stacking.strip().upper()
    if len(stacking) == 0 or any(ch not in "ABC" for ch in stacking):
        raise ValueError("stacking must be non-empty and contain only A/B/C")

    nlayer = len(stacking)

    a_hex = a / np.sqrt(2.0)
    d111 = a / np.sqrt(3.0)

    a1 = np.array([a_hex, 0.0, 0.0])
    a2 = np.array([0.5 * a_hex, 0.5 * np.sqrt(3.0) * a_hex, 0.0])
    a3 = np.array([0.0, 0.0, nlayer * d111])

    shift = (a1 + a2) / 3.0
    shifts = {
        "A": 0.0 * shift,
        "B": 1.0 * shift,
        "C": 2.0 * shift,
    }

    positions = []
    symbols = []

    for k, letter in enumerate(stacking):
        lateral = shifts[letter]
        z = k * d111
        for i in range(nx):
            for j in range(ny):
                positions.append(i * a1 + j * a2 + lateral + [0.0, 0.0, z])
                symbols.append("Cu")

    atoms = Atoms(
        symbols=symbols,
        positions=np.array(positions),
        cell=np.vstack([nx * a1, ny * a2, a3]),
        pbc=(True, True, True),
    )
    atoms.wrap(eps=1e-12)
    return atoms


def enumerate_stackings(nlayer: int):
    """
    Enumerate all stackings of length = nlayer with:
      1) adjacent layers different
      2) first layer different from last layer
    """
    if nlayer < 2:
        raise ValueError("nlayer must be >= 2")

    letters = ("A", "B", "C")
    results = []

    def dfs(prefix):
        k = len(prefix)
        if k == nlayer:
            if prefix[0] != prefix[-1]:
                results.append("".join(prefix))
            return

        for ch in letters:
            if k > 0 and ch == prefix[-1]:
                continue
            prefix.append(ch)
            dfs(prefix)
            prefix.pop()

    dfs([])
    return results


# ===================== 新增：等价去重 =====================

_MAP = {"A": 0, "B": 1, "C": 2}
_INV = {0: "A", 1: "B", 2: "C"}


def _rotate(seq, r):
    n = len(seq)
    r %= n
    return seq[r:] + seq[:r]


def canonical_stacking(stacking: str) -> str:
    """
    Return a canonical representative under equivalence operations:
      1) cyclic shift along z (periodic)
      2) global relabeling A/B/C -> (A,B,C)+k mod 3 (in-plane translation)
      3) reversal (mirror along z), plus cyclic shifts and relabeling

    The canonical form is the lexicographically smallest string among all transforms.
    """
    s = stacking.strip().upper()
    n = len(s)
    seq = [_MAP[ch] for ch in s]

    candidates = []

    # original direction
    for r in range(n):
        rot = _rotate(seq, r)
        for k in (0, 1, 2):  # global relabel
            cand = [ (x + k) % 3 for x in rot ]
            candidates.append("".join(_INV[x] for x in cand))

    # reversed direction
    rev = list(reversed(seq))
    for r in range(n):
        rot = _rotate(rev, r)
        for k in (0, 1, 2):
            cand = [ (x + k) % 3 for x in rot ]
            candidates.append("".join(_INV[x] for x in cand))

    return min(candidates)


def deduplicate_stackings(stackings):
    """
    Deduplicate by canonical form; keep one representative per equivalence class.
    Return: (unique_stackings, mapping)
      - unique_stackings: list of kept stackings
      - mapping: dict canonical -> representative
    """
    rep_by_can = {}
    for s in stackings:
        can = canonical_stacking(s)
        if can not in rep_by_can:
            rep_by_can[can] = s  # keep the first encountered as representative

    # stable order: sort by canonical (or keep insertion order;这里用 canonical 排序更“规整”)
    unique = [rep_by_can[k] for k in sorted(rep_by_can.keys())]
    return unique, rep_by_can


# =======================================================

if __name__ == "__main__":
    # ======= 你只需要改这里 =======
    nlayer = 4      # 层数（任意 >=2）
    nx, ny = 4, 4   # 面内尺寸
    # =============================

    out_dir = f"output_{nlayer}layers"
    os.makedirs(out_dir, exist_ok=True)

    stackings_all = enumerate_stackings(nlayer)
    stackings_unique, _ = deduplicate_stackings(stackings_all)

    print(f"Total valid stackings ({nlayer} layers): {len(stackings_all)}")
    print(f"Unique (after symmetry/periodic dedup): {len(stackings_unique)}")

    for idx, stacking in enumerate(stackings_unique, start=1):
        atoms = build_cu111_from_stacking(
            stacking=stacking,
            nx=nx,
            ny=ny,
            a=3.615,
        )

        out_path = os.path.join(
            out_dir,
            f"Cu_111_{nlayer}L_{idx:05d}_{stacking}.vasp"
        )

        write(out_path, atoms, format="vasp")

    print("Done.")
    print(f"Wrote {len(stackings_unique)} unique .vasp files into: {out_dir}/")
