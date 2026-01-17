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


# ===================== 等价 + 跨层去重核心 =====================

_MAP = {"A": 0, "B": 1, "C": 2}
_INV = {0: "A", 1: "B", 2: "C"}


def _rotate(seq, r):
    n = len(seq)
    r %= n
    return seq[r:] + seq[:r]


def primitive_root_string(s: str) -> str:
    """
    如果 s 是某个更短字符串重复得到的，返回最短周期；否则返回 s 本身。
    例：ABABAB -> AB；ABCABC -> ABC
    """
    n = len(s)
    for p in range(1, n + 1):
        if n % p != 0:
            continue
        if s == s[:p] * (n // p):
            return s[:p]
    return s


def all_equivalent_transforms(stacking: str):
    """
    生成同层长度下的所有等价变换：
      - z 循环平移
      - A/B/C 全局循环换名
      - 反向（倒序） + 上述两种
    """
    s = stacking.strip().upper()
    n = len(s)
    seq = [_MAP[ch] for ch in s]

    out = []

    # original
    for r in range(n):
        rot = _rotate(seq, r)
        for k in (0, 1, 2):
            cand = [ (x + k) % 3 for x in rot ]
            out.append("".join(_INV[x] for x in cand))

    # reversed
    rev = list(reversed(seq))
    for r in range(n):
        rot = _rotate(rev, r)
        for k in (0, 1, 2):
            cand = [ (x + k) % 3 for x in rot ]
            out.append("".join(_INV[x] for x in cand))

    return out


def canonical_same_length(stacking: str) -> str:
    """同一层数长度下的 canonical（不改变长度）。"""
    return min(all_equivalent_transforms(stacking))


def canonical_reduced_key(stacking: str) -> str:
    """
    跨层 canonical key：
      1) 先枚举所有同长度等价变换
      2) 每个变换做 primitive 约化（允许长度变短）
      3) 对约化后的短串再做同长度 canonical
      4) 取全局最小作为 key

    这样 ABAB、ABABAB、ABABABABAB 都会归并到 key=AB（或其 canonical 等价）。
    """
    best = None
    for t in all_equivalent_transforms(stacking):
        root = primitive_root_string(t)
        can_root = canonical_same_length(root)
        if best is None or can_root < best:
            best = can_root
    return best


def enumerate_unique_stackings_up_to(nlayer_max: int):
    """
    遍历 2..nlayer_max 所有 stackings，跨层去重：
    - 若某结构等价于更短周期（例如 ABAB -> AB），只保留最短的那个（自然会在更小层数先被加入）
    返回：
      kept: dict[key] = representative_stacking (代表结构，层数最少)
      per_layer: dict[nlayer] = list of kept stackings with this layer count (方便输出)
    """
    if nlayer_max < 2:
        raise ValueError("nlayer_max must be >= 2")

    kept = {}          # key -> representative
    per_layer = {}     # nlayer -> [representatives added at this nlayer]

    for n in range(2, nlayer_max + 1):
        per_layer[n] = []
        stackings = enumerate_stackings(n)

        for s in stackings:
            key = canonical_reduced_key(s)

            # 如果这个 key 已经由更短层数结构占领了，就跳过（保留最短）
            if key in kept:
                continue

            # 否则收下：此时 n 一定是该 key 的最小可达层数（因为我们从小到大遍历）
            kept[key] = s
            per_layer[n].append(s)

    return kept, per_layer


# =======================================================

if __name__ == "__main__":
    # ======= 你只需要改这里 =======
    nlayer_max = 12    # 会遍历 2..nlayer_max
    nx, ny = 4, 4
    # =============================

    out_dir = f"output_2to{nlayer_max}layers_unique"
    os.makedirs(out_dir, exist_ok=True)

    kept, per_layer = enumerate_unique_stackings_up_to(nlayer_max)

    total_kept = sum(len(v) for v in per_layer.values())
    print(f"Traverse layers: 2..{nlayer_max}")
    print(f"Total unique after cross-layer dedup: {total_kept}")

    # 写文件：按层数分子目录更清晰
    file_count = 0
    for n in range(2, nlayer_max + 1):
        reps = per_layer[n]
        if not reps:
            continue

        layer_dir = os.path.join(out_dir, f"{n}layers")
        os.makedirs(layer_dir, exist_ok=True)

        for idx, stacking in enumerate(sorted(reps), start=1):
            atoms = build_cu111_from_stacking(
                stacking=stacking,
                nx=nx,
                ny=ny,
                a=3.615,
            )

            out_path = os.path.join(
                layer_dir,
                f"Cu_111_{n}L_{idx:05d}_{stacking}.vasp"
            )
            write(out_path, atoms, format="vasp")
            file_count += 1

        print(f"  {n} layers kept: {len(reps)}")

    print("Done.")
    print(f"Wrote {file_count} unique .vasp files into: {out_dir}/")
