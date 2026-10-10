#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    python rename_by_formula_ase.py <input_dir> <output_dir>

Description:
    Rename POSCAR-like .vasp files by chemical formula (parsed from symbols+counts),
    and copy to <output_dir>.

    This version does NOT use ase.formula (avoids AssertionError on some ASE versions).
    It also ignores any extra trailing blocks after coordinates.

    If multiple structures share the same formula:
        Hf10Al6.vasp, Hf10Al6_2.vasp, ...

Example:
    python rename_by_formula_ase.py success renamed_vasp
"""

import os
import sys
import shutil


def parse_symbols_counts_vasp5(path: str):
    """Parse VASP5 POSCAR header (symbols line + counts line)."""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    if len(lines) < 7:
        raise ValueError(f"Too short: only {len(lines)} non-empty lines")

    symbols = lines[5].split()
    counts_str = lines[6].split()

    if not symbols:
        raise ValueError("Empty symbols line (line 6)")

    try:
        counts = [int(x) for x in counts_str]
    except Exception as e:
        raise ValueError(f"Counts line (line 7) not all int: {counts_str}") from e

    if len(symbols) != len(counts):
        raise ValueError(f"Symbols/counts mismatch: {symbols} vs {counts}")

    return symbols, counts


def formula_from_symbols_counts(symbols, counts, mode="vasp"):
    """
    Build formula string.
    mode="vasp": keep VASP order (as in POSCAR symbols line)
    mode="alpha": sort by element symbol alphabetically
    """
    pairs = [(s, n) for s, n in zip(symbols, counts) if n != 0]

    if mode == "alpha":
        pairs.sort(key=lambda x: x[0])

    # concatenate like Hf10Al6 (no '1' shown)
    out = []
    for s, n in pairs:
        out.append(s)
        if n != 1:
            out.append(str(n))
    return "".join(out)


def main(input_dir, output_dir):
    if not os.path.isdir(input_dir):
        print(f"❌ 输入目录不存在: {input_dir}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    formula_count = {}
    ok = fail = 0

    for fname in sorted(os.listdir(input_dir)):
        if not fname.endswith(".vasp"):
            continue

        src = os.path.join(input_dir, fname)

        try:
            symbols, counts = parse_symbols_counts_vasp5(src)
            formula = formula_from_symbols_counts(symbols, counts, mode="vasp")  # 或 "alpha"
        except Exception as e:
            fail += 1
            print(f"❌ 读取失败: {fname}  错误: {type(e).__name__}: {repr(e)}")
            continue

        n = formula_count.get(formula, 0) + 1
        formula_count[formula] = n

        new_name = f"{formula}.vasp" if n == 1 else f"{formula}_{n}.vasp"
        dst = os.path.join(output_dir, new_name)

        shutil.copy2(src, dst)
        ok += 1
        print(f"✔ {fname} -> {new_name}")

    print(f"\n🎉 完成！成功: {ok} 失败: {fail} 输出目录：{output_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python rename_by_formula_ase.py <input_dir> <output_dir>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])

