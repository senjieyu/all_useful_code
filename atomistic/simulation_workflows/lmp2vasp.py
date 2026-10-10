#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch convert LAMMPS data files (.lmp) -> VASP .vasp files.

Usage:
  python data2vasp.py <input_dir> <output_dir>

Example:
  python data2vasp.py ./kspacing_test10 ./vasp_out
"""

import sys
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        print("Usage: python data2vasp.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = Path(sys.argv[1]).expanduser().resolve()
    output_dir = Path(sys.argv[2]).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- dependency ----
    try:
        from ase.io.lammpsdata import read_lammps_data
        from ase.io import write
        from ase.data import atomic_numbers
    except ImportError:
        print("ERROR: ASE is not installed. Install with: pip install ase")
        sys.exit(1)

    # ---- element mapping ----
    # 默认：把所有 LAMMPS atom type 都当成 Cu
    # 多元素示例：{1: atomic_numbers["Cu"], 2: atomic_numbers["O"]}
    Z_of_type = {1: atomic_numbers["Cu"]}

    # ---- find input files ----
    # 只读 .lmp（如果你也想兼容 .data，可以把第二行打开）
    files = sorted(input_dir.glob("*.lmp"))
    # files = sorted(set(input_dir.glob("*.lmp")) | set(input_dir.glob("*.data")))

    if not files:
        print(f"No .lmp files found in: {input_dir}")
        sys.exit(1)

    ok, fail = 0, 0
    for fpath in files:
        outpath = output_dir / (fpath.stem + ".vasp")

        try:
            atoms = read_lammps_data(
                str(fpath),
                style="atomic",       # 适配 atom_style atomic
                Z_of_type=Z_of_type,
                sort_by_id=True
            )
            # 写 VASP POSCAR (vasp5=True 输出元素行；direct=True 输出分数坐标)
            write(str(outpath), atoms, format="vasp", vasp5=True, direct=True)

            ok += 1
            print(f"[OK]   {fpath.name} -> {outpath.name}")
        except Exception as e:
            fail += 1
            print(f"[FAIL] {fpath.name}: {e}")

    print(f"\nDone. Success: {ok}, Failed: {fail}")
    print(f"Output folder: {output_dir}")

if __name__ == "__main__":
    main()
