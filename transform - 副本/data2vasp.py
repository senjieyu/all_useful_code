#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch convert LAMMPS .data files -> VASP .vasp files.

Usage:
  python data2vasp.py <input_dir> <output_dir> --elements Cu Ni ...
"""

import argparse
import sys
import glob
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Batch convert LAMMPS .data files to VASP .vasp files.")
    parser.add_argument("input_dir", type=str, help="Input directory containing .data files")
    parser.add_argument("output_dir", type=str, help="Output directory for .vasp files")
    parser.add_argument("--elements", nargs="+", default=["Cu"], help="Element symbols corresponding to LAMMPS atom types 1, 2, ... (default: Cu)")
    
    args = parser.parse_args()

    input_dir = Path(args.input_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    
    if not input_dir.exists():
        print(f"Error: Input directory '{input_dir}' does not exist.")
        sys.exit(1)
        
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
    # Map input elements to Z_of_type
    Z_of_type = {}
    for i, symbol in enumerate(args.elements):
        if symbol not in atomic_numbers:
             print(f"Error: Unknown element symbol '{symbol}'")
             sys.exit(1)
        Z_of_type[i + 1] = atomic_numbers[symbol]
    
    print(f"Using element mapping: {dict(zip(Z_of_type.keys(), args.elements))}")

    pattern = str(input_dir / "*.data")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"No .data files found in: {input_dir}")
        sys.exit(1)

    ok, fail = 0, 0
    for f in files:
        fpath = Path(f)
        outpath = output_dir / (fpath.stem + ".vasp")

        try:
            atoms = read_lammps_data(
                str(fpath),
                style="atomic",
                Z_of_type=Z_of_type,
                sort_by_id=True
            )
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
