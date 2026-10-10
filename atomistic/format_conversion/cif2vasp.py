#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch convert CIF files (.cif) -> VASP POSCAR files (.vasp).

Usage:
  python batch_cif2vasp.py <input_dir> <output_dir>
"""

import os
import sys
import argparse
from ase.io import read, write

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main():
    parser = argparse.ArgumentParser(description="Recursively convert CIF files to VASP POSCAR format.")
    parser.add_argument("input_dir", help="Input directory containing .cif files")
    parser.add_argument("output_dir", help="Output directory for .vasp files")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    if not os.path.isdir(input_dir):
        print(f"Error: input_dir not found -> {input_dir}")
        sys.exit(1)

    count = 0
    fail = 0

    for root, dirs, files in os.walk(input_dir):
        for fname in files:
            if fname.lower().endswith(".cif"):
                cif_path = os.path.join(root, fname)

                # Keep directory structure
                rel_path = os.path.relpath(root, input_dir)
                out_subdir = os.path.join(output_dir, rel_path)
                mkdir(out_subdir)

                # Output filename
                base = os.path.splitext(fname)[0]
                out_path = os.path.join(out_subdir, base + ".vasp")

                try:
                    atoms = read(cif_path)
                    write(out_path, atoms, format="vasp")
                    count += 1
                    print(f"[OK] {cif_path}  ->  {out_path}")
                except Exception as e:
                    fail += 1
                    print(f"[FAIL] {cif_path}  ({e})")

    print("\n========== Done ==========")
    print(f"Converted: {count}")
    print(f"Failed   : {fail}")
    print(f"Output dir: {output_dir}")

if __name__ == "__main__":
    main()
