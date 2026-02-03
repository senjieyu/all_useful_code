#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
from tqdm import tqdm
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator

def find_poscar_like_files(root_dir):
    """Recursively find POSCAR/CONTCAR/.vasp files."""
    patterns = ("POSCAR", "CONTCAR")
    matched = []

    for dirpath, _, filenames in os.walk(root_dir):
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            if fn in patterns or fn.lower().endswith(".vasp"):
                matched.append(full)

    return sorted(matched)

def main():
    parser = argparse.ArgumentParser(description="Convert VASP POSCAR/CONTCAR/.vasp files to extended XYZ format.")
    parser.add_argument("input_dir", help="Input directory containing VASP files")
    parser.add_argument("out_xyz", help="Output .xyz filename")
    parser.add_argument("--zero-props", action="store_true", help="Set energy, forces, and stress to zero (default: False)")
    
    args = parser.parse_args()

    files = find_poscar_like_files(args.input_dir)
    print(f"Found {len(files)} VASP structure files")

    if not files:
        print("❌ No POSCAR/CONTCAR or *.vasp files found.")
        sys.exit(1)

    # Remove output file if it exists to avoid appending to old data
    if os.path.exists(args.out_xyz):
        os.remove(args.out_xyz)

    count = 0
    for idx, fp in tqdm(list(enumerate(files)), total=len(files)):
        try:
            atoms = read(fp)
            
            atoms.info["label"] = idx
            atoms.info["label_name"] = os.path.basename(fp)

            if args.zero_props:
                forces = np.zeros((len(atoms), 3))
                stress6 = np.zeros(6)
                atoms.calc = SinglePointCalculator(atoms, energy=0.0, forces=forces, stress=stress6)
                atoms.info["virial"] = np.zeros(9)

            # Write incrementally
            write(args.out_xyz, atoms, format="extxyz", append=True)
            count += 1
        except Exception as e:
            print(f"Error reading {fp}: {e}")

    print(f"✅ Wrote: {args.out_xyz} (frames={count})")

if __name__ == "__main__":
    main()
