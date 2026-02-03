#!/usr/bin/env python3
from ase.io import iread, write
import os
import sys
import argparse
from collections import Counter, defaultdict

def main():
    parser = argparse.ArgumentParser(description="Convert XYZ trajectory to individual VASP POSCAR files with statistics.")
    parser.add_argument("infile", help="Input XYZ file")
    parser.add_argument("--outdir", default="xyz2poscar", help="Output directory (default: xyz2poscar)")
    
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    formula_counts = Counter()
    formula_written = defaultdict(int)
    total_structures = 0

    print(f"Reading from {args.infile}...")
    
    # Use iread for memory efficiency
    for at in iread(args.infile):
        # Determine formula (metal mode preferred for consistency with original script)
        f = at.get_chemical_formula(mode="metal")
        
        formula_counts[f] += 1
        formula_written[f] += 1
        total_structures += 1
        
        idx = formula_written[f]
        filename = f"{f}_{idx:04d}.vasp"
        path = os.path.join(args.outdir, filename)
        
        try:
            write(path, at, format="vasp")
        except Exception as e:
            print(f"Error writing {path}: {e}")

    # Write summary
    summary_path = os.path.join(args.outdir, "summary.txt")
    unique_formulas = sorted(formula_counts.keys())

    with open(summary_path, "w", encoding="utf-8") as fp:
        fp.write("=== Chemical Formulas List ===\n")
        fp.write(", ".join(unique_formulas) + "\n\n")

        fp.write("=== Structure Counts per Formula ===\n")
        for f in unique_formulas:
            fp.write(f"{f}\t{formula_counts[f]}\n")

        fp.write("\n=== Total Structures ===\n")
        fp.write(str(total_structures) + "\n")

    print(f"Done. Processed {total_structures} structures into '{args.outdir}/'. Summary at '{summary_path}'.")

if __name__ == "__main__":
    main()
