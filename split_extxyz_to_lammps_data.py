#!/usr/bin/env python3
import argparse
import os
from ase.io import read, write


def main():
    ap = argparse.ArgumentParser(
        description="Split extended XYZ (extxyz) into per-frame LAMMPS data files (atomic)."
    )
    ap.add_argument("infile", help="Input extxyz file, e.g. 109.xyz")
    ap.add_argument("-o", "--outdir", default="structures_data", help="Output directory")
    ap.add_argument("--no-wrap", action="store_true", help="Do not wrap atoms into cell")
    ap.add_argument(
        "--vacuum",
        type=float,
        default=None,
        help="If set, center each frame and add vacuum (Angstrom). Use for non-PBC clusters/surfaces.",
    )
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    frames = read(args.infile, index=":")
    print("Number of frames:", len(frames))

    for i, atoms in enumerate(frames):
        # Option: add vacuum & center (typically for non-PBC)
        if args.vacuum is not None:
            atoms.center(vacuum=args.vacuum)

        # Default: wrap into cell (good for PBC; fixes negative coords)
        if not args.no_wrap and args.vacuum is None:
            atoms.wrap()

        outname = os.path.join(args.outdir, f"str_{i:03d}.data")
        write(outname, atoms, format="lammps-data", atom_style="atomic")

    print(f"Done. Wrote {len(frames)} files into: {args.outdir}")


if __name__ == "__main__":
    main()

