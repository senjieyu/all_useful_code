#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rand_dis_battery_data_prism_correct_v2.py

Purpose
-------
Read XYZ/ExtXYZ (multi-frame supported), generate EOS scaled structures, and write
LAMMPS triclinic (prism) data files using the robust rotation/projection method
from transform.py to avoid issues with Hexagonal lattices.

Changes vs Original
-------------------
- Replaced `build_lammps_prism_from_cell` with `_build_lammps_prism` (from transform.py).
- Updated coordinate calculation: strictly project Cartesian coords onto the
  LAMMPS prism basis vectors (e1, e2, e3) instead of estimating from fractional coords.
- Ensures 'a' vector is along x-axis, 'b' in xy-plane.
"""

import os
import sys
import argparse
import numpy as np
from ase.io import iread
from ase.build import make_supercell
from ase.data import atomic_numbers, atomic_masses


# ----------------------------
# Utilities from transform.py (Core Logic)
# ----------------------------

def _build_lammps_prism(lattice):
    """
    Build LAMMPS prism parameters and basis vectors from an arbitrary lattice.
    Logic copied from transform.py to ensure robust Hexagonal handling.

    Returns:
        lx, ly, lz, xy, xz, yz: Box dimensions and tilts
        e1, e2, e3: Normalized basis vectors aligned with LAMMPS requirements
    """
    a, b, c = lattice[0], lattice[1], lattice[2]

    # 1. Align a along x-axis (e1)
    lx = np.linalg.norm(a)
    e1 = a / lx

    # 2. Align b in xy-plane (e2)
    xy = np.dot(b, e1)
    b_perp = b - xy * e1
    ly = np.linalg.norm(b_perp)

    if ly < 1e-12:
        # Fallback for degenerate cells (shouldn't happen in real crystals)
        e2 = np.array([0., 1., 0.]) if abs(e1[1]) < 0.9 else np.array([0., 0., 1.])
        e2 -= np.dot(e2, e1) * e1
        e2 /= np.linalg.norm(e2)
        ly = 0.0
    else:
        e2 = b_perp / ly

    # 3. c vector determines the rest (e3 is simply cross product for orthogonality check)
    e3 = np.cross(e1, e2)
    e3 /= np.linalg.norm(e3)

    xz = np.dot(c, e1)
    yz = np.dot(c, e2)

    # lz is the height perpendicular to xy plane
    lz = np.linalg.norm(c - xz * e1 - yz * e2)

    return lx, ly, lz, xy, xz, yz, e1, e2, e3


# ----------------------------
# Existing Utilities
# ----------------------------
def parse_sc_matrix(sc_str: str) -> np.ndarray:
    rows = [r.strip() for r in sc_str.split(";") if r.strip()]
    if len(rows) != 3:
        raise ValueError("Supercell matrix must have 3 rows separated by ';'")
    mat = []
    for r in rows:
        parts = r.split()
        if len(parts) != 3:
            raise ValueError("Each row must have 3 numbers")
        mat.append([float(x) for x in parts])
    return np.array(mat, dtype=float)


def auto_specorder(atoms):
    """Stable element order by first appearance."""
    syms = atoms.get_chemical_symbols()
    seen = []
    for s in syms:
        if s not in seen:
            seen.append(s)
    return seen


def masses_from_specorder(specorder):
    masses = []
    for sym in specorder:
        znum = atomic_numbers.get(sym, None)
        masses.append(float(atomic_masses[znum]) if znum else 0.0)
    return masses


# ----------------------------
# Core writer (Updated logic)
# ----------------------------
def write_lammps_data_unified(
        outpath: str,
        atoms,
        specorder,
        wrap_atoms: bool = True,
        title: str = ""
):
    """
    Write LAMMPS data file using rotation/projection logic from transform.py.
    This ensures that even if the input cell is rotated (common in hex),
    the output aligns strictly with LAMMPS requirements (a // x).
    """

    # 1. Handle wrapping if requested
    # We do this before extracting the cell to ensure atoms are inside the conceptual unit cell
    if wrap_atoms:
        atoms.wrap()

    # 2. Get Cell and Prism parameters
    lattice = atoms.get_cell()
    lx, ly, lz, xy, xz, yz, e1, e2, e3 = _build_lammps_prism(lattice)

    # 3. Project positions to the prism basis
    # In LAMMPS: x = r . e1, y = r . e2, z = r . e3
    # This effectively rotates the atoms so they match the (lx, ly, lz) box definition.
    positions = atoms.get_positions()
    basis_matrix = np.vstack([e1, e2, e3]).T
    projected_pos = positions @ basis_matrix

    # 4. Handle element mapping (specorder)
    sym_list = atoms.get_chemical_symbols()
    sym_to_type = {s: i + 1 for i, s in enumerate(specorder)}

    try:
        types = [sym_to_type[s] for s in sym_list]
    except KeyError as e:
        raise ValueError(f"Element {e} not in specorder={specorder}. Use --specorder to include all elements.")

    masses = masses_from_specorder(specorder)
    natoms = len(atoms)
    ntypes = len(specorder)

    # 5. Write File
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    with open(outpath, "w", encoding="utf-8") as f:
        f.write((title if title else "LAMMPS data (Unified Transform)") + "\n\n")
        f.write(f"{natoms} atoms\n")
        f.write(f"{ntypes} atom types\n\n")

        # Bounds: transform.py uses 0.0 as origin, which is cleaner
        f.write(f"0.0 {lx:.10f} xlo xhi\n")
        f.write(f"0.0 {ly:.10f} ylo yhi\n")
        f.write(f"0.0 {lz:.10f} zlo zhi\n")
        f.write(f"{xy:.10f} {xz:.10f} {yz:.10f} xy xz yz\n\n")

        f.write("Masses\n\n")
        for t, m in enumerate(masses, start=1):
            f.write(f"{t} {m:.10f}  # {specorder[t - 1]}\n")
        f.write("\n")

        f.write("Atoms # atomic\n\n")
        # Format: id type x y z
        for i in range(natoms):
            f.write(
                f"{i + 1} {types[i]} {projected_pos[i, 0]:.10f} {projected_pos[i, 1]:.10f} {projected_pos[i, 2]:.10f}\n")


# ----------------------------
# Main workflow (EOS generator)
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate EOS-scaled LAMMPS data files (Robust Hexagonal Support)."
    )
    ap.add_argument("input_xyz", help="input xyz/extxyz (multi-frame ok)")
    ap.add_argument("out_dir", help="output base directory")

    ap.add_argument(
        "--vol-scale",
        default="0.98,1.00,1.02,1.04,1.06,1.08,1.10,1.12",
        help="comma-separated V/V0 list (default 0.98..1.12)"
    )
    ap.add_argument(
        "--supercell",
        default="1 0 0; 0 1 0; 0 0 1",
        help='3x3 supercell matrix like "2 0 0; 0 2 0; 0 0 2"'
    )
    ap.add_argument("--subdir", default="eos_data", help="subdir under struct_i")
    ap.add_argument("--prefix", default="data_V", help="output filename prefix")
    ap.add_argument("--suffix", default=".lmp", help="output filename suffix")

    ap.add_argument(
        "--specorder",
        default=None,
        help='element order, e.g. "Cu,Ni,O". Default: auto by first appearance.'
    )

    ap.add_argument("--no-wrap", action="store_true", help="disable atomic wrapping")
    # Note: --no-shift-bounds removed as new logic always defaults to 0.0 origin

    args = ap.parse_args()

    in_file = os.path.abspath(args.input_xyz)
    out_base = os.path.abspath(args.out_dir)

    try:
        sc = parse_sc_matrix(args.supercell)
    except Exception as e:
        print(f"[ERROR] bad --supercell: {e}")
        sys.exit(2)

    try:
        vol_scale = [float(x) for x in args.vol_scale.split(",") if x.strip()]
    except Exception as e:
        print(f"[ERROR] bad --vol-scale: {e}")
        sys.exit(2)

    n_written = 0
    n_frames = 0

    print(f"Reading from: {in_file}")

    for frame_i, atoms in enumerate(iread(in_file, index=":")):
        n_frames += 1

        struct_dir = os.path.join(out_base, f"struct_{frame_i}")
        eos_dir = os.path.join(struct_dir, args.subdir)

        # Make supercell first
        a0 = make_supercell(atoms, sc)

        # Determine specorder
        if args.specorder:
            specorder = [s.strip() for s in args.specorder.split(",") if s.strip()]
        else:
            specorder = auto_specorder(a0)

        # Get scaled positions to preserve relative geometry during scaling
        spos = a0.get_scaled_positions()
        cell0 = a0.get_cell().array

        for vs in vol_scale:
            ls = vs ** (1.0 / 3.0)  # linear scale factor
            a_scaled = a0.copy()

            # Apply scaling to cell
            new_cell = cell0 * ls
            a_scaled.set_cell(new_cell, scale_atoms=False)

            # Restore relative atom positions
            a_scaled.set_scaled_positions(spos)

            outname = os.path.join(eos_dir, f"{args.prefix}{vs:.4f}{args.suffix}")

            write_lammps_data_unified(
                outname,
                a_scaled,
                specorder=specorder,
                wrap_atoms=(not args.no_wrap),
                title=f"struct_{frame_i} EOS V/V0={vs:.4f}"
            )
            n_written += 1

    print(f"Frames processed: {n_frames}")
    print(f"Data files written: {n_written}")
    if n_frames > 0:
        example_path = os.path.join(out_base, "struct_0", args.subdir, f"{args.prefix}1.0000{args.suffix}")
        print(f"Example output: {example_path}")


if __name__ == "__main__":
    main()

