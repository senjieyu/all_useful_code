#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vasp2data.py
Convert VASP POSCAR/CONTCAR -> LAMMPS data (triclinic prism)
- Supports VASP4/VASP5
- Supports Selective Dynamics (assigns atom types by flags if enabled)
- KEY FIX: rotate/project atomic Cartesian coordinates into the same prism basis
  used to define (lx,ly,lz,xy,xz,yz), so box/coords are consistent.

Usage:
  python vasp2data.py input.vasp output.data
  python vasp2data.py input_folder output_folder
"""

import sys
import os
import numpy as np
import argparse

def _parse_poscar(lines):
    if len(lines) < 8:
        raise ValueError("POSCAR too short")

    # scale
    try:
        scale = float(lines[1].strip())
    except ValueError as e:
        raise ValueError("Failed to read scale on line 2") from e

    # lattice (3x3)
    lattice = []
    for i in range(2, 5):
        lattice.append([float(x) for x in lines[i].split()])
    lattice = np.array(lattice, dtype=float) * scale  # rows: a,b,c

    # VASP4/5 detection for natoms line
    # VASP5 has element symbols on line 6, counts on line 7
    # VASP4 has counts on line 6
    def _is_int_list(s):
        try:
            _ = [int(x) for x in s.split()]
            return True
        except Exception:
            return False

    if _is_int_list(lines[5]):
        natoms_line_idx = 5
        start_idx = 6
    else:
        natoms_line_idx = 6
        start_idx = 7

    natoms_list = [int(x) for x in lines[natoms_line_idx].split()]
    total_atoms = int(sum(natoms_list))

    # Selective Dynamics?
    selective = False
    if start_idx < len(lines) and lines[start_idx].strip().upper().startswith("S"):
        selective = True
        start_idx += 1

    # coord type line
    if start_idx >= len(lines):
        raise ValueError("Missing coordinate type line (Direct/Cartesian)")
    coord_type = lines[start_idx].strip().lower()
    start_idx += 1

    return lattice, scale, natoms_list, total_atoms, selective, coord_type, start_idx


def _build_lammps_prism_from_lattice(lattice):
    """
    Build LAMMPS triclinic prism parameters and the orthonormal basis (e1,e2,e3)
    consistent with those parameters.

    lattice rows: a,b,c in Cartesian coordinates
    Returns: (lx,ly,lz,xy,xz,yz,e1,e2,e3)
    """
    a = lattice[0].astype(float)
    b = lattice[1].astype(float)
    c = lattice[2].astype(float)

    lx = np.linalg.norm(a)
    if lx < 1e-12:
        raise ValueError("Lattice vector a has near-zero length")

    e1 = a / lx  # x-axis along a

    # b projected onto e1 gives xy; remaining part defines y-axis
    xy = float(np.dot(b, e1))
    b_perp = b - xy * e1
    ly = float(np.linalg.norm(b_perp))
    if ly < 1e-12:
        # Degenerate (a || b). Still define some e2 orthogonal to e1.
        tmp = np.array([0.0, 1.0, 0.0])
        if np.linalg.norm(np.cross(tmp, e1)) < 1e-8:
            tmp = np.array([0.0, 0.0, 1.0])
        e2 = tmp - np.dot(tmp, e1) * e1
        e2 = e2 / (np.linalg.norm(e2) + 1e-12)
        ly = 0.0
    else:
        e2 = b_perp / ly

    # Right-handed z-axis
    e3 = np.cross(e1, e2)
    e3 = e3 / (np.linalg.norm(e3) + 1e-12)

    # c components in this basis => xz, yz, lz
    xz = float(np.dot(c, e1))
    yz = float(np.dot(c, e2))
    c_perp = c - xz * e1 - yz * e2
    lz = float(np.linalg.norm(c_perp))

    # lz can be 0 for 2D slabs; that's allowed but be aware in LAMMPS
    return lx, ly, lz, xy, xz, yz, e1, e2, e3


def _cart_from_direct(frac, lattice):
    # frac: (3,) , lattice rows are a,b,c
    # Direct coords: r = f1*a + f2*b + f3*c
    return frac @ lattice


def convert_poscar_to_data(poscar_file, output_file, mass=63.546):
    if not os.path.exists(poscar_file):
        print(f"Error: {poscar_file} not found.")
        return

    try:
        with open(poscar_file, "r") as f:
            lines = f.readlines()

        lattice, scale, natoms_list, total_atoms, selective, coord_type, start_idx = _parse_poscar(lines)

        # Build prism and basis
        lx, ly, lz, xy, xz, yz, e1, e2, e3 = _build_lammps_prism_from_lattice(lattice)

        atoms_data = []
        atom_id = 1
        current_line = start_idx

        # Read coords
        for n in natoms_list:
            for _ in range(n):
                if current_line >= len(lines):
                    break
                parts = lines[current_line].split()
                if len(parts) < 3:
                    current_line += 1
                    continue

                # type assignment based on selective flags
                atom_type = 1
                if selective and len(parts) >= 6:
                    flags = [p.upper() for p in parts[3:6]]
                    if flags == ["F", "F", "F"]:
                        atom_type = 2  # FFF
                    elif flags == ["F", "F", "T"]:
                        atom_type = 3  # FFT

                vals = np.array([float(parts[0]), float(parts[1]), float(parts[2])], dtype=float)

                # Convert to Cartesian in the original POSCAR frame
                if coord_type.startswith("d"):  # Direct
                    r_cart = _cart_from_direct(vals, lattice)
                else:  # Cartesian
                    # In VASP, if coord_type is Cartesian, positions are also scaled by 'scale'
                    r_cart = vals * scale

                # === KEY FIX: project into the same prism basis used for lx,ly,lz,xy,xz,yz ===
                x = float(np.dot(r_cart, e1))
                y = float(np.dot(r_cart, e2))
                z = float(np.dot(r_cart, e3))

                atoms_data.append((atom_id, atom_type, x, y, z))
                atom_id += 1
                current_line += 1

        # Write LAMMPS data
        with open(output_file, "w") as f:
            f.write(f"Converted from {os.path.basename(poscar_file)} (Triclinic)\n\n")
            f.write(f"{total_atoms} atoms\n")
            
            # Determine max atom type
            max_type = max([t for _, t, _, _, _ in atoms_data]) if atoms_data else 1
            f.write(f"{max_type} atom types\n\n")

            f.write(f"0.0 {lx:.6f} xlo xhi\n")
            f.write(f"0.0 {ly:.6f} ylo yhi\n")
            f.write(f"0.0 {lz:.6f} zlo zhi\n")
            f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n\n")

            f.write("Masses\n\n")
            for t in range(1, max_type + 1):
                f.write(f"{t} {mass}\n")
            f.write("\n")

            f.write("Atoms # atomic\n\n")
            for atom_id, atom_type, x, y, z in atoms_data:
                f.write(f"{atom_id} {atom_type} {x:.6f} {y:.6f} {z:.6f}\n")

        print(f"Converted: {os.path.basename(poscar_file)} -> {os.path.basename(output_file)}")

    except Exception as e:
        print(f"Error converting {os.path.basename(poscar_file)}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Convert VASP POSCAR to LAMMPS data.")
    parser.add_argument("input", help="Input file or folder")
    parser.add_argument("output", help="Output file or folder")
    parser.add_argument("--mass", type=float, default=63.546, help="Atom mass (default: 63.546 for Cu)")
    
    args = parser.parse_args()

    input_arg = args.input
    output_arg = args.output

    # Batch mode
    if os.path.isdir(input_arg):
        if not os.path.exists(output_arg):
            os.makedirs(output_arg)
            print(f"Created output directory: {output_arg}")

        print(f"Processing directory: {input_arg} -> {output_arg}")

        file_list = sorted(os.listdir(input_arg))
        count = 0
        for filename in file_list:
            src_file = os.path.join(input_arg, filename)
            if not os.path.isfile(src_file) or filename.startswith("."):
                continue

            base_name = os.path.splitext(filename)[0]
            dst_file = os.path.join(output_arg, base_name + ".data")

            convert_poscar_to_data(src_file, dst_file, args.mass)
            count += 1

        print(f"Done! Processed {count} files.")

    # Single file mode
    else:
        out_dir = os.path.dirname(output_arg)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        convert_poscar_to_data(input_arg, output_arg, args.mass)


if __name__ == "__main__":
    main()
