#!/usr/bin/env python
import re
import os
import sys
import numpy as np
import argparse
from ase.data import chemical_symbols

def parse_lattice(lattice_str):
    """Parse Lattice=\"...\" into 3x3 matrix"""
    nums = [float(x) for x in lattice_str.strip().split()]
    if len(nums) != 9:
        raise ValueError(f"Lattice not 9 numbers: {lattice_str}")
    return np.array(nums).reshape(3, 3)

def parse_n_floats(s, n):
    """Parse n floats from string"""
    nums = [float(x) for x in s.strip().split()]
    if len(nums) != n:
        raise ValueError(f"Expected {n} numbers, got {len(nums)}: {s}")
    return nums

def build_type_map(all_elements):
    """Map elements to integer types based on atomic number (ASE order)"""
    elements = sorted(all_elements, key=lambda x: chemical_symbols.index(x))
    return {e: i for i, e in enumerate(elements)}

def convert(fin, fout=None):
    if fout is None:
        root, _ = os.path.splitext(fin)
        fout = root + ".cfg"

    # ---------- Pass 1: Collect all elements ----------
    all_elements = set()

    with open(fin, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue

            try:
                nat = int(line)
            except ValueError:
                continue # Skip garbage lines between frames if any

            header = f.readline()
            if not header:
                break

            for _ in range(nat):
                atom_line = f.readline()
                if not atom_line:
                    break
                sp = atom_line.split()
                if sp:
                    all_elements.add(sp[0])

    type_map = build_type_map(all_elements)
    print("Element to type map:", type_map)

    # ---------- Pass 2: Convert ----------
    with open(fin, "r") as f_in, open(fout, "w") as f_out:
        frame = 0
        while True:
            line = f_in.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue

            try:
                nat = int(line)
            except ValueError:
                # Robustness: try to read until next integer or EOF
                continue

            header = f_in.readline()
            if not header:
                break
            header = header.strip()

            # Lattice
            lat_match = re.search(r'Lattice="([^"]+)"', header)
            if not lat_match:
                print(f"Warning: Frame {frame} missing Lattice, skipping.")
                # consume atoms
                for _ in range(nat): f_in.readline()
                continue
                
            lattice = parse_lattice(lat_match.group(1))
            volume = abs(np.linalg.det(lattice))

            # Energy
            en_match = re.search(r'energy=([^\s]+)', header)
            energy = float(en_match.group(1)) if en_match else 0.0

            # Stress
            st_match = re.search(r'stress="([^"]+)"', header)
            if st_match:
                stress_mat = np.array(parse_n_floats(st_match.group(1), 9)).reshape(3, 3)
                virial_mat = -stress_mat * volume
            else:
                virial_mat = np.zeros((3,3))

            # Atoms
            species = []
            pos = []
            forces = []
            
            valid_frame = True
            for _ in range(nat):
                atom_line = f_in.readline()
                if not atom_line:
                    valid_frame = False
                    break
                sp = atom_line.split()
                s = sp[0]
                # Try to parse coords/forces. If missing forces, assume 0
                try:
                    x, y, z = map(float, sp[1:4])
                    if len(sp) >= 7:
                        fx, fy, fz = map(float, sp[4:7])
                    else:
                        fx, fy, fz = 0.0, 0.0, 0.0
                except ValueError:
                    valid_frame = False
                    break
                    
                species.append(s)
                pos.append([x, y, z])
                forces.append([fx, fy, fz])

            if not valid_frame:
                print(f"Error reading atoms in frame {frame}")
                break

            pos = np.array(pos)
            forces = np.array(forces)

            # Write CFG
            f_out.write("BEGIN_CFG\n")
            f_out.write(" Size\n")
            f_out.write(f"  {nat:6d} \n")
            f_out.write(" Supercell \n")
            f_out.write("{:15.10f} {:15.10f} {:15.10f}\n".format(*lattice[0]))
            f_out.write("{:15.10f} {:15.10f} {:15.10f}\n".format(*lattice[1]))
            f_out.write("{:15.10f} {:15.10f} {:15.10f}\n".format(*lattice[2]))
            f_out.write(
                "AtomData:  id type       cartes_x      cartes_y      cartes_z     fx          fy          fz\n"
            )
            for i in range(nat):
                f_out.write(
                    " {:6d} {:6d} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(
                        i + 1,
                        type_map[species[i]],
                        pos[i, 0],
                        pos[i, 1],
                        pos[i, 2],
                        forces[i, 0],
                        forces[i, 1],
                        forces[i, 2],
                    )
                )
            f_out.write("Energy \n")
            f_out.write(f"\t{energy} \n")
            f_out.write("PlusStress:  xx          yy          zz          yz          xz          xy \n")
            f_out.write(
                f"\t{virial_mat[0,0]}  \t{virial_mat[1,1]}  \t{virial_mat[2,2]}  "
                f"\t{virial_mat[1,2]}  \t{virial_mat[0,2]}  \t{virial_mat[0,1]} \n"
            )
            f_out.write("END_CFG \n")

            frame += 1

    print(f"Generated: {fout}, processed {frame} frames")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert extxyz to DeepMD .cfg format.")
    parser.add_argument("fin", help="Input .extxyz file")
    parser.add_argument("fout", nargs="?", help="Output .cfg file (optional)")
    args = parser.parse_args()

    convert(args.fin, args.fout)
