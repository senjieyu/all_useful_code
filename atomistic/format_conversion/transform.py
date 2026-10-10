#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unified Format Transformation Script
Supports conversion between XYZ, VASP (POSCAR/CONTCAR), LAMMPS Data, CIF, and DeepMD CFG formats.
"""

import os
import sys
import argparse
import glob
import re
import numpy as np
from pathlib import Path
from tqdm import tqdm

# Try importing ASE
try:
    from ase import Atoms
    from ase.io import read, write, iread
    from ase.io.lammpsdata import read_lammps_data
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.data import atomic_numbers, chemical_symbols
    HAS_ASE = True
except ImportError:
    HAS_ASE = False

# ==========================================
# Shared Utilities
# ==========================================

def get_files(input_path, patterns):
    """
    Get list of files from input_path (file or directory) matching patterns.
    """
    input_path = Path(input_path)
    if input_path.is_file():
        return [input_path]
    
    if not input_path.is_dir():
        return []
    
    files = []
    if isinstance(patterns, str):
        patterns = [patterns]
        
    for pattern in patterns:
        # Recursive search if pattern contains **
        if "**" in pattern:
             files.extend(list(input_path.glob(pattern)))
        else:
             # Search in directory (non-recursive by default unless glob does it)
             # But here we want to be flexible. Let's use walk for robust finding if needed,
             # or just glob. The previous scripts used os.walk sometimes.
             # Let's use recursive glob for simplicity and power.
             files.extend(list(input_path.rglob(pattern)))
             
    return sorted(list(set(files)))

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

# ==========================================
# Core Logic: VASP -> Data (Complex Triclinic)
# ==========================================

def _build_lammps_prism(lattice):
    a, b, c = lattice[0], lattice[1], lattice[2]
    lx = np.linalg.norm(a)
    e1 = a / lx
    xy = np.dot(b, e1)
    b_perp = b - xy * e1
    ly = np.linalg.norm(b_perp)
    if ly < 1e-12:
        e2 = np.array([0., 1., 0.]) if abs(e1[1]) < 0.9 else np.array([0., 0., 1.])
        e2 -= np.dot(e2, e1) * e1
        e2 /= np.linalg.norm(e2)
        ly = 0.0
    else:
        e2 = b_perp / ly
    e3 = np.cross(e1, e2)
    e3 /= np.linalg.norm(e3)
    xz = np.dot(c, e1)
    yz = np.dot(c, e2)
    lz = np.linalg.norm(c - xz*e1 - yz*e2)
    return lx, ly, lz, xy, xz, yz, e1, e2, e3

def convert_atoms_to_data(atoms, output_file, mass=None):
    try:
        # 1. Get Cell and Prism
        lattice = atoms.get_cell()
        lx, ly, lz, xy, xz, yz, e1, e2, e3 = _build_lammps_prism(lattice)
        
        # 2. Get Coords and project to prism basis
        positions = atoms.get_positions()
        # Projection: x = r . e1, y = r . e2, z = r . e3
        # We can do this vectorized:
        # R = [x, y, z]
        # X_lammps = R . [e1, e2, e3]^T
        basis_matrix = np.vstack([e1, e2, e3]).T
        projected_pos = positions @ basis_matrix
        
        # 3. Get Masses and Types
        # ASE atoms have atomic numbers. We need to map them to 1..N types.
        # We sort by atomic number to be deterministic.
        unique_numbers = sorted(set(atoms.get_atomic_numbers()))
        type_map = {z: i+1 for i, z in enumerate(unique_numbers)}
        
        # If mass is provided as float, use it for all.
        # Otherwise try to get from atoms.
        masses = {}
        if mass is not None:
            for t in range(1, len(unique_numbers) + 1):
                masses[t] = mass
        else:
            # Use ASE masses
            # Note: atomic_masses in ASE is a list indexed by Z
            # Or use atoms.get_masses() if available, but that is per atom.
            # We need per type.
            # Let's use the standard mass for the element Z.
            from ase.data import atomic_masses
            for z in unique_numbers:
                masses[type_map[z]] = atomic_masses[z]

        # 4. Write
        total_atoms = len(atoms)
        with open(output_file, "w") as f:
            f.write(f"Converted from ASE Atoms object (Triclinic)\n\n")
            f.write(f"{total_atoms} atoms\n")
            f.write(f"{len(unique_numbers)} atom types\n\n")
            
            f.write(f"0.0 {lx:.6f} xlo xhi\n0.0 {ly:.6f} ylo yhi\n0.0 {lz:.6f} zlo zhi\n")
            f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n\n")
            
            f.write("Masses\n\n")
            for t in sorted(masses.keys()):
                f.write(f"{t} {masses[t]:.4f}\n")
            f.write("\nAtoms # atomic\n\n")
            
            atomic_nums = atoms.get_atomic_numbers()
            for i, (z, pos) in enumerate(zip(atomic_nums, projected_pos)):
                atom_type = type_map[z]
                f.write(f"{i+1} {atom_type} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
        return True
    except Exception as e:
        print(f"Error converting atoms to data: {e}")
        return False

def handle_vasp2data(args):
    inputs = get_files(args.input, ["POSCAR*", "CONTCAR*", "*.vasp"])
    ensure_dir(args.output)
    
    for f in tqdm(inputs):
        try:
            # Use ASE to read, which handles element recognition automatically
            atoms = read(f)
            outname = Path(args.output) / (f.stem + ".data")
            convert_atoms_to_data(atoms, outname, args.mass)
        except Exception as e:
            print(f"Error {f}: {e}")

# ==========================================
# Core Logic: ExtXYZ -> CFG
# ==========================================

def convert_extxyz_to_cfg(fin, fout):
    # Simplified version of the logic
    try:
        from ase.data import chemical_symbols
        def parse_lattice(s): return np.array([float(x) for x in s.split()]).reshape(3,3)
        
        all_elements = set()
        # Pass 1: scan elements
        with open(fin, 'r') as f:
            while True:
                line = f.readline()
                if not line: break
                try: nat = int(line.strip())
                except: continue
                f.readline() # header
                for _ in range(nat):
                    line = f.readline()
                    if line: all_elements.add(line.split()[0])
                    
        elements = sorted(all_elements, key=lambda x: chemical_symbols.index(x))
        type_map = {e: i for i, e in enumerate(elements)}
        
        # Pass 2: convert
        with open(fin, 'r') as f_in, open(fout, 'w') as f_out:
            frame = 0
            while True:
                line = f_in.readline()
                if not line: break
                line = line.strip()
                if not line: continue
                try: nat = int(line)
                except: continue
                
                header = f_in.readline().strip()
                lat_match = re.search(r'Lattice="([^"]+)"', header)
                if not lat_match: 
                    for _ in range(nat): f_in.readline()
                    continue
                lattice = parse_lattice(lat_match.group(1))
                
                en_match = re.search(r'energy=([^\s]+)', header)
                energy = float(en_match.group(1)) if en_match else 0.0
                
                st_match = re.search(r'stress="([^"]+)"', header)
                virial = np.zeros(9)
                if st_match:
                    stress = np.array([float(x) for x in st_match.group(1).split()])
                    vol = abs(np.linalg.det(lattice))
                    virial = -stress * vol # approximate 9 components
                    # Note: stress usually 9 floats in extxyz
                
                atom_lines = []
                for _ in range(nat): atom_lines.append(f_in.readline())
                
                f_out.write("BEGIN_CFG\n Size\n")
                f_out.write(f"  {nat:6d} \n Supercell \n")
                for i in range(3): f_out.write(f"{lattice[i][0]:15.10f} {lattice[i][1]:15.10f} {lattice[i][2]:15.10f}\n")
                f_out.write("AtomData:  id type       cartes_x      cartes_y      cartes_z     fx          fy          fz\n")
                
                for i, al in enumerate(atom_lines):
                    parts = al.split()
                    s = parts[0]
                    xyz = [float(x) for x in parts[1:4]]
                    f = [0.,0.,0.]
                    if len(parts) >= 7: f = [float(x) for x in parts[4:7]]
                    f_out.write(f" {i+1:6d} {type_map.get(s, 0):6d} {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f} {f[0]:12.6f} {f[1]:12.6f} {f[2]:12.6f}\n")
                    
                f_out.write(f"Energy \n\t{energy} \n")
                f_out.write("PlusStress:  xx          yy          zz          yz          xz          xy \n")
                # Need to check order, typically xx yy zz yz xz xy for VASP/LAMMPS interaction
                # If stress was 3x3:
                if len(virial) == 9:
                    # 0 1 2
                    # 3 4 5
                    # 6 7 8
                    # xx=0, yy=4, zz=8, yz=5, xz=2, xy=1 (assumed symmetric)
                    v = virial
                    f_out.write(f"\t{v[0]} \t{v[4]} \t{v[8]} \t{v[5]} \t{v[2]} \t{v[1]} \n")
                else:
                    f_out.write("\t0.0 \t0.0 \t0.0 \t0.0 \t0.0 \t0.0 \n")
                f_out.write("END_CFG \n")
                frame += 1
        return True
    except Exception as e:
        print(f"Error converting to CFG: {e}")
        return False

# ==========================================
# Handlers
# ==========================================

def handle_xyz2poscar(args):
    ensure_dir(args.outdir)
    inputs = get_files(args.input, ["*.xyz", "*.extxyz", "*.dump*"])
    
    for infile in inputs:
        print(f"Processing {infile}...")
        try:
            # Check format or try auto-detect
            # If dump file, might need element mapping
            kwargs = {}
            if "dump" in infile.name.lower() and args.elements:
                # Build map if elements provided
                # Note: ASE read(format='lammps-dump-text') reads types as atomic numbers if not specified
                # But here we want to map types 1..N to elements
                pass 
                
            count = 0
            # Use iread for streaming
            # Note: iread might not work well with all dump formats in auto-detect mode, 
            # but standard XYZ/ExtXYZ is fine.
            # For dump files, we often need to specify format.
            fmt = None
            if "dump" in infile.name.lower():
                fmt = 'lammps-dump-text'

            # Determine element map if needed (mostly for dump)
            map_dic = {}
            if args.elements:
                for i, sym in enumerate(args.elements):
                    map_dic[i+1] = sym
            
            for atoms in iread(str(infile), index=":", format=fmt):
                # If element mapping is needed (e.g. dump file has types but no symbols)
                if map_dic:
                    try:
                        an = atoms.get_atomic_numbers()
                        # If atomic numbers are just 1, 2, 3... and we want to map them
                        # Check if they look like types (1..N)
                        if max(an) <= len(args.elements):
                            new_syms = [map_dic.get(z, 'X') for z in an]
                            atoms.set_chemical_symbols(new_syms)
                    except Exception:
                        pass

                # Naming strategy
                if args.naming == "simple":
                    # 0.vasp, 1.vasp
                    outname = f"{count}.vasp"
                else:
                    # {filename}_{formula}_{index}.vasp
                    formula = atoms.get_chemical_formula(mode="metal")
                    outname = f"{infile.stem}_{formula}_{count:04d}.vasp"
                    
                write(Path(args.outdir) / outname, atoms, format="vasp")
                count += 1
            print(f"  -> Wrote {count} frames to {args.outdir}")
        except Exception as e:
            print(f"Error processing {infile}: {e}")

def handle_dump2xyz(args):
    ensure_dir(Path(args.output).parent)
    inputs = get_files(args.input, ["*.dump*", "dump*"])
    
    # Element map
    map_dic = {}
    if args.elements:
        for i, sym in enumerate(args.elements):
            map_dic[i+1] = sym
    
    total_frames = 0
    # Determine output mode: single file (if output ends with .xyz) or directory
    out_path = Path(args.output)
    is_single_file = out_path.suffix.lower() in ['.xyz', '.extxyz']
    
    if is_single_file:
        if out_path.exists(): os.remove(out_path)
        print(f"Converting dump files to {out_path}...")
        
        for infile in inputs:
            try:
                for atoms in iread(str(infile), index=":", format='lammps-dump-text'):
                    if map_dic:
                        an = atoms.get_atomic_numbers()
                        new_syms = [map_dic.get(z, 'X') for z in an]
                        atoms.set_chemical_symbols(new_syms)
                    
                    write(out_path, atoms, format='extxyz', append=True)
                    total_frames += 1
            except Exception as e:
                print(f"Error converting {infile}: {e}")
    else:
        ensure_dir(out_path)
        print(f"Converting dump files to directory {out_path}...")
        for infile in inputs:
            try:
                outname = out_path / (infile.name + ".xyz")
                temp_frames = []
                for atoms in iread(str(infile), index=":", format='lammps-dump-text'):
                    if map_dic:
                        an = atoms.get_atomic_numbers()
                        new_syms = [map_dic.get(z, 'X') for z in an]
                        atoms.set_chemical_symbols(new_syms)
                    temp_frames.append(atoms)
                write(outname, temp_frames, format='extxyz')
                total_frames += len(temp_frames)
                print(f"  {infile.name} -> {outname.name}")
            except Exception as e:
                print(f"Error converting {infile}: {e}")
                
    print(f"Done. Processed {total_frames} frames.")

def handle_xyz2dump(args):
    ensure_dir(Path(args.output).parent)
    inputs = get_files(args.input, ["*.xyz", "*.extxyz"])
    
    out_path = Path(args.output)
    is_single_file = out_path.suffix.lower() in ['.dump', '.lammpstrj'] or not out_path.suffix
    
    # If output looks like a directory but user wants single file? 
    # Usually dump files are single large files.
    
    if len(inputs) == 1 and is_single_file:
        # One to One
        infile = inputs[0]
        print(f"Converting {infile} to {out_path}...")
        try:
            atoms_list = read(str(infile), index=":")
            write(out_path, atoms_list, format='lammps-dump-text')
        except Exception as e:
            print(f"Error: {e}")
    else:
        # Many to Many or Many to One
        if is_single_file:
            print(f"Merging {len(inputs)} files to {out_path}...")
            if out_path.exists(): os.remove(out_path)
            for infile in tqdm(inputs):
                try:
                    atoms_list = read(str(infile), index=":")
                    write(out_path, atoms_list, format='lammps-dump-text', append=True)
                except Exception as e:
                    print(f"Error {infile}: {e}")
        else:
            ensure_dir(out_path)
            print(f"Converting {len(inputs)} files to {out_path}...")
            for infile in tqdm(inputs):
                try:
                    atoms_list = read(str(infile), index=":")
                    outname = out_path / (infile.stem + ".dump")
                    write(outname, atoms_list, format='lammps-dump-text')
                except Exception as e:
                    print(f"Error {infile}: {e}")

def handle_dump2poscar(args):
    # Wrapper reusing xyz2poscar logic but strictly for dump files
    # args.input -> directory or file
    # args.outdir -> output directory
    # args.elements -> required
    args.naming = args.naming if hasattr(args, 'naming') else 'descriptive'
    handle_xyz2poscar(args)

def handle_poscar2dump(args):
    # Similar to poscar2xyz but output dump
    inputs = get_files(args.input, ["POSCAR*", "CONTCAR*", "*.vasp"])
    if not inputs:
        print("No input files found.")
        return

    output_path = Path(args.output)
    # If merged, single file. Else directory.
    
    if args.merge:
        if output_path.is_dir():
            outfile = output_path / "merged.dump"
        else:
            outfile = output_path
            ensure_dir(outfile.parent)
            
        if outfile.exists(): os.remove(outfile)
        print(f"Merging {len(inputs)} files into {outfile}...")
        
        for f in tqdm(inputs):
            try:
                atoms = read(f)
                write(outfile, atoms, format="lammps-dump-text", append=True)
            except Exception as e:
                print(f"Error {f}: {e}")
    else:
        ensure_dir(output_path)
        print(f"Converting {len(inputs)} files to directory {output_path}...")
        for f in tqdm(inputs):
            try:
                atoms = read(f)
                outname = output_path / (f.stem + ".dump")
                write(outname, atoms, format="lammps-dump-text")
            except Exception as e:
                print(f"Error {f}: {e}")


def handle_poscar2xyz(args):
    # If args.merge is True, output to a single file in outdir (or named outdir if it looks like a file)
    # If args.merge is False, output individual files to outdir
    
    inputs = get_files(args.input, ["POSCAR*", "CONTCAR*", "*.vasp"])
    if not inputs:
        print("No input files found.")
        return

    output_path = Path(args.output)
    
    if args.merge:
        # Determine output filename
        if output_path.suffix.lower() in ['.xyz', '.extxyz']:
            # User provided full filename
            outfile = output_path
            ensure_dir(outfile.parent)
        else:
            # User provided directory
            ensure_dir(output_path)
            outfile = output_path / "merged.xyz"
        
        if outfile.exists(): os.remove(outfile)
        
        print(f"Merging {len(inputs)} files into {outfile}...")
        count = 0
        for f in tqdm(inputs):
            try:
                atoms = read(f)
                atoms.info['origin_file'] = f.name
                write(outfile, atoms, format="extxyz", append=True)
                count += 1
            except Exception as e:
                print(f"Skipping {f}: {e}")
        print(f"Done. Wrote {count} frames.")
        
    else:
        ensure_dir(output_path)
        print(f"Converting {len(inputs)} files to {output_path}...")
        for f in tqdm(inputs):
            try:
                atoms = read(f)
                outname = output_path / (f.stem + ".xyz")
                write(outname, atoms, format="extxyz")
            except Exception as e:
                print(f"Error {f}: {e}")

def handle_data2vasp(args):
    if not HAS_ASE:
        print("Error: ASE not installed.")
        return
        
    inputs = get_files(args.input, ["*.data", "*.lmp"])
    ensure_dir(args.output)
    
    # Map elements
    Z_of_type = {}
    for i, sym in enumerate(args.elements):
        if sym in atomic_numbers:
            Z_of_type[i+1] = atomic_numbers[sym]
        else:
            print(f"Warning: Unknown element {sym}")
            
    print(f"Element mapping: {args.elements}")
    
    for f in tqdm(inputs):
        try:
            atoms = read_lammps_data(str(f), style="atomic", Z_of_type=Z_of_type, sort_by_id=True)
            outname = Path(args.output) / (f.stem + ".vasp")
            write(outname, atoms, format="vasp", vasp5=True, direct=True)
        except Exception as e:
            print(f"Error {f}: {e}")

def handle_cif2vasp(args):
    inputs = get_files(args.input, ["*.cif"])
    ensure_dir(args.output)
    
    for f in tqdm(inputs):
        try:
            atoms = read(f)
            outname = Path(args.output) / (f.stem + ".vasp")
            write(outname, atoms, format="vasp")
        except Exception as e:
            print(f"Error {f}: {e}")

def handle_xyz2data(args):
    ensure_dir(args.outdir)
    inputs = get_files(args.input, ["*.xyz", "*.extxyz"])
    
    for infile in inputs:
        print(f"Processing {infile}...")
        count = 0
        for i, atoms in enumerate(iread(str(infile), index=":")):
            if args.vacuum: atoms.center(vacuum=args.vacuum)
            if not args.no_wrap and not args.vacuum: atoms.wrap()
            
            outname = Path(args.outdir) / f"{infile.stem}_{i:04d}.data"
            write(outname, atoms, format="lammps-data", atom_style="atomic")
            count += 1
        print(f"  -> Wrote {count} files.")

def handle_extxyz2cfg(args):
    convert_extxyz_to_cfg(args.input, args.output)

def handle_merge(args):
    # Merge XYZ files
    inputs = get_files(args.input, ["*.xyz", "*.extxyz"])
    # Filter out output file if it is in the list
    out_path = Path(args.output).resolve()
    inputs = [f for f in inputs if f.resolve() != out_path]
    
    if not inputs:
        print("No files to merge.")
        return

    print(f"Merging {len(inputs)} files to {out_path}...")
    with open(out_path, "wb") as outfile:
        for f in tqdm(inputs):
            with open(f, "rb") as infile:
                data = infile.read()
                outfile.write(data)
                if not data.endswith(b"\n"):
                    outfile.write(b"\n")
    print("Done.")

# ==========================================
# Main Dispatcher
# ==========================================

def main():
    parser = argparse.ArgumentParser(
        description="""
Powerful Format Conversion Tool for Atomistic Simulations.
Supports conversion between XYZ, VASP (POSCAR), LAMMPS Data/Dump, CIF, and DeepMD formats.

Examples:
  # XYZ -> POSCAR (split frames)
  python transform.py xyz2poscar traj.xyz out_dir

  # POSCAR -> XYZ (merge all)
  python transform.py poscar2xyz vasp_dir merged.xyz --merge

  # Dump -> XYZ (with element mapping)
  python transform.py dump2xyz dump.lammpstrj out.xyz --elements Cu Ni

  # VASP -> LAMMPS Data (auto mass)
  python transform.py vasp2data vasp_dir data_dir

  # LAMMPS Data -> VASP (specify elements for types 1, 2...)
  python transform.py data2vasp data_dir vasp_dir --elements Ti Al
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- xyz2poscar ---
    p_x2p = subparsers.add_parser("xyz2poscar", help="Convert XYZ/Dump file(s) to individual POSCARs")
    p_x2p.add_argument("input", help="Input XYZ/Dump file or directory")
    p_x2p.add_argument("outdir", help="Output directory")
    p_x2p.add_argument("--naming", choices=["descriptive", "simple"], default="descriptive", help="Naming style (default: descriptive)")
    p_x2p.add_argument("--elements", nargs="+", help="Elements for dump file mapping (if input is dump)")
    p_x2p.set_defaults(func=handle_xyz2poscar)

    # --- poscar2xyz ---
    p_p2x = subparsers.add_parser("poscar2xyz", help="Convert POSCARs to XYZ")
    p_p2x.add_argument("input", help="Input directory or file")
    p_p2x.add_argument("output", help="Output directory (or filename if --merge)")
    p_p2x.add_argument("--merge", action="store_true", help="Merge all outputs into one XYZ file")
    p_p2x.set_defaults(func=handle_poscar2xyz)
    
    # --- dump2xyz ---
    p_dump2xyz = subparsers.add_parser("dump2xyz", help="Convert LAMMPS Dump to XYZ")
    p_dump2xyz.add_argument("input", help="Input Dump file or directory")
    p_dump2xyz.add_argument("output", help="Output XYZ file or directory")
    p_dump2xyz.add_argument("--elements", nargs="+", required=True, help="Element mapping (e.g. Cu Ni)")
    p_dump2xyz.set_defaults(func=handle_dump2xyz)

    # --- xyz2dump ---
    p_xyz2dump = subparsers.add_parser("xyz2dump", help="Convert XYZ to LAMMPS Dump")
    p_xyz2dump.add_argument("input", help="Input XYZ file or directory")
    p_xyz2dump.add_argument("output", help="Output Dump file or directory")
    p_xyz2dump.set_defaults(func=handle_xyz2dump)

    # --- dump2poscar ---
    p_d2p = subparsers.add_parser("dump2poscar", help="Convert Dump to POSCARs (alias to xyz2poscar)")
    p_d2p.add_argument("input", help="Input Dump file or directory")
    p_d2p.add_argument("outdir", help="Output directory")
    p_d2p.add_argument("--naming", choices=["descriptive", "simple"], default="simple", help="Naming style (default: simple)")
    p_d2p.add_argument("--elements", nargs="+", required=True, help="Element mapping")
    p_d2p.set_defaults(func=handle_dump2poscar)

    # --- poscar2dump ---
    p_p2d = subparsers.add_parser("poscar2dump", help="Convert POSCARs to Dump")
    p_p2d.add_argument("input", help="Input directory or file")
    p_p2d.add_argument("output", help="Output directory (or filename if --merge)")
    p_p2d.add_argument("--merge", action="store_true", help="Merge all outputs into one Dump file")
    p_p2d.set_defaults(func=handle_poscar2dump)


    # --- data2vasp ---
    p_d2v = subparsers.add_parser("data2vasp", help="Convert LAMMPS data to POSCAR")
    p_d2v.add_argument("input", help="Input directory or file")
    p_d2v.add_argument("output", help="Output directory")
    p_d2v.add_argument("--elements", nargs="+", default=["Cu"], help="Element symbols for atom types (1 2 ...)")
    p_d2v.set_defaults(func=handle_data2vasp)

    # --- vasp2data ---
    p_v2d = subparsers.add_parser("vasp2data", help="Convert POSCAR to LAMMPS data")
    p_v2d.add_argument("input", help="Input directory or file")
    p_v2d.add_argument("output", help="Output directory")
    p_v2d.add_argument("--mass", type=float, default=None, help="Force atomic mass (optional, auto-detected if not set)")
    p_v2d.set_defaults(func=handle_vasp2data)

    # --- cif2vasp ---
    p_c2v = subparsers.add_parser("cif2vasp", help="Convert CIF to POSCAR")
    p_c2v.add_argument("input", help="Input directory or file")
    p_c2v.add_argument("output", help="Output directory")
    p_c2v.set_defaults(func=handle_cif2vasp)
    
    # --- xyz2data ---
    p_x2d = subparsers.add_parser("xyz2data", help="Split XYZ trajectory into LAMMPS data files")
    p_x2d.add_argument("input", help="Input XYZ file or directory")
    p_x2d.add_argument("outdir", help="Output directory")
    p_x2d.add_argument("--no-wrap", action="store_true", help="Do not wrap atoms")
    p_x2d.add_argument("--vacuum", type=float, help="Add vacuum and center")
    p_x2d.set_defaults(func=handle_xyz2data)

    # --- extxyz2cfg ---
    p_e2c = subparsers.add_parser("extxyz2cfg", help="Convert ExtXYZ to DeepMD CFG")
    p_e2c.add_argument("input", help="Input file")
    p_e2c.add_argument("output", help="Output file")
    p_e2c.set_defaults(func=handle_extxyz2cfg)
    
    # --- merge ---
    p_merge = subparsers.add_parser("merge", help="Merge multiple XYZ files")
    p_merge.add_argument("input", help="Input directory")
    p_merge.add_argument("output", help="Output filename")
    p_merge.set_defaults(func=handle_merge)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
