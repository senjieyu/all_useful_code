# -*- coding: utf-8 -*-
# core_generator.py

import os
import shutil
import subprocess
import tempfile
import random

import math
import time

try:
    import numpy as np
    from scipy.spatial import cKDTree

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("[WARN] Numpy/Scipy not found. Merging atoms will be slower/less accurate.")


# ------------------------
# 基础工具
# ------------------------
def find_atomsk():
    atomsk = shutil.which("atomsk") or shutil.which("atomsk.exe")
    if atomsk is None:
        raise RuntimeError("找不到 atomsk，请确认已安装并加入 PATH")
    return atomsk


def ensure_dir(d):
    os.makedirs(d, exist_ok=True)


def _sanitize_miller(s: str) -> str:
    """
    让 Atomsk 0.13.1 能稳定识别 orient 的 Miller 向量
    """
    s = str(s).strip()
    s = s.replace(",", " ")
    while "  " in s:
        s = s.replace("  ", " ")

    if (s.startswith("[") and s.endswith("]")) or (s.startswith("<") and s.endswith(">")):
        br_l, br_r = s[0], s[-1]
        core = s[1:-1].strip().replace(" ", "")
        if core in ("000", "0", ""):
            raise ValueError(f"非法 Miller 向量：{s}")
        return f"{br_l}{core}{br_r}"

    parts = s.split()
    if len(parts) == 3 and all(p.lstrip("+-").isdigit() for p in parts):
        core = "".join(parts)
        if core in ("000",):
            raise ValueError(f"非法 Miller 向量：{s}")
        return f"[{core}]"

    return s


def run_cmd_stream(cmd, cwd=None, timeout_s=None, log_path=None):
    """
    实时输出子进程 stdout/stderr 到终端，并可写入 log 文件。
    """
    print("\n[RUN]", " ".join(cmd))
    log_f = open(log_path, "w", encoding="utf-8") if log_path else None

    try:
        p = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,  # 行缓冲
            universal_newlines=True,
        )

        for line in iter(p.stdout.readline, ""):
            print(line, end="")
            if log_f:
                log_f.write(line)
        p.stdout.close()

        ret = p.wait(timeout=timeout_s)
        if ret != 0:
            raise RuntimeError(f"命令失败（返回码 {ret}）：{' '.join(cmd)}")
    except subprocess.TimeoutExpired:
        p.kill()
        raise RuntimeError(f"命令超时（>{timeout_s}s）：{' '.join(cmd)}")
    finally:
        if log_f:
            log_f.close()


def write_poly_txt(path, cfg, box_dims=None, n_grains=None):
    box = box_dims if box_dims else cfg.POLY_BOX
    ngrains = n_grains if n_grains is not None else cfg.POLY_GRAINS
    mode = cfg.TWIN_MODE

    with open(path, "w", encoding="utf-8") as f:
        f.write(f"box {box[0]} {box[1]} {box[2]}\n")

    # 'NT' mode or 'GNT' mode logic for write_poly_txt
    if mode == 'NT' or mode == 'GNT':
        # "Mode 6 based NT": Random position + Parallel Twin Planes
        # 1. Random Position (like all_Equal)
        ngrains = n_grains if n_grains is not None else cfg.POLY_GRAINS

        # Check for random orientation switches (Global)
        is_random_x = getattr(cfg, 'NT_GRAIN_ORIENT_RANDOM_X', True)
        is_random_y = getattr(cfg, 'NT_GRAIN_ORIENT_RANDOM_Y', True)
        is_random_z = getattr(cfg, 'NT_GRAIN_ORIENT_RANDOM_Z', True)

        # For GNT, we might pass specific random switch via cfg (ZoneConfig)
        if mode == 'GNT' and isinstance(cfg, ZoneConfig):
            is_random_x = cfg.NT_GRAIN_ORIENT_RANDOM_X
            is_random_y = cfg.NT_GRAIN_ORIENT_RANDOM_Y
            is_random_z = cfg.NT_GRAIN_ORIENT_RANDOM_Z

        print(f"[INFO] Generating {ngrains} grains for {mode} mode.")
        print(f"       Random Orient -> X:{is_random_x}, Y:{is_random_y}, Z:{is_random_z}")

        for _ in range(ngrains):
            x = random.uniform(0, box[0])
            y = random.uniform(0, box[1])
            z = random.uniform(0, box[2])

            # Atomsk node format: x y z alpha beta gamma
            # alpha: rotation around Z
            # beta: rotation around new Y
            # gamma: rotation around new Z

            # Implementation logic:
            # 1. Start with Y=[111] (Twin Normal)
            # 2. To randomize "around Y" (stacking axis rotation), we rotate beta?
            #    No, in Z-Y-Z Euler:
            #    - alpha (Z) -> tilts the Y axis away from vertical?
            #    - beta (Y) -> rotates around the Y axis.
            #    - gamma (Z) -> rotates around the new Z.

            # Wait, standard Euler (Z-Y-Z):
            # 1. Rotate around Z by alpha
            # 2. Rotate around NEW Y by beta
            # 3. Rotate around NEW Z by gamma

            # If we want to keep Y axis vertical (parallel twin planes), we must NOT rotate around X or Z axes that would tilt Y.
            # So alpha=0, beta=random, gamma=0?
            # Let's test:
            # - alpha=0, beta=theta, gamma=0 -> Rotate around Y by theta. Y axis stays vertical. Correct.

            # If we want 3D random orientation:
            # We need full random Euler angles.

            if is_random_x and is_random_y and is_random_z:
                # Full 3D random
                alpha = random.uniform(0, 360)
                beta = random.uniform(0, 360)
                gamma = random.uniform(0, 360)
            else:
                # Partial random
                # This is tricky because Euler angles are coupled.
                # If we only want "Random Y rotation" (stacking axis rotation), we use beta.
                # If we add "Random X", it's not directly mapped to one Euler angle in Z-Y-Z.

                # Let's fallback to:
                # If ANY of X or Z is enabled, we assume the user wants the twin planes to TILT.
                # If only Y is enabled, we keep twin planes parallel.

                if is_random_x or is_random_z:
                    # Tilt allowed.
                    # We map the switches roughly to Euler angles to give "some" randomness.
                    alpha = random.uniform(0, 360) if is_random_z else 0.0
                    # Beta in Z-Y-Z tilts the Z axis?
                    # Atomsk documentation says: "The orientation of the grain is defined by 3 Euler angles (phi, theta, psi) in the Z-Y-Z convention."
                    # Actually, for cubic crystal:
                    # [100] [010] [001] aligned with box X Y Z.
                    # We built the seed such that Y=[111] (Twin Normal).
                    # So the crystal is ALREADY rotated in the seed.
                    # The node rotation is ADDITIONAL rotation.

                    # If we apply 0 0 0, the seed orientation is preserved (Y=[111]).

                    # If we want to rotate around Y (twin normal), we need to rotate around the global Y axis.
                    # In Z-Y-Z convention, rotating around Y is the middle angle (theta).
                    # So 0 theta 0 => rotates around Y.

                    # So:
                    # Random Y (Stacking Axis) -> Beta
                    # Random Z (Tilt) -> Alpha
                    # Random X (Tilt) -> Gamma? (Gamma is around new Z, so effectively another tilt or in-plane rotation depending on beta).

                    alpha = random.uniform(0, 360) if is_random_z else 0.0
                    beta = random.uniform(0, 360) if is_random_y else 0.0
                    gamma = random.uniform(0, 360) if is_random_x else 0.0
                else:
                    # Only Y (or none)
                    alpha = 0.0
                    beta = random.uniform(0, 360) if is_random_y else 0.0
                    gamma = 0.0

            f.write(f"node {x:.4f} {y:.4f} {z:.4f} {alpha:.2f} {beta:.2f} {gamma:.2f}\n")

    elif mode == 'all_Equal':
        # Explicit random orientation for all_Equal mode to ensure 3D rotation
        print(f"[INFO] Generating {ngrains} grains with EXPLICIT RANDOM 3D orientation.")
        for _ in range(ngrains):
            x = random.uniform(0, box[0])
            y = random.uniform(0, box[1])
            z = random.uniform(0, box[2])
            # Random Euler angles (degrees)
            alpha = random.uniform(0, 360)
            beta = random.uniform(0, 360)
            gamma = random.uniform(0, 360)
            f.write(f"node {x:.4f} {y:.4f} {z:.4f} {alpha:.2f} {beta:.2f} {gamma:.2f}\n")
        else:
            # Default random mode
            f.write(f"random {int(ngrains)}\n")


def parse_lammps_atoms_count(lmp_path: str) -> int:
    """
    快速读取 LAMMPS data 头部，解析 'XXXX atoms'
    """
    with open(lmp_path, "r", errors="ignore") as f:
        for _ in range(120):
            line = f.readline()
            if not line:
                break
            s = line.strip().split()
            if len(s) == 2 and s[1].lower() == "atoms":
                return int(s[0])
    raise RuntimeError(f"无法从 {lmp_path} 头部解析 atoms 数（文件可能未写完/格式不对）")


def remove_doubles_average(input_path, output_path, dist_thresh):
    """
    Reads a LAMMPS data file, finds overlapping atoms (< dist_thresh),
    merges them by averaging their positions, and writes a new LAMMPS data file.

    Using cKDTree for performance.
    """
    if not HAS_NUMPY:
        raise RuntimeError("Numpy/Scipy required for average merging.")

    # Explicit check for zero/negative threshold to disable merging
    if dist_thresh <= 1e-6:
        print(f"[INFO] Distance threshold is {dist_thresh} (<= 0). Skipping merge.")
        shutil.copy(input_path, output_path)
        return

    print(f"[INFO] Merging overlapping atoms (dist < {dist_thresh}) with averaging...")

    # 1. Read LAMMPS file (simplified reader for 'atomic' or similar style)
    # We assume standard header followed by 'Atoms' section.

    header_lines = []
    atoms_lines = []

    with open(input_path, 'r') as f:
        lines = f.readlines()

    atom_start = -1
    box_bounds = []

    for i, line in enumerate(lines):
        if "Atoms" in line:
            atom_start = i + 2
            break
        header_lines.append(line)
        if "xlo xhi" in line or "ylo yhi" in line or "zlo zhi" in line:
            parts = line.split()
            box_bounds.append((float(parts[0]), float(parts[1])))

    if atom_start == -1 or atom_start >= len(lines):
        # Fallback: Atomsk lmp output sometimes puts 'Atoms' earlier or has different structure?
        # Let's try to find 'Atoms' again, maybe case sensitive?
        # Usually it is 'Atoms' or 'Atoms # atomic'
        for i, line in enumerate(lines):
            if line.strip().startswith("Atoms"):
                atom_start = i + 2
                break
        if atom_start == -1:
            raise ValueError(f"Could not find Atoms section in LAMMPS file: {input_path}")

    # Parse atoms
    # We need: id, type, x, y, z.
    # Usually: id type x y z (5 cols) or id type q x y z (6 cols)
    # We will try to detect x,y,z columns.

    # Let's assume standard 5 column for now (from Atomsk lmp output)
    # If not, we scan the first line.

    raw_data = []  # store [id, type, x, y, z, original_line_suffix]

    first_line_parts = lines[atom_start].split()
    n_cols = len(first_line_parts)

    # Heuristic for x,y,z index
    # If 5 cols: id type x y z
    # If 6 cols (charge): id type q x y z
    if n_cols == 5:
        x_idx, y_idx, z_idx = 2, 3, 4
        type_idx = 1
    elif n_cols >= 6:
        # Assume last 3 are x y z? Or 2,3,4?
        # Atomsk often writes: id type q x y z nx ny nz
        # Safest is usually 2,3,4 for atomic, or check if column 2 is charge (float) or x (coord).
        # Let's assume 2,3,4 if 6 cols?
        # Actually for 'charge' style: id type q x y z
        # so indices 3,4,5.
        # How to distinguish?
        # Let's check the header for "atom types" line? No.
        # Let's assume if 5 cols -> 2,3,4. If >5, try 2,3,4 first, check bounds?
        # Let's stick to 2,3,4 unless we are sure.
        # Wait, if `atomsk ... lmp` produced it, it defaults to 'atomic' (5 cols) usually unless charges are present.
        # The previous steps generated 5 cols (id type x y z) or 4 cols CFG.
        # The merged.cfg was converted to temp_merged.lmp.
        x_idx, y_idx, z_idx = 2, 3, 4
        type_idx = 1
    else:
        # Fallback
        x_idx, y_idx, z_idx = 2, 3, 4
        type_idx = 1

    coords = []
    types = []
    ids = []

    for line in lines[atom_start:]:
        parts = line.split()
        if len(parts) < 3: continue
        try:
            x = float(parts[x_idx])
            y = float(parts[y_idx])
            z = float(parts[z_idx])
            tid = int(parts[type_idx])
            ids.append(parts[0])
            coords.append([x, y, z])
            types.append(tid)
        except:
            continue

    coords = np.array(coords)
    n_atoms = len(coords)

    # 2. Find duplicates using KDTree
    tree = cKDTree(coords)

    # query_pairs returns set of (i, j) with i < j and dist < r
    pairs = tree.query_pairs(r=dist_thresh)

    if not pairs:
        print("[INFO] No overlapping atoms found.")
        # Must copy to output path because the rest of the script expects output_path to exist
        shutil.copy(input_path, output_path)
        return

    print(f"[INFO] Found {len(pairs)} overlapping pairs. Merging...")

    # 3. Merge logic (Disjoint Set Union / Connected Components)
    # We want to group all connected atoms (A-B, B-C) into one cluster and average them.

    parent = list(range(n_atoms))

    def find(i):
        if parent[i] == i: return i
        parent[i] = find(parent[i])
        return parent[i]

    def union(i, j):
        root_i = find(i)
        root_j = find(j)
        if root_i != root_j:
            parent[root_j] = root_i

    for i, j in pairs:
        union(i, j)

    # Group by root
    clusters = {}
    for i in range(n_atoms):
        root = find(i)
        if root not in clusters:
            clusters[root] = []
        clusters[root].append(i)

    # 4. Create new atom list
    new_coords = []
    new_types = []

    for root, indices in clusters.items():
        if len(indices) == 1:
            idx = indices[0]
            new_coords.append(coords[idx])
            new_types.append(types[idx])
        else:
            # Average coordinates
            cluster_coords = coords[indices]
            # Simple average (center of mass assuming equal mass)
            avg_pos = np.mean(cluster_coords, axis=0)

            # Type? Majority vote or keep first?
            # Let's keep the type of the first one (usually they are same element in this script)
            t_counts = {}
            for idx in indices:
                t = types[idx]
                t_counts[t] = t_counts.get(t, 0) + 1
            # Get max count type
            main_type = max(t_counts, key=t_counts.get)

            new_coords.append(avg_pos)
            new_types.append(main_type)

    # 5. Write output
    with open(output_path, 'w') as f:
        # Header
        for line in header_lines:
            if "atoms" in line and len(line.split()) == 2:
                f.write(f"{len(new_coords)} atoms\n")
            else:
                f.write(line)

        f.write("\nAtoms\n\n")

        for i, (xyz, t) in enumerate(zip(new_coords, new_types)):
            # id type x y z
            f.write(f"{i + 1} {t} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")

    print(f"[INFO] Merged {n_atoms} -> {len(new_coords)} atoms.")


class ZoneConfig:
    """Helper to mock a Config object for a specific GNT Zone."""

    def __init__(self, parent_cfg, zone_dict):
        # 1. Inherit global settings
        self.ELEMENT = parent_cfg.ELEMENT
        self.LATTICE_A = parent_cfg.LATTICE_A
        self.TWIN_ORIENT = parent_cfg.TWIN_ORIENT
        self.STACK_AXIS = getattr(parent_cfg, "STACK_AXIS", "Y")
        self.SEED_DUP_X = parent_cfg.SEED_DUP_X
        self.SEED_DUP_Z = parent_cfg.SEED_DUP_Z

        # 2. Zone specific settings
        # We set TWIN_MODE to 'NT' to ensure write_poly_txt generates
        # nodes rotated ONLY around Y axis (keeping twins parallel).
        self.TWIN_MODE = 'NT'

        # 3. Map parameters for build_layered_seed
        # These are accessed based on the mode passed to build_layered_seed override

        # 'equal' mode
        self.EQUAL_THICKNESS = zone_dict.get('equal_thickness', 3)

        # 'random' mode
        self.RANDOM_MIN_LAYERS = zone_dict.get('min_layers', 2)
        self.RANDOM_MAX_LAYERS = zone_dict.get('max_layers', 8)
        self.RANDOM_LAYER_COUNT = zone_dict.get('layer_count', 10)

        # 'ratio' mode
        self.RATIO_TWIN_FRACTION = zone_dict.get('ratio_fraction', 0.3)
        self.RATIO_TOTAL_LAYERS = zone_dict.get('total_layers', 20)
        self.RATIO_MIN_THICKNESS = zone_dict.get('min_thickness', 2)

        # 'all_Equal' mode params
        self.ALL_EQUAL_MATRIX_THICKNESS = zone_dict.get('matrix_thickness', 5)
        self.ALL_EQUAL_TWIN_THICKNESS = zone_dict.get('twin_thickness', 1)
        self.ALL_EQUAL_LAYER_COUNT = zone_dict.get('layer_count', 10)

        # Random Orientation Switch
        self.NT_GRAIN_ORIENT_RANDOM_X = zone_dict.get('random_orient_X', True)
        self.NT_GRAIN_ORIENT_RANDOM_Y = zone_dict.get('random_orient_Y', True)
        self.NT_GRAIN_ORIENT_RANDOM_Z = zone_dict.get('random_orient_Z', True)

        # Backward compatibility map 'random_orient' to all if present
        if 'random_orient' in zone_dict:
            ro = zone_dict['random_orient']
            # If user explicitly set random_orient=True/False, apply to Y (classic)
            # or apply to all? Let's assume they mean "Randomize Orientation" broadly.
            # But wait, old default was Y-rotation only for NT.
            # Let's trust the specific keys first.
            if 'random_orient_Y' not in zone_dict:
                self.NT_GRAIN_ORIENT_RANDOM_Y = ro

        # Twin Enable Switch
        self.ENABLE_TWIN = zone_dict.get('enable_twin', True)


class Node:
    def __init__(self, x, y, z, a, b, g, zone_idx):
        self.x = x
        self.y = y
        self.z = z
        self.a = a
        self.b = b
        self.g = g
        self.zone_idx = zone_idx


def parse_cfg_and_filter(cfg_path, valid_zone_idx, all_nodes, default_box=None, target_box=None):
    """
    Parse Atomsk CFG file and keep atoms belonging to grains in valid_zone_idx.
    Returns a list of atom tuples: (atom_type, x, y, z) in ABSOLUTE Angstroms.
    If target_box is provided, atoms will be wrapped into [0, target_box].
    """
    valid_atoms = []

    if default_box is None:
        default_box = [1.0, 1.0, 1.0]

    with open(cfg_path, 'r') as f:
        lines = f.readlines()

    # Parse Header
    entry_count = 0
    aux_map = {}  # aux_name -> index offset from coords

    # Parse Box from H matrix
    h_matrix = {}  # (row, col) -> val

    header_end_idx = 0
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("Number of particles ="):
            n_particles = int(line.split('=')[1].strip())
        elif line.startswith("entry_count ="):
            entry_count = int(line.split('=')[1].strip())
        elif line.startswith("auxiliary["):
            # format: auxiliary[0] = grainID
            parts = line.split('=')
            idx = int(parts[0].split('[')[1].split(']')[0])
            name = parts[1].strip()
            aux_map[name] = idx
        elif line.startswith("H0("):
            # H0(1,1) = 40.0
            p1 = line.split('=')
            val = float(p1[1].strip())
            key = p1[0].strip()  # H0(1,1)
            # Parse indices
            r = int(key[3])
            c = int(key[5])
            h_matrix[(r, c)] = val
        elif ".NO_VELOCITY." in line:
            pass
        elif len(line.split()) == 1 and line.replace('.', '').isdigit():
            # Mass line
            pass

        # Check for start of data
        try:
            parts = line.split()
            if len(parts) == entry_count:
                # Found data start?
                # Atomsk usually puts data after all headers.
                # But we need to process data lines only after we parsed H matrix.
                pass
        except:
            pass

    # Determine Box Dimensions
    # If H matrix was not found, fallback to default_box
    if not h_matrix:
        lx, ly, lz = default_box
        print(f"[WARN] H-matrix not found in {os.path.basename(cfg_path)}. Using default box: {lx}, {ly}, {lz}")
    else:
        lx = h_matrix.get((1, 1), default_box[0])
        ly = h_matrix.get((2, 2), default_box[1])
        lz = h_matrix.get((3, 3), default_box[2])

    print(f"[DEBUG] CFG Parsing: {os.path.basename(cfg_path)}")
    print(f"        Box Dimensions used: lx={lx}, ly={ly}, lz={lz}")

    if 'grainID' not in aux_map:
        print(f"[WARN] No grainID found in {cfg_path}, skipping filtering.")
        return []

    grain_col_idx = 3 + aux_map['grainID']

    for line in lines:
        parts = line.split()
        if len(parts) != entry_count:
            continue

        # Check if valid data line (floats)
        try:
            # Parse reduced coords
            sx = float(parts[0])
            sy = float(parts[1])
            sz = float(parts[2])

            # Check grainID
            gid = int(float(parts[grain_col_idx]))

            # Atomsk GrainID is 1-based index corresponding to node order.
            node_idx = gid - 1

            if 0 <= node_idx < len(all_nodes):
                if all_nodes[node_idx].zone_idx == valid_zone_idx:
                    # Keep this atom
                    # Convert to Absolute
                    # Simple orthorhombic conversion:
                    x = sx * lx
                    y = sy * ly
                    z = sz * lz

                    # Apply full H matrix if needed?
                    # x = H11*sx + H12*sy + H13*sz
                    # Usually H12, H13 are 0.
                    # Let's assume standard rectangular box for GNT.

                    # Wrap if needed
                    if target_box:
                        x = x % target_box[0]
                        y = y % target_box[1]
                        z = z % target_box[2]

                    # Type? CFG doesn't explicitly store type ID usually, but Atomsk might?
                    # Atomsk lmp -> CFG usually loses type ID unless stored in aux?
                    # Or it assumes type 1?
                    # Let's assume type 1 (Cu) for now.
                    valid_atoms.append((1, x, y, z))
        except ValueError:
            continue

    if valid_atoms:
        xs = [a[1] for a in valid_atoms]
        ys = [a[2] for a in valid_atoms]
        zs = [a[3] for a in valid_atoms]
        print(
            f"        Parsed Coords Range: X[{min(xs):.2f}, {max(xs):.2f}], Y[{min(ys):.2f}, {max(ys):.2f}], Z[{min(zs):.2f}, {max(zs):.2f}]")
    else:
        print("        No valid atoms extracted.")

    return valid_atoms


def write_merged_lammps(output_path, atom_tuples, box_dims, element):
    """
    Write a LAMMPS data file directly from merged atoms.
    atom_tuples: list of (type, x, y, z)
    """
    with open(output_path, 'w') as f:
        f.write("LAMMPS data file generated by Trae Agent\n\n")
        f.write(f"{len(atom_tuples)} atoms\n")
        f.write("1 atom types\n\n")

        f.write(f"0.0 {box_dims[0]:.6f} xlo xhi\n")
        f.write(f"0.0 {box_dims[1]:.6f} ylo yhi\n")
        f.write(f"0.0 {box_dims[2]:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        # Dummy mass
        f.write("1 63.546\n\n")

        f.write("Atoms # atomic\n\n")

        for i, (t, x, y, z) in enumerate(atom_tuples):
            # id type x y z
            f.write(f"{i + 1} {t} {x:.6f} {y:.6f} {z:.6f}\n")


# ------------------------
# 孪晶层生成逻辑
# ------------------------
def build_simple_seed(atomsk, cfg, matrix_h, twin_h, tmpdir, timeout_s):
    """
    Generate a simple 1-period seed: [Matrix(matrix_h) + Twin(twin_h)]
    """
    element = cfg.ELEMENT
    a = cfg.LATTICE_A
    orient = cfg.TWIN_ORIENT
    stack_axis = getattr(cfg, "STACK_AXIS", "Y")

    dup_x = cfg.SEED_DUP_X
    dup_z = cfg.SEED_DUP_Z

    o1, o2, o3 = [_sanitize_miller(x) for x in orient]

    # 1. Base Units
    unit_m_base = os.path.join(tmpdir, "base_m.xsf")
    run_cmd_stream(
        [atomsk, "--create", "fcc", str(a), element, "orient", o1, o2, o3, unit_m_base, "-overwrite"],
        cwd=tmpdir, timeout_s=timeout_s, log_path=None
    )

    unit_t_base = os.path.join(tmpdir, "base_t.xsf")
    run_cmd_stream(
        [atomsk, unit_m_base, "-mirror", "0", stack_axis, "-wrap", unit_t_base, "-overwrite"],
        cwd=tmpdir, timeout_s=timeout_s, log_path=None
    )

    # 2. Duplicate to thickness
    # Matrix
    m_layer = os.path.join(tmpdir, "layer_m.xsf")
    if stack_axis == "Y":
        d_args_m = [str(dup_x), str(matrix_h), str(dup_z)]
        d_args_t = [str(dup_x), str(twin_h), str(dup_z)]
    elif stack_axis == "Z":
        d_args_m = [str(dup_x), str(dup_z), str(matrix_h)]
        d_args_t = [str(dup_x), str(dup_z), str(twin_h)]
    else:
        d_args_m = [str(matrix_h), str(dup_x), str(dup_z)]
        d_args_t = [str(twin_h), str(dup_x), str(dup_z)]

    run_cmd_stream([atomsk, unit_m_base, "-duplicate"] + d_args_m + [m_layer, "-overwrite"], cwd=tmpdir,
                   timeout_s=timeout_s)

    # Twin
    t_layer = os.path.join(tmpdir, "layer_t.xsf")
    run_cmd_stream([atomsk, unit_t_base, "-duplicate"] + d_args_t + [t_layer, "-overwrite"], cwd=tmpdir,
                   timeout_s=timeout_s)

    # 3. Merge
    seed_out = os.path.join(tmpdir, f"seed_m{matrix_h}_t{twin_h}.xsf")
    run_cmd_stream([atomsk, "--merge", stack_axis, "2", m_layer, t_layer, seed_out, "-overwrite"], cwd=tmpdir,
                   timeout_s=timeout_s)

    return seed_out


def merge_lammps_files(file_list, intervals, box_dims, output_path):
    """
    Merge multiple LAMMPS data files by slicing atoms based on Y-intervals.
    file_list: list of paths to .lmp files
    intervals: list of (y_min, y_max) tuples corresponding to each file
    box_dims: [x, y, z] box dimensions
    """
    print(f"[INFO] Merging {len(file_list)} blocks into final structure...")

    # We read all atoms, filter them, and write a new file.
    # Since all files share the same box and atom types, we can just copy the header from the first one.

    all_atoms = []  # list of lines
    atom_id_counter = 1

    header_lines = []
    box_bounds = []

    for i, (fpath, (ymin, ymax)) in enumerate(zip(file_list, intervals)):
        print(f"  - Processing Block {i}: {os.path.basename(fpath)} (Keep {ymin:.1f} < y < {ymax:.1f})")

        with open(fpath, 'r') as f:
            lines = f.readlines()

        # Parse header if first file
        atom_start_idx = -1
        for idx, line in enumerate(lines):
            if "Atoms" in line:
                atom_start_idx = idx + 2  # Skip "Atoms" and blank line
                break
            if i == 0:
                header_lines.append(line)

        if atom_start_idx == -1:
            print(f"[WARN] No Atoms section found in {fpath}")
            continue

        # Filter atoms
        # LAMMPS atom line format: id type x y z ... (or id molecule type x y z ...)
        # We need to find which column is y. Usually for atom_style atomic: id type x y z
        # For full: id mol type q x y z
        # Let's assume standard atomic/charge/full where x,y,z are consecutive floats.
        # We'll try to detect by looking at the line.

        for line in lines[atom_start_idx:]:
            parts = line.split()
            if len(parts) < 5: continue

            # Heuristic: try to find 3 floats
            # Usually column 2,3,4 (0-indexed) or 3,4,5
            # Atomsk usually outputs 'id type x y z' or 'id type q x y z'
            # Let's parse all as floats and check bounds

            try:
                # Atomsk default for 'atomsk --polycrystal ... lmp' is usually 'atomic' or 'charge'
                # atomic: id type x y z
                # charge: id type q x y z
                # We can check the header for "atom types" etc but not atom style.
                # Let's assume the Y coordinate is the 4th column (index 3) if atomic, or 5th (index 4) if charge?
                # Safer: check x,y,z against box dims?
                # Actually, let's just assume Atomsk lmp output: id type x y z

                # Check column count
                # 5 columns -> id type x y z
                # 6 columns -> id type q x y z ?? No, charge is 6 usually?
                # Let's use the provided file structure if we can see it.
                # Assuming 'id type x y z' (5 cols) or 'id type x y z nx ny nz' (8 cols)

                # We will re-index IDs anyway.
                # Let's parse coords.
                # Try index 2,3,4 (atomic)
                y_val = float(parts[3])  # default guess

                # Check if this makes sense
                # If y_val is within [0, box_y], good.

                if ymin <= y_val <= ymax:
                    # Keep this atom
                    # Re-assign ID
                    parts[0] = str(atom_id_counter)
                    atom_id_counter += 1
                    all_atoms.append(" ".join(parts) + "\n")

            except ValueError:
                continue

    # Write output
    with open(output_path, 'w') as f:
        # Update atom count in header
        for line in header_lines:
            if "atoms" in line and len(line.split()) == 2:
                f.write(f"{len(all_atoms)} atoms\n")
            else:
                f.write(line)

        f.write("\nAtoms\n\n")
        for line in all_atoms:
            f.write(line)

    return len(all_atoms)


def generate_layer_sequence(cfg, mode_override=None):
    mode = mode_override if mode_override else cfg.TWIN_MODE
    layers = []

    # === SPECIAL HANDLING FOR EXACT LAYER COUNT (4-5 layers) ===
    # Since Atomsk requires duplicating Unit Cells (3 layers per cell),
    # we cannot directly get 4 or 5 atomic layers via -duplicate.
    # We must use -cut to trim the excess layers.
    # However, this function only returns the duplication count (integers).
    #
    # Workaround:
    # If the user wants 4 or 5 layers, they probably mean 4 or 5 ATOMIC layers.
    # 1 Unit Cell = 3 layers.
    # 2 Unit Cells = 6 layers.
    # To get 4 or 5, we have to generate 2 Unit Cells (6 layers) and rely on the user/post-process to cut it,
    # OR we can try to create a fractional unit cell seed? Atomsk supports -cut.
    #
    # But for now, since we are constrained to integer unit cells in this script:
    # 1 => 3 layers
    # 2 => 6 layers
    #
    # If the user REALLY wants 4-5 layers, they have to accept 3 (1 cell) or 6 (2 cells).
    # 1.5 cells is not supported by -duplicate.
    #
    # Alternative: Use 2 cells (6 layers) as the closest approximation that covers the range.

    if mode == 'random':
        count = getattr(cfg, "RANDOM_LAYER_COUNT", 10)
        min_l = cfg.RANDOM_MIN_LAYERS
        max_l = cfg.RANDOM_MAX_LAYERS
        for _ in range(count):
            layers.append(random.randint(min_l, max_l))
    elif mode == 'equal':
        thick = cfg.EQUAL_THICKNESS
        layers = [thick, thick, thick, thick]
    elif mode == 'all_Equal':
        count = getattr(cfg, "ALL_EQUAL_LAYER_COUNT", 10)
        m_h = getattr(cfg, "ALL_EQUAL_MATRIX_THICKNESS", 5)
        t_h = getattr(cfg, "ALL_EQUAL_TWIN_THICKNESS", 1)
        # Force strict alternation: M, T, M, T...
        # Ensure it starts with M and the total length is respected.
        for i in range(count):
            if i % 2 == 0:
                layers.append(m_h)
            else:
                layers.append(t_h)
    elif mode == 'ratio':
        total_layers = getattr(cfg, "RATIO_TOTAL_LAYERS", 20)
        frac = cfg.RATIO_TWIN_FRACTION
        min_h = cfg.RATIO_MIN_THICKNESS
        target_twin_h = int(total_layers * frac)
        target_matrix_h = total_layers - target_twin_h
        curr_m = 0
        curr_t = 0
        while (curr_m < target_matrix_h) or (curr_t < target_twin_h):
            rem_m = target_matrix_h - curr_m
            if rem_m > 0:
                h = random.randint(min_h, max(min_h, int(rem_m / 2) + min_h))
                h = min(h, rem_m)
                layers.append(h)
                curr_m += h
            elif curr_t < target_twin_h:
                pass
            rem_t = target_twin_h - curr_t
            if rem_t > 0:
                h = random.randint(min_h, max(min_h, int(rem_t / 2) + min_h))
                h = min(h, rem_t)
                layers.append(h)
                curr_t += h
        layers = [x for x in layers if x > 0]
    else:
        layers = [5, 5, 5, 5]
    return layers


def build_layered_seed(atomsk, cfg, tmpdir, timeout_s, out_dir, mode_override=None):
    # Create a unique sub-directory for this seed generation to avoid file conflicts
    # This prevents issues where 'unit_matrix_base.xsf' from a previous run (or zone)
    # interferes with the current one, or where existing output files trigger interactive prompts.
    seed_tmp_dir = os.path.join(tmpdir, f"seed_build_{time.time_ns()}")
    os.makedirs(seed_tmp_dir, exist_ok=True)

    element = cfg.ELEMENT
    a = cfg.LATTICE_A
    orient = cfg.TWIN_ORIENT
    stack_axis = getattr(cfg, "STACK_AXIS", "Y")
    dup_x = cfg.SEED_DUP_X
    dup_z = cfg.SEED_DUP_Z
    o1, o2, o3 = [_sanitize_miller(x) for x in orient]

    unit_m_base = os.path.join(seed_tmp_dir, "unit_matrix_base.xsf")
    run_cmd_stream(
        [atomsk, "--create", "fcc", str(a), element, "orient", o1, o2, o3, "-wrap", unit_m_base, "-overwrite"],
        cwd=seed_tmp_dir, timeout_s=timeout_s, log_path=None
    )

    unit_t_base = os.path.join(seed_tmp_dir, "unit_twin_base.xsf")
    run_cmd_stream(
        [atomsk, unit_m_base, "-mirror", "0", stack_axis, "-wrap", unit_t_base, "-overwrite"],
        cwd=seed_tmp_dir, timeout_s=timeout_s, log_path=None
    )

    layers = generate_layer_sequence(cfg, mode_override=mode_override)
    print(f"[INFO] Layer sequence (thickness in unit cells along {stack_axis}): {layers}")

    # Check for Twin Enable Switch
    # Default to True if not present (backward compatibility)
    # For GNT zones, it is in cfg.ENABLE_TWIN (via ZoneConfig)
    # For NT global, it is in cfg.NT_ENABLE_TWIN (via Config)
    enable_twin = getattr(cfg, "ENABLE_TWIN", getattr(cfg, "NT_ENABLE_TWIN", True))

    if not enable_twin:
        print("[INFO] Twin generation DISABLED. Generating single crystal seed with total thickness.")
        # Optimization: If twin is disabled, we can just use the wrapped unit matrix base as seed.
        # Atomsk --polycrystal will automatically tile it to fill the grains.
        # This avoids creating huge seed files via -duplicate and avoids potential overlap issues.
        # The only downside is if the user wanted a specific non-tiled seed size, but for polycrystal it doesn't matter.
        return unit_m_base

        # Original Logic (Disabled for Optimization)
        # total_h = sum(layers)
        # seed_out_tmp = os.path.join(seed_tmp_dir, "seed_single_crystal.xsf")
        #
        # if stack_axis == "Y":
        #     d_args = [str(dup_x), str(total_h), str(dup_z)]
        # elif stack_axis == "Z":
        #     d_args = [str(dup_x), str(dup_z), str(total_h)]
        # else:
        #     d_args = [str(total_h), str(dup_x), str(dup_z)]
        #
        # run_cmd_stream(
        #     [atomsk, unit_m_base, "-duplicate"] + d_args + [seed_out_tmp, "-overwrite"],
        #     cwd=seed_tmp_dir, timeout_s=timeout_s,
        #     log_path=None
        # )
        # return seed_out_tmp

    layer_files = []
    for i, h in enumerate(layers):
        is_twin = (i % 2 != 0)
        base_file = unit_t_base if is_twin else unit_m_base
        layer_fname = f"layer_{i:03d}.xsf"
        layer_path = os.path.join(seed_tmp_dir, layer_fname)

        if stack_axis == "Y":
            d_args = [str(dup_x), str(h), str(dup_z)]
        elif stack_axis == "Z":
            d_args = [str(dup_x), str(dup_z), str(h)]
        else:
            d_args = [str(h), str(dup_x), str(dup_z)]

        run_cmd_stream(
            [atomsk, base_file, "-duplicate"] + d_args + [layer_path, "-overwrite"],
            cwd=seed_tmp_dir, timeout_s=timeout_s,
            log_path=None
        )
        layer_files.append(layer_fname)

    seed_out_tmp = os.path.join(seed_tmp_dir, "layered_seed.xsf")
    merge_cmd = [atomsk, "--merge", stack_axis, str(len(layer_files))] + layer_files + [seed_out_tmp, "-overwrite"]
    run_cmd_stream(
        merge_cmd,
        cwd=seed_tmp_dir, timeout_s=timeout_s, log_path=None
    )

    # Copy the final seed back to the main tmpdir (or return the full path)
    # We return the full path in the sub-tmpdir
    return seed_out_tmp


def build_poly(atomsk, seed_xsf, poly_txt, tmpdir, timeout_s, out_dir):
    out = os.path.join(tmpdir, "poly.lmp")
    run_cmd_stream(
        [atomsk, "--polycrystal", seed_xsf, poly_txt, out, "-wrap", "-overwrite"],
        cwd=tmpdir,
        timeout_s=timeout_s,
        log_path=os.path.join(out_dir, "log_step4_polycrystal.txt"),
    )
    return out


def run_generation(cfg):
    # 生成带时间戳的输出目录名
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    base_out_dir = getattr(cfg, "OUT_DIR", "./output_project")
    mode = cfg.TWIN_MODE
    out_dir = f"{base_out_dir}_{mode}_{timestamp}"

    ensure_dir(out_dir)
    # Convert out_dir to absolute path to avoid issues when running subprocesses in tmpdir
    out_dir = os.path.abspath(out_dir)
    print(f"[INFO] Output directory: {out_dir}")

    # --- New Feature: Auto-Save Run Parameters ---
    params_log_path = os.path.join(out_dir, "run_parameters.txt")
    with open(params_log_path, 'w', encoding='utf-8') as f:
        f.write("=" * 40 + "\n")
        f.write(f"Run Timestamp: {timestamp}\n")
        f.write(f"Output Directory: {out_dir}\n")
        f.write("=" * 40 + "\n\n")

        # Dump global config attributes
        f.write("[Global Configuration]\n")
        for key in dir(cfg):
            if key.isupper():
                val = getattr(cfg, key)
                f.write(f"{key} = {val}\n")
        f.write("\n")

        # Dump GNT Zones detail if applicable
        if cfg.TWIN_MODE == 'GNT':
            f.write("[GNT Zones Detail]\n")
            zones = getattr(cfg, 'GNT_ZONES', [])
            for i, z in enumerate(zones):
                f.write(f"  Zone {i}:\n")
                for k, v in z.items():
                    f.write(f"    {k}: {v}\n")
                f.write("\n")

    print(f"[INFO] Parameters saved to: {params_log_path}")
    # ---------------------------------------------

    atomsk = find_atomsk()
    timeout_s = cfg.ATOMSK_TIMEOUT
    max_atoms = cfg.MAX_ATOMS

    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"[INFO] Working in temp dir: {tmpdir}")

        mode = cfg.TWIN_MODE

        if mode == 'GNT':
            # --- New Logic for GNT (Ghost Node / Multi-Pass Strategy) ---
            print("[INFO] Running GNT Mode (Multi-Pass Ghost Node Strategy)...")
            zones = getattr(cfg, 'GNT_ZONES', [])
            if not zones:
                raise ValueError("GNT_ZONES not defined in setting.py")

            total_ratio = sum(z['ratio'] for z in zones)
            total_h = cfg.POLY_BOX[1]
            box_dims = cfg.POLY_BOX

            # 1. Generate ALL nodes for ALL zones
            all_nodes = []  # List of Node objects

            current_y_base = 0.0

            print(f"[INFO] Generating global node distribution...")
            for i, z_conf in enumerate(zones):
                h_ratio = z_conf['ratio']
                zone_h = total_h * (h_ratio / total_ratio)
                y_min = current_y_base
                y_max = current_y_base + zone_h
                current_y_base += zone_h

                n_grains = z_conf['n_grains']

                print(f"  - Zone {i}: Y [{y_min:.1f}, {y_max:.1f}], Grains: {n_grains}")

                for _ in range(n_grains):
                    x = random.uniform(0, box_dims[0])
                    y = random.uniform(y_min, y_max)
                    z = random.uniform(0, box_dims[2])

                    # Random Orientation Logic (based on Zone setting)
                    is_random = z_conf.get('random_orient', True)
                    # Support independent 3D rotation switches
                    rx = z_conf.get('random_orient_X', True)
                    ry = z_conf.get('random_orient_Y', True)
                    rz = z_conf.get('random_orient_Z', True)

                    alpha = 0.0
                    beta = 0.0
                    gamma = 0.0

                    if rx and ry and rz:
                        # Full 3D random
                        alpha = random.uniform(0, 360)
                        beta = random.uniform(0, 360)
                        gamma = random.uniform(0, 360)
                    else:
                        # Fallback logic (similar to write_poly_txt)
                        if rx or rz:
                            alpha = random.uniform(0, 360) if rz else 0.0
                            beta = random.uniform(0, 360) if ry else 0.0
                            gamma = random.uniform(0, 360) if rx else 0.0
                        else:
                            # Only Y (or none)
                            alpha = 0.0
                            beta = random.uniform(0, 360) if ry else 0.0
                            gamma = 0.0

                    all_nodes.append(Node(x, y, z, alpha, beta, gamma, i))

            print(f"[INFO] Total Nodes: {len(all_nodes)}")

            # 2. Write poly_all.txt
            poly_all_path = os.path.join(tmpdir, "poly_all.txt")
            with open(poly_all_path, 'w') as f:
                f.write(f"box {box_dims[0]} {box_dims[1]} {box_dims[2]}\n")
                for n in all_nodes:
                    f.write(f"node {n.x:.4f} {n.y:.4f} {n.z:.4f} {n.a:.2f} {n.b:.2f} {n.g:.2f}\n")

            # 3. Multi-Pass Generation
            collected_atom_lines = []

            for i, z_conf in enumerate(zones):
                print(f"\n[STEP] Processing Zone {i} Pass...")

                # A. Build Seed for this Zone
                z_cfg = ZoneConfig(cfg, z_conf)
                # Ensure ZoneConfig has the new random orient attributes before accessing
                # (Though they are initialized in ZoneConfig.__init__, access here is just for printing)
                rx = getattr(z_cfg, 'NT_GRAIN_ORIENT_RANDOM_X', True)
                ry = getattr(z_cfg, 'NT_GRAIN_ORIENT_RANDOM_Y', True)
                rz = getattr(z_cfg, 'NT_GRAIN_ORIENT_RANDOM_Z', True)
                print(f"       Config -> Twin: {z_cfg.ENABLE_TWIN}, Random Orient: X:{rx} Y:{ry} Z:{rz}")

                seed_mode = z_conf.get('mode', 'equal')
                seed_i = build_layered_seed(atomsk, z_cfg, tmpdir, timeout_s, out_dir, mode_override=seed_mode)

                # B. Run Polycrystal (Full Box, All Nodes)
                # Output to CFG to preserve GrainID
                temp_cfg = os.path.join(tmpdir, f"temp_pass_{i}.cfg")

                # Use stream logging but maybe suppress the massive "Generating grain #..." output?
                # No, atomsk output is hard to filter without losing error messages.
                # We can just run it.
                # Note: We REMOVE "-wrap" here to avoid Atomsk NaN errors and handle wrapping in Python manually.
                print(f"       [RUN] Atomsk Polycrystal for Zone {i}...")
                run_cmd_stream(
                    [atomsk, "--polycrystal", seed_i, poly_all_path, temp_cfg, "-overwrite"],
                    cwd=tmpdir, timeout_s=timeout_s, log_path=None
                )

                # C. Filter Atoms
                # Only keep atoms belonging to grains defined in Zone i
                if os.path.exists(temp_cfg):
                    valid_atoms = parse_cfg_and_filter(temp_cfg, i, all_nodes, default_box=box_dims,
                                                       target_box=box_dims)
                    print(f"       Extracted {len(valid_atoms)} atoms for Zone {i}.")
                    collected_atom_lines.extend(valid_atoms)
                else:
                    print(f"[WARN] {temp_cfg} not found. Maybe Atomsk failed to write it?")

            # 4. Merge and Output
            print(f"\n[INFO] Merging {len(collected_atom_lines)} atoms...")

            # Convert to final LMP with cleanup
            try:
                raw_val = getattr(cfg, 'MERGE_DOUBLES_DIST', "NOT_SET")
                print(f"[DEBUG] cfg.MERGE_DOUBLES_DIST = {raw_val}")
            except:
                pass

            remove_dist = getattr(cfg, 'MERGE_DOUBLES_DIST', 1.2)
            print(f"[INFO] Removing overlapped atoms (dist < {remove_dist} A)...")

            poly_out = os.path.join(out_dir, "poly.lmp")

            # Write merged atoms directly to LAMMPS data file
            temp_merged_lmp = os.path.join(tmpdir, "temp_merged.lmp")
            write_merged_lammps(temp_merged_lmp, collected_atom_lines, box_dims, cfg.ELEMENT)

            # Then remove doubles (now using Average Merge)
            if HAS_NUMPY:
                remove_doubles_average(temp_merged_lmp, poly_out, remove_dist)
            else:
                # Fallback to Atomsk delete
                run_cmd_stream(
                    [atomsk, temp_merged_lmp, "-remove-doubles", str(remove_dist), poly_out, "-overwrite"],
                    cwd=tmpdir, timeout_s=timeout_s, log_path=None
                )

            n_final = parse_lammps_atoms_count(poly_out)
            print(f"\n[INFO] Final GNT structure atoms: {n_final}")

        elif mode == 'NT':
            # --- Logic for NT (Uniform Grain, Adjustable Twin) ---
            internal_mode = getattr(cfg, 'NT_INTERNAL_MODE', 'all_Equal')
            print(f"[INFO] Running NT mode with internal structure: {internal_mode}")

            # 1. Build Layered Seed using the chosen internal mode
            seed = build_layered_seed(atomsk, cfg, tmpdir, timeout_s, out_dir, mode_override=internal_mode)

            # Save seed backup
            seed_final = os.path.join(out_dir, "seed_layered.xsf")
            shutil.copy(seed, seed_final)

            # 2. Generate Poly.txt (Grid-based nodes)
            poly_txt = os.path.join(tmpdir, "poly.txt")
            write_poly_txt(poly_txt, cfg)
            shutil.copy(poly_txt, os.path.join(out_dir, "poly.txt"))

            # 3. Polycrystal Generation
            poly_lmp = build_poly(atomsk, seed, poly_txt, tmpdir, timeout_s, out_dir)

            # 4. Cleanup and Output
            remove_dist = getattr(cfg, 'MERGE_DOUBLES_DIST', 1.2)
            print(f"[INFO] Removing overlapped atoms (dist < {remove_dist} A)...")
            poly_out = os.path.join(out_dir, "poly.lmp")

            # Need temp LMP first because remove_doubles_average reads LMP
            temp_poly_lmp = os.path.join(tmpdir, "temp_poly.lmp")
            shutil.copy(poly_lmp, temp_poly_lmp)

            if HAS_NUMPY:
                remove_doubles_average(temp_poly_lmp, poly_out, remove_dist)
            else:
                run_cmd_stream(
                    [atomsk, poly_lmp, "-remove-doubles", str(remove_dist), poly_out, "-overwrite"],
                    cwd=tmpdir, timeout_s=timeout_s, log_path=None
                )

            n_final = parse_lammps_atoms_count(poly_out)
            print(f"\n[INFO] NT structure atoms: {n_final}")

        else:
            # --- Original Logic for Random/Equal/Ratio/all_Equal ---
            # 1. 生成多层孪晶种子
            seed = build_layered_seed(atomsk, cfg, tmpdir, timeout_s, out_dir)

            # 保存种子备份
            seed_final = os.path.join(out_dir, "seed_layered.xsf")
            shutil.copy(seed, seed_final)

            # SANITY CHECK: Try converting seed to lmp to check validity
            print("[DEBUG] Verifying seed file validity...")
            seed_check = os.path.join(tmpdir, "seed_check.lmp")
            run_cmd_stream([atomsk, seed, seed_check, "-overwrite"], cwd=tmpdir, timeout_s=timeout_s, log_path=None)
            print("[DEBUG] Seed file is valid.")

            # 2. 多晶生成
            poly_txt = os.path.join(tmpdir, "poly.txt")
            write_poly_txt(poly_txt, cfg)
            shutil.copy(poly_txt, os.path.join(out_dir, "poly.txt"))

            poly_lmp = build_poly(atomsk, seed, poly_txt, tmpdir, timeout_s, out_dir)

            # 3. Cleanup and Output
            remove_dist = getattr(cfg, 'MERGE_DOUBLES_DIST', 1.2)
            print(f"[INFO] Removing overlapped atoms (dist < {remove_dist} A)...")
            poly_out = os.path.join(out_dir, "poly.lmp")

            temp_poly_lmp = os.path.join(tmpdir, "temp_poly.lmp")
            shutil.copy(poly_lmp, temp_poly_lmp)

            if HAS_NUMPY:
                remove_doubles_average(temp_poly_lmp, poly_out, remove_dist)
            else:
                run_cmd_stream(
                    [atomsk, poly_lmp, "-remove-doubles", str(remove_dist), poly_out, "-overwrite"],
                    cwd=tmpdir, timeout_s=timeout_s, log_path=None
                )

            n_final = parse_lammps_atoms_count(poly_out)
            print(f"\n[INFO] poly.lmp real atoms = {n_final}")

        # Final check
        if n_final > max_atoms:
            print(f"[WARN] Atom count {n_final} exceeds limit {max_atoms}!")

    print("\n[DONE] Generation Complete!")
    print(f"  Poly: {poly_out}")
    print(f"  Logs: {out_dir}")
