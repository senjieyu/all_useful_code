# -*- coding: utf-8 -*-
# core_generator.py

import os
import shutil
import subprocess
import tempfile
import random

import math
import time

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

def parse_cfg_and_filter(cfg_path, valid_zone_idx, all_nodes):
    """
    Parse Atomsk CFG file and keep atoms belonging to grains in valid_zone_idx.
    Returns a list of atom lines (in LAMMPS data format style or similar).
    """
    valid_atoms = []
    
    with open(cfg_path, 'r') as f:
        lines = f.readlines()
    
    # Parse Header
    entry_count = 0
    aux_map = {} # aux_name -> index offset from coords
    
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
        elif ".NO_VELOCITY." in line:
            pass
        elif len(line.split()) == 1 and line.replace('.','').isdigit():
             # Mass line usually? No, Atomsk CFG:
             # Mass
             # Element
             # Coordinates...
             pass
        
        # Check for start of data
        # In Atomeye CFG, after header lines, we have Mass, Element, then data.
        # But we don't know exactly where header ends.
        # Usually checking for first line with 'entry_count' floats.
        try:
            parts = line.split()
            if len(parts) == entry_count:
                # This is likely a data line.
                # But wait, Mass is a single float. Element is a string.
                # So data lines are those with `entry_count` items.
                pass
        except:
            pass

    # Simplified parsing for Atomsk CFG
    # We look for lines with `entry_count` columns.
    # And we need to find which column is grainID.
    # Standard: s_x s_y s_z [aux...]
    # So grainID index = 3 + aux_map['grainID']
    
    if 'grainID' not in aux_map:
        print(f"[WARN] No grainID found in {cfg_path}, skipping filtering (keeping all? No, unsafe).")
        return []
        
    grain_col_idx = 3 + aux_map['grainID']
    
    for line in lines:
        parts = line.split()
        if len(parts) != entry_count:
            continue
            
        # Check if valid data line (floats)
        try:
            # Check grainID
            gid = int(float(parts[grain_col_idx]))
            
            # Atomsk GrainID is 1-based index corresponding to node order.
            # Map to 0-based index
            node_idx = gid - 1
            
            if 0 <= node_idx < len(all_nodes):
                if all_nodes[node_idx].zone_idx == valid_zone_idx:
                    # Keep this atom
                    # We store it as: (type, x, y, z) - wait, we need type?
                    # CFG doesn't store type ID, it implies type by blocks.
                    # But Atomsk usually writes 1 block for mono-element.
                    # Let's assume single element for now (Cu/Al).
                    # We need to convert reduced coords to absolute?
                    # Atomsk CFG output is usually in Angstroms if H is identity?
                    # Wait, H0(1,1)=100.
                    # If data is 0.5, it is reduced.
                    # We need to read H matrix to convert.
                    # For simplicity, let's assume we just want to merge them later.
                    # We can store the whole line and write back to a CFG?
                    valid_atoms.append(line)
        except ValueError:
            continue
            
    return valid_atoms

def write_merged_cfg(output_path, atom_lines, box_dims, element):
    """
    Write a valid CFG file with merged atoms.
    """
    with open(output_path, 'w') as f:
        f.write(f"Number of particles = {len(atom_lines)}\n")
        f.write("A = 1.0 Angstrom (basic length-scale)\n")
        f.write(f"H0(1,1) = {box_dims[0]:.6f}\n")
        f.write(f"H0(1,2) = 0.0\n")
        f.write(f"H0(1,3) = 0.0\n")
        f.write(f"H0(2,1) = 0.0\n")
        f.write(f"H0(2,2) = {box_dims[1]:.6f}\n")
        f.write(f"H0(2,3) = 0.0\n")
        f.write(f"H0(3,1) = 0.0\n")
        f.write(f"H0(3,2) = 0.0\n")
        f.write(f"H0(3,3) = {box_dims[2]:.6f}\n")
        f.write(".NO_VELOCITY.\n")
        f.write(f"entry_count = 4\n") # x, y, z, grainID
        f.write("auxiliary[0] = grainID\n")
        
        # Mass and Element (Approximation)
        # We need atomic mass.
        # cfg.ELEMENT -> mass?
        # Let's just put a dummy mass or try to look it up. 
        # Atomsk might complain if mass is wrong? No, it uses Element name.
        # But CFG format requires Mass line.
        # Let's use 1.0 for dummy if unknown, or hardcode common ones.
        mass_map = {"Cu": 63.546, "Al": 26.98, "Fe": 55.845, "Ni": 58.693}
        mass = mass_map.get(element, 1.0)
        
        f.write(f"{mass:.4f}\n")
        f.write(f"{element}\n")
        
        for line in atom_lines:
            f.write(line)
            
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
        
    run_cmd_stream([atomsk, unit_m_base, "-duplicate"] + d_args_m + [m_layer, "-overwrite"], cwd=tmpdir, timeout_s=timeout_s)
    
    # Twin
    t_layer = os.path.join(tmpdir, "layer_t.xsf")
    run_cmd_stream([atomsk, unit_t_base, "-duplicate"] + d_args_t + [t_layer, "-overwrite"], cwd=tmpdir, timeout_s=timeout_s)
    
    # 3. Merge
    seed_out = os.path.join(tmpdir, f"seed_m{matrix_h}_t{twin_h}.xsf")
    run_cmd_stream([atomsk, "--merge", stack_axis, "2", m_layer, t_layer, seed_out, "-overwrite"], cwd=tmpdir, timeout_s=timeout_s)
    
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
    
    all_atoms = [] # list of lines
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
                atom_start_idx = idx + 2 # Skip "Atoms" and blank line
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
                y_val = float(parts[3]) # default guess
                
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
                h = random.randint(min_h, max(min_h, int(rem_m/2) + min_h))
                h = min(h, rem_m)
                layers.append(h)
                curr_m += h
            elif curr_t < target_twin_h:
                pass 
            rem_t = target_twin_h - curr_t
            if rem_t > 0:
                h = random.randint(min_h, max(min_h, int(rem_t/2) + min_h))
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
        [atomsk, "--create", "fcc", str(a), element, "orient", o1, o2, o3, unit_m_base, "-overwrite"],
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
        total_h = sum(layers)
        seed_out_tmp = os.path.join(seed_tmp_dir, "seed_single_crystal.xsf")
        
        if stack_axis == "Y":
            d_args = [str(dup_x), str(total_h), str(dup_z)]
        elif stack_axis == "Z":
            d_args = [str(dup_x), str(dup_z), str(total_h)]
        else:
            d_args = [str(total_h), str(dup_x), str(dup_z)]
            
        run_cmd_stream(
            [atomsk, unit_m_base, "-duplicate"] + d_args + [seed_out_tmp, "-overwrite"],
            cwd=seed_tmp_dir, timeout_s=timeout_s, 
            log_path=None
        )
        return seed_out_tmp
    
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
            all_nodes = [] # List of Node objects
            
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
                print(f"       [RUN] Atomsk Polycrystal for Zone {i}...")
                run_cmd_stream(
                    [atomsk, "--polycrystal", seed_i, poly_all_path, temp_cfg, "-wrap", "-overwrite"],
                    cwd=tmpdir, timeout_s=timeout_s, log_path=None
                )
                
                # C. Filter Atoms
                # Only keep atoms belonging to grains defined in Zone i
                valid_atoms = parse_cfg_and_filter(temp_cfg, i, all_nodes)
                print(f"       Extracted {len(valid_atoms)} atoms for Zone {i}.")
                collected_atom_lines.extend(valid_atoms)
                
            # 4. Merge and Output
            print(f"\n[INFO] Merging {len(collected_atom_lines)} atoms...")
            merged_cfg = os.path.join(tmpdir, "merged.cfg")
            write_merged_cfg(merged_cfg, collected_atom_lines, box_dims, cfg.ELEMENT)
            
            # Convert to final LMP
            poly_out = os.path.join(out_dir, "poly.lmp")
            run_cmd_stream(
                [atomsk, merged_cfg, poly_out, "-overwrite"],
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
            
            poly_out = os.path.join(out_dir, "poly.lmp")
            shutil.copy(poly_lmp, poly_out)
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
            
            poly_out = os.path.join(out_dir, "poly.lmp")
            shutil.copy(poly_lmp, poly_out)
            n_final = parse_lammps_atoms_count(poly_out)
            print(f"\n[INFO] poly.lmp real atoms = {n_final}")

        # Final check
        if n_final > max_atoms:
             print(f"[WARN] Atom count {n_final} exceeds limit {max_atoms}!")

    print("\n[DONE] Generation Complete!")
    print(f"  Poly: {poly_out}")
    print(f"  Logs: {out_dir}")
