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

def write_poly_txt(path, cfg):
    box = cfg.POLY_BOX
    ngrains = cfg.POLY_GRAINS
    mode = cfg.TWIN_MODE

    with open(path, "w", encoding="utf-8") as f:
        f.write(f"box {box[0]} {box[1]} {box[2]}\n")
        
        if mode == 'NT':
            # "Mode 6 based NT": Random position + Parallel Twin Planes
            # 1. Random Position (like all_Equal)
            print(f"[INFO] Generating {ngrains} grains with Random Position & Parallel Twin Planes (NT mode).")
            
            for _ in range(ngrains):
                x = random.uniform(0, box[0])
                y = random.uniform(0, box[1])
                z = random.uniform(0, box[2])
                
                # 2. Parallel Twin Planes
                # Twin Plane Normal is Y axis [111].
                # To keep planes parallel, we must NOT tilt Y axis.
                # We can only rotate AROUND Y axis.
                # Atomsk Euler (Z-Y-Z): 0, theta, 0 => Rotate around Y by theta.
                
                alpha = 0.0
                beta = random.uniform(0, 360)
                gamma = 0.0
                
                f.write(f"node {x:.4f} {y:.4f} {z:.4f} {alpha:.2f} {beta:.2f} {gamma:.2f}\n")
                        
        elif mode == 'GNT':
            # Explicit nodes with fixed orientation
            # Angle 0 0 0 ensures the seed orientation is preserved (Y=[111] aligned with Box Y)
            print(f"[INFO] Generating {ngrains} grains with fixed orientation (0,0,0) for {mode} mode.")
            
            for _ in range(ngrains):
                x = random.uniform(0, box[0])
                z = random.uniform(0, box[2])
                
                if getattr(cfg, 'GNT_GRAIN_GRADIENT', False):
                    # Gradient density: High density at Top (Y=H), Low density at Bottom (Y=0)
                    u = random.random()
                    y = box[1] * math.sqrt(u)
                else:
                    y = random.uniform(0, box[1])
                    
                # Format: node x y z alpha beta gamma
                f.write(f"node {x:.4f} {y:.4f} {z:.4f} 0.0 0.0 0.0\n")
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
    element = cfg.ELEMENT
    a = cfg.LATTICE_A
    orient = cfg.TWIN_ORIENT
    stack_axis = getattr(cfg, "STACK_AXIS", "Y")
    dup_x = cfg.SEED_DUP_X
    dup_z = cfg.SEED_DUP_Z
    o1, o2, o3 = [_sanitize_miller(x) for x in orient]
    unit_m_base = os.path.join(tmpdir, "unit_matrix_base.xsf")
    run_cmd_stream(
        [atomsk, "--create", "fcc", str(a), element, "orient", o1, o2, o3, unit_m_base, "-overwrite"],
        cwd=tmpdir, timeout_s=timeout_s, log_path=os.path.join(out_dir, "log_step1_create.txt")
    )
    unit_t_base = os.path.join(tmpdir, "unit_twin_base.xsf")
    run_cmd_stream(
        [atomsk, unit_m_base, "-mirror", "0", stack_axis, "-wrap", unit_t_base, "-overwrite"],
        cwd=tmpdir, timeout_s=timeout_s, log_path=os.path.join(out_dir, "log_step2_mirror.txt")
    )
    layers = generate_layer_sequence(cfg, mode_override=mode_override)
    print(f"[INFO] Layer sequence (thickness in unit cells along {stack_axis}): {layers}")
    layer_files = []
    for i, h in enumerate(layers):
        is_twin = (i % 2 != 0)
        base_file = unit_t_base if is_twin else unit_m_base
        layer_fname = f"layer_{i:03d}.xsf"
        layer_path = os.path.join(tmpdir, layer_fname)
        if stack_axis == "Y":
            d_args = [str(dup_x), str(h), str(dup_z)]
        elif stack_axis == "Z":
            d_args = [str(dup_x), str(dup_z), str(h)]
        else:
            d_args = [str(h), str(dup_x), str(dup_z)]
        run_cmd_stream(
            [atomsk, base_file, "-duplicate"] + d_args + [layer_path, "-overwrite"],
            cwd=tmpdir, timeout_s=timeout_s, 
            log_path=None
        )
        layer_files.append(layer_fname)
    seed_out = os.path.join(tmpdir, "layered_seed.xsf")
    merge_cmd = [atomsk, "--merge", stack_axis, str(len(layer_files))] + layer_files + [seed_out, "-overwrite"]
    run_cmd_stream(
        merge_cmd,
        cwd=tmpdir, timeout_s=timeout_s, log_path=os.path.join(out_dir, "log_step3_merge.txt")
    )
    return seed_out

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
    out_dir = f"{base_out_dir}_{timestamp}"
    
    ensure_dir(out_dir)
    print(f"[INFO] Output directory: {out_dir}")
    
    atomsk = find_atomsk()
    timeout_s = cfg.ATOMSK_TIMEOUT
    max_atoms = cfg.MAX_ATOMS
    
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"[INFO] Working in temp dir: {tmpdir}")
        
        mode = cfg.TWIN_MODE
        
        if mode == 'GNT':
            # --- New Logic for GNT ---
            spacings = getattr(cfg, 'GNT_SPACINGS', [10, 6, 2])
            twin_h = getattr(cfg, 'GNT_TWIN_THICKNESS', 1)

            # 1. Generate Poly.txt (Nodes)
            poly_txt = os.path.join(tmpdir, "poly.txt")
            write_poly_txt(poly_txt, cfg)
            shutil.copy(poly_txt, os.path.join(out_dir, "poly.txt"))
            
            # 2. Generate Blocks
            block_files = []
            box_h = cfg.POLY_BOX[1]
            n_blocks = len(spacings)
            block_h = box_h / n_blocks
            intervals = []
            
            for i, sp in enumerate(spacings):
                print(f"\n[STEP] Generating Block {i+1}/{n_blocks} with spacing {sp}...")
                
                # a. Build Seed
                seed_i = build_simple_seed(atomsk, cfg, sp, twin_h, tmpdir, timeout_s)
                
                # b. Build Polycrystal (Full Box)
                block_name = f"block_{i}.lmp"
                block_path = os.path.join(tmpdir, block_name)
                
                # We use the SAME poly.txt so grain boundaries match
                run_cmd_stream(
                    [atomsk, "--polycrystal", seed_i, poly_txt, block_path, "-wrap", "-overwrite"],
                    cwd=tmpdir, timeout_s=timeout_s, log_path=None
                )
                
                block_files.append(block_path)
                
                # Define Interval
                ymin = i * block_h
                ymax = (i + 1) * block_h
                intervals.append((ymin, ymax))
                
            # 3. Merge and Slice
            poly_out = os.path.join(out_dir, "poly.lmp")
            n_final = merge_lammps_files(block_files, intervals, cfg.POLY_BOX, poly_out)
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
