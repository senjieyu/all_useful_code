# -*- coding: utf-8 -*-
# core_generator.py

import os
import shutil
import subprocess
import tempfile
import random

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

def write_poly_txt(path, box, ngrains):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"box {box[0]} {box[1]} {box[2]}\n")
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
def generate_layer_sequence(cfg):
    """
    根据配置生成每一层的厚度序列 (List[int])。
    假设偶数索引为 Matrix，奇数索引为 Twin。
    """
    mode = cfg.TWIN_MODE
    layers = []
    
    if mode == 'random':
        # 这里改为 "LAYER_COUNT" 而不是 "PAIR_COUNT" 以更灵活
        count = getattr(cfg, "RANDOM_LAYER_COUNT", 10) 
        min_l = cfg.RANDOM_MIN_LAYERS
        max_l = cfg.RANDOM_MAX_LAYERS
        for _ in range(count):
            layers.append(random.randint(min_l, max_l))
            
    elif mode == 'equal':
        # 等间距，生成一些层，交替即可
        # 数量任意，只要能构成一个周期即可，但为了 Atomsk 处理方便，生成 2~4 层
        thick = cfg.EQUAL_THICKNESS
        layers = [thick, thick, thick, thick]
            
    elif mode == 'ratio':
        total_layers = getattr(cfg, "RATIO_TOTAL_LAYERS", 20)
        frac = cfg.RATIO_TWIN_FRACTION
        min_h = cfg.RATIO_MIN_THICKNESS
        
        target_twin_h = int(total_layers * frac)
        target_matrix_h = total_layers - target_twin_h
        
        curr_m = 0
        curr_t = 0
        
        # 这里的 total_layers 是指厚度总和 (unit cells) 还是层数？
        # 假设是总厚度
        while (curr_m < target_matrix_h) or (curr_t < target_twin_h):
            # Matrix Layer
            rem_m = target_matrix_h - curr_m
            if rem_m > 0:
                h = random.randint(min_h, max(min_h, int(rem_m/2) + min_h))
                h = min(h, rem_m)
                layers.append(h)
                curr_m += h
            elif curr_t < target_twin_h:
                # 只有 twin 还没满，matrix 满了 -> 补个 0 或者跳过
                pass 
            
            # Twin Layer
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

def build_layered_seed(atomsk, cfg, tmpdir, timeout_s, out_dir):
    """
    生成多层孪晶种子 (按照用户要求的 Al 孪晶多晶建模步骤修改)
    Stacking Axis: Y
    """
    element = cfg.ELEMENT
    a = cfg.LATTICE_A
    orient = cfg.TWIN_ORIENT
    stack_axis = getattr(cfg, "STACK_AXIS", "Y") # Y axis is [111]
    
    dup_x = cfg.SEED_DUP_X
    dup_z = cfg.SEED_DUP_Z
    
    o1, o2, o3 = [_sanitize_miller(x) for x in orient]

    # 1. 基础单元 (Unit Cell Matrix) - 1x1x1
    # atomsk --create fcc 4.05 Al orient [11-2] [111] [-110] base.xsf
    unit_m_base = os.path.join(tmpdir, "unit_matrix_base.xsf")
    run_cmd_stream(
        [atomsk, "--create", "fcc", str(a), element, "orient", o1, o2, o3, unit_m_base],
        cwd=tmpdir, timeout_s=timeout_s, log_path=os.path.join(out_dir, "log_step1_create.txt")
    )
    
    # 2. 生成 Mirror Unit (Unit Twin)
    # atomsk base.xsf -mirror 0 Y -wrap mirror.xsf
    unit_t_base = os.path.join(tmpdir, "unit_twin_base.xsf")
    run_cmd_stream(
        [atomsk, unit_m_base, "-mirror", "0", stack_axis, "-wrap", unit_t_base],
        cwd=tmpdir, timeout_s=timeout_s, log_path=os.path.join(out_dir, "log_step2_mirror.txt")
    )

    # 3. 生成层序列
    layers = generate_layer_sequence(cfg)
    print(f"[INFO] Layer sequence (thickness in unit cells along {stack_axis}): {layers}")
    
    layer_files = []
    
    # 4. 生成每一层文件 (Duplicate)
    # duplicate X Y Z -> 如果 stack_axis 是 Y，则 thickness 放在 Y 位置
    # duplicate 1 3 1
    for i, h in enumerate(layers):
        is_twin = (i % 2 != 0)
        base_file = unit_t_base if is_twin else unit_m_base
        layer_fname = f"layer_{i:03d}.xsf"
        layer_path = os.path.join(tmpdir, layer_fname)
        
        # 构造 duplicate 参数
        # 假设 stack_axis == "Y"
        if stack_axis == "Y":
            d_args = [str(dup_x), str(h), str(dup_z)]
        elif stack_axis == "Z":
            d_args = [str(dup_x), str(dup_z), str(h)] # 这里的变量名有点乱，总之非 stack 方向用 dup_x/z
        else: # X
            d_args = [str(h), str(dup_x), str(dup_z)]
            
        run_cmd_stream(
            [atomsk, base_file, "-duplicate"] + d_args + [layer_path],
            cwd=tmpdir, timeout_s=timeout_s, 
            log_path=None
        )
        layer_files.append(layer_fname)

    # 5. 合并 (Merge)
    # atomsk --merge Y N file1 file2 ... out
    seed_out = os.path.join(tmpdir, "layered_seed.xsf")
    merge_cmd = [atomsk, "--merge", stack_axis, str(len(layer_files))] + layer_files + [seed_out]
    
    run_cmd_stream(
        merge_cmd,
        cwd=tmpdir, timeout_s=timeout_s, log_path=os.path.join(out_dir, "log_step3_merge.txt")
    )
    
    return seed_out

def build_poly(atomsk, seed_xsf, poly_txt, tmpdir, timeout_s, out_dir):
    out = os.path.join(tmpdir, "poly.lmp")
    run_cmd_stream(
        [atomsk, "--polycrystal", seed_xsf, poly_txt, out, "-wrap"],
        cwd=tmpdir,
        timeout_s=timeout_s,
        log_path=os.path.join(out_dir, "log_step4_polycrystal.txt"),
    )
    return out

import time

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
        
        # 1. 生成多层孪晶种子
        seed = build_layered_seed(atomsk, cfg, tmpdir, timeout_s, out_dir)
        
        # 保存种子备份
        seed_final = os.path.join(out_dir, "seed_layered.xsf")
        shutil.copy(seed, seed_final)
        
        # 2. 多晶生成
        poly_txt = os.path.join(tmpdir, "poly.txt")
        write_poly_txt(poly_txt, cfg.POLY_BOX, cfg.POLY_GRAINS)
        shutil.copy(poly_txt, os.path.join(out_dir, "poly.txt"))
        
        poly_lmp = build_poly(atomsk, seed, poly_txt, tmpdir, timeout_s, out_dir)
        
        # 3. 检查原子数
        if not os.path.exists(poly_lmp):
            raise RuntimeError(f"生成失败：未找到输出文件 {poly_lmp}。请检查 Atomsk 日志（如内存不足、参数错误等）。")
            
        n_atoms = parse_lammps_atoms_count(poly_lmp)
        print(f"\n[INFO] poly.lmp real atoms = {n_atoms}")
        if n_atoms > max_atoms:
            raise RuntimeError(
                f"真实原子数 {n_atoms} 超过 MAX_ATOMS={max_atoms}。\n"
                f"建议：减小 POLY_BOX / POLY_GRAINS / SEED_DUP_X/Y 或层厚度。"
            )
            
        poly_out = os.path.join(out_dir, "poly.lmp")
        shutil.copy(poly_lmp, poly_out)
        
    print("\n[DONE] Generation Complete!")
    print(f"  Seed: {seed_final}")
    print(f"  Poly: {poly_out}")
    print(f"  Logs: {out_dir}")
