# -*- coding: utf-8 -*-
import os
import model_construct as mc
import cons_batch_generate as batch

# ===== 新增：多晶 & block 模块 =====
import make_poly_and_blocks as mpb


def format_vec(v):
    return f"[{v[0]:6.3f}, {v[1]:6.3f}, {v[2]:6.3f}]"


def print_job_info(struct_type, mat_conf, tasks, dirs, v1, v2, info):
    """
    打印详细的任务仪表盘 (带晶向指数)
    """
    print("\n" + "=" * 80)
    print(f"{'STRUCTURAL GENERATION JOB DASHBOARD':^80}")
    print("=" * 80)

    # 1. 基础信息
    print(f"\n[Basic Information]")
    print(f"  - Element:       {mat_conf['element']}")
    print(f"  - Phase Struct:  {struct_type}")
    print(f"  - Lattice:       a={mat_conf['lattice_a']}")
    print(f"  - Supercell:     {mat_conf['size']}")
    print(f"  - Vacuum:        {mat_conf['vacuum']} Angstrom")

    # 2. 几何与滑移系信息
    print(f"\n[Slip Systems & Crystallography]")
    print(f"  ------------------------------------------------------------")
    print(f"  * Main Slip Path (v1):")
    print(f"      Direction : {info.get('main_label', 'N/A')}")
    print(f"      Type      : {info.get('main_desc', 'N/A')}")
    print(f"      Vector    : {format_vec(v1)}")
    print(f"  ------------------------------------------------------------")
    print(f"  * Secondary Slip Path (v2):")
    print(f"      Direction : {info.get('sec_label', 'N/A')}")
    print(f"      Type      : {info.get('sec_desc', 'N/A')}")
    print(f"      Vector    : {format_vec(v2)}")
    print(f"  ------------------------------------------------------------")

    # 3. 输出路径
    print(f"\n[Output Directories]")
    print(f"  - Batch Files:   {os.path.abspath(dirs['out'])}")

    # 4. 任务详情
    print(f"\n[Task Configuration]")
    total_est = 0

    def est_files(t_name, t_cfg):
        if not t_cfg['enable']:
            return 0, "DISABLED"
        base = 0
        if t_name == "1d_main":
            base = t_cfg['n_points']
        elif t_name == "1d_sec":
            base = t_cfg['n_points']
        elif t_name == "2d_surf":
            base = t_cfg['grid_n'] ** 2
        elif t_name == "key_pts":
            base = len(t_cfg['fractions'])
        elif t_name == "strain":
            base = len(t_cfg['strains']) * len(t_cfg['key_fracs'])
        elif t_name == "twin":
            base = 1
        total = base * (1 + t_cfg.get('n_perturb', 0))
        return total, f"Active  | Est: {total:<3} files"

    for key, cfg in tasks.items():
        c, s = est_files(key, cfg)
        total_est += c
        print(f"  - {key:<10}: {s}")

    print(f"\n  >> Total Estimated Files: {total_est}")
    print("=" * 80 + "\n")


# ======================================================================
# ============================== MAIN ===================================
# ======================================================================
if __name__ == "__main__":

    # --- 1. 用户配置 ---
    STRUCTURE_TYPE = "FCC"  # FCC / BCC / HCP

    MAT_CONFIG = {
        "element": "Cu",
        "lattice_a": 3.615,
        "lattice_c": None,
        "size": (3, 3, 12),
        "amp": 0.02,
        "vacuum": 0.0,
    }

    TASKS = {
        "1d_main": {"enable": True, "n_points": 21, "n_perturb": 0},
        "1d_sec": {"enable": True, "n_points": 21, "n_perturb": 0},
        "2d_surf": {"enable": True, "grid_n": 6, "n_perturb": 0},
        "key_pts": {"enable": True, "fractions": [0.0, 1.0], "n_perturb": 0},
        "strain": {"enable": True, "strains": [-0.01], "key_fracs": [0.0], "n_perturb": 0},
        "twin": {"enable": True, "n_perturb": 0}
    }

    # --- 2. 路径初始化 ---
    ROOT_DIR = os.getcwd()
    BASIC_MODEL_DIR = os.path.join(ROOT_DIR, "Basic_model1")
    BATCH_OUT_DIR = os.path.join(ROOT_DIR, "test_batch_file1")
    BATCH_RELAX_DIR = os.path.join(BATCH_OUT_DIR, "to_relax1")

    for d in [BASIC_MODEL_DIR, BATCH_OUT_DIR, BATCH_RELAX_DIR]:
        os.makedirs(d, exist_ok=True)

    dirs = {
        'base_dir': BASIC_MODEL_DIR,
        'out': BATCH_OUT_DIR,
        'relax': BATCH_RELAX_DIR
    }

    for key in TASKS:
        TASKS[key]['amp'] = MAT_CONFIG['amp']
        TASKS[key]['element'] = MAT_CONFIG['element']
        TASKS[key]['struct_type'] = STRUCTURE_TYPE

    # --- 3. 构建模型 ---
    print("Initializing Model Construction...")
    builder = mc.StructureFactory()

    if STRUCTURE_TYPE == "FCC":
        base_atoms, v_main, v_sec, info = builder.build_fcc(
            MAT_CONFIG['element'],
            MAT_CONFIG['lattice_a'],
            MAT_CONFIG['size'],
            vacuum=MAT_CONFIG['vacuum']
        )
    elif STRUCTURE_TYPE == "BCC":
        base_atoms, v_main, v_sec, info = builder.build_bcc(
            MAT_CONFIG['element'],
            MAT_CONFIG['lattice_a'],
            MAT_CONFIG['size'],
            vacuum=MAT_CONFIG['vacuum']
        )
    elif STRUCTURE_TYPE == "HCP":
        c_val = MAT_CONFIG['lattice_c'] if MAT_CONFIG['lattice_c'] else 1.633 * MAT_CONFIG['lattice_a']
        base_atoms, v_main, v_sec, info = builder.build_hcp(
            MAT_CONFIG['element'],
            MAT_CONFIG['lattice_a'],
            c_val,
            MAT_CONFIG['size'],
            vacuum=MAT_CONFIG['vacuum']
        )
    else:
        raise ValueError(f"Unknown Structure Type: {STRUCTURE_TYPE}")

    mc.AtomUtils.write_vasp_file(
        base_atoms,
        f"Base_{STRUCTURE_TYPE}_{MAT_CONFIG['element']}",
        BASIC_MODEL_DIR
    )

    # --- 4. 打印面板 ---
    print_job_info(STRUCTURE_TYPE, MAT_CONFIG, TASKS, dirs, v_main, v_sec, info)

    # --- 5. 执行原有任务 ---
    input("Press Enter to start generation...")

    total_gen = 0
    if TASKS['1d_main']['enable']:
        total_gen += batch.generate_1d_main(base_atoms, v_main, TASKS['1d_main'], dirs)
    if TASKS['1d_sec']['enable']:
        total_gen += batch.generate_1d_sec(base_atoms, v_sec, TASKS['1d_sec'], dirs)
    if TASKS['2d_surf']['enable']:
        raw_v1 = base_atoms.get_cell()[0] / MAT_CONFIG['size'][0]
        raw_v2 = base_atoms.get_cell()[1] / MAT_CONFIG['size'][1]
        total_gen += batch.generate_2d_surface(
            base_atoms, raw_v1, raw_v2, TASKS['2d_surf'], dirs
        )
    if TASKS['key_pts']['enable']:
        total_gen += batch.generate_key_points(base_atoms, v_main, TASKS['key_pts'], dirs)
    if TASKS['strain']['enable']:
        total_gen += batch.generate_strain(base_atoms, v_main, TASKS['strain'], dirs)
    if TASKS['twin']['enable']:
        total_gen += batch.generate_twin(base_atoms, TASKS['twin'], dirs)

    print(f"\nDone. Total files: {total_gen}")

    # ==================================================================
    # ============ 新增模块：多晶 + block（不影响原逻辑） ============
    # ==================================================================
    POLY_BLOCK_CONFIG = {
        "enable": True,
        "element": MAT_CONFIG["element"],
        "structure": STRUCTURE_TYPE,
        "lattice_a": MAT_CONFIG["lattice_a"],

        # --- 孪晶晶胞 ---
        "use_twin_unit": True,
        "twin_orient": ("[11-2]", "[111]", "[-110]"),
        "twin_duplicate": (1, 3, 1),

        # --- 多晶 ---
        "poly_box": (10, 10, 10),   #   150 150 150 Å
        "poly_grains": 2,# 5
        "max_poly_atoms": 300000,

        # --- 切块 ---
        "block_atom_range": (100, 300),
        "block_window_min": (12, 12, 12),
        "block_window_max": (20, 20, 20),
        "n_blocks": 40 ,

        # --- 随机性 ---
        "seed": 42,
        "block_perturb_amp": 0.0,
    }

    if POLY_BLOCK_CONFIG["enable"]:
        print("\n" + "=" * 80)
        print("[Polycrystal] Generating Voronoi polycrystal + blocks using Atomsk")
        print("=" * 80)

        res = mpb.run_from_config(POLY_BLOCK_CONFIG, BASIC_MODEL_DIR)

        print("[Polycrystal] Output summary:")
        for k, v in res.items():
            print(f"  - {k:15s}: {v}")
