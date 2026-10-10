# -*- coding: utf-8 -*-
# main.py (NO pause / NO input; stage monitor + heartbeat; safe Ctrl+C handling)
import os
import time
import threading
from datetime import datetime
from contextlib import contextmanager

import model_construct as mc
import cons_batch_generate as batch

import make_poly_and_blocks as mpb
import basic as basicgen


# ======================================================================
# ========================= 1) USER CONFIG =============================
# 说明：使用者建议只改这一块
# ======================================================================

PROJECT_NAME = "AOTU_run"
RUN_NAME = None  # None -> 自动时间戳

# 监控输出参数
MONITOR = {
    "enable_heartbeat": True,   # True：长任务时每隔一段时间输出“还在跑”
    "heartbeat_sec": 10,        # 心跳间隔（秒）
}

STRUCTURE_TYPE = "FCC"  # "FCC" / "BCC" / "HCP"

MAT_CONFIG = {
    "element": "Cu",
    "lattice_a": 3.615,
    "lattice_c": None,        # HCP 才用；None -> c=1.633*a
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
    "twin": {"enable": True, "n_perturb": 0},
}

POLY_BLOCK_CONFIG = {
    "enable": True,
    "element": MAT_CONFIG["element"],
    "structure": STRUCTURE_TYPE,
    "lattice_a": MAT_CONFIG["lattice_a"],

    "use_twin_unit": True,
    "twin_orient": ("[11-2]", "[111]", "[-110]"),
    "twin_duplicate": (1, 3, 1),

    "poly_box": (10, 10, 10),
    "poly_grains": 2,
    "max_poly_atoms": 300000,

    "block_atom_range": (100, 300),
    "block_window_min": (12, 12, 12),
    "block_window_max": (20, 20, 20),
    "n_blocks": 40,

    "seed": 42,
    "block_perturb_amp": 0.0,
}

# basic：用你那份 basic.py（超胞 -> del atoms[idx] -> 写出 -> 读回校验）
BASIC_GEN_CONFIG = {
    "enable": True,
    "structure": STRUCTURE_TYPE,
    "element": MAT_CONFIG["element"],
    "lattice_a": MAT_CONFIG["lattice_a"],
    "lattice_c": MAT_CONFIG["lattice_c"],
    "supercell": (3, 3, 3),
    "vacancy_index": 1,
    "output_dir": None,  # main 注入
    "name_prefix": "Basic",
    "verbose": True,
}


# ======================================================================
# ========================= 2) UTILITIES ===============================
# ======================================================================

def _now_tag():
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)
    return p


def log(msg: str):
    print(msg, flush=True)


class Heartbeat:
    def __init__(self, name: str, heartbeat_sec: int = 10, enable: bool = True):
        self.name = name
        self.heartbeat_sec = max(1, int(heartbeat_sec))
        self.enable = bool(enable)
        self._stop_evt = threading.Event()
        self._thread = None

    def start(self):
        if not self.enable:
            return self
        self._stop_evt.clear()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def _run(self):
        t0 = time.time()
        while not self._stop_evt.is_set():
            time.sleep(self.heartbeat_sec)
            if self._stop_evt.is_set():
                break
            dt = time.time() - t0
            log(f"    ... [Heartbeat] {self.name} running ({dt:.1f}s elapsed) ...")

    def stop(self):
        if not self.enable:
            return
        self._stop_evt.set()
        if self._thread is not None:
            self._thread.join(timeout=1.0)


@contextmanager
def stage(name: str, heartbeat: bool = False):
    log(f"\n>>> [START] {name}")
    hb = None
    t0 = time.time()
    try:
        if heartbeat and MONITOR.get("enable_heartbeat", True):
            hb = Heartbeat(name, MONITOR.get("heartbeat_sec", 10), enable=True).start()
        yield
        dt = time.time() - t0
        log(f">>> [DONE ] {name} ({dt:.2f}s)")
    except Exception as e:
        dt = time.time() - t0
        log(f">>> [FAIL ] {name} ({dt:.2f}s) : {e}")
        raise
    finally:
        if hb is not None:
            hb.stop()


def format_vec(v):
    return f"[{v[0]:6.3f}, {v[1]:6.3f}, {v[2]:6.3f}]"


def print_job_info(struct_type, mat_conf, tasks, dirs, v1, v2, info):
    log("\n" + "=" * 80)
    log(f"{'STRUCTURAL GENERATION JOB DASHBOARD':^80}")
    log("=" * 80)

    log("\n[Basic Information]")
    log(f"  - Element:       {mat_conf['element']}")
    log(f"  - Phase Struct:  {struct_type}")
    log(f"  - Lattice:       a={mat_conf['lattice_a']}")
    log(f"  - Supercell:     {mat_conf['size']}")
    log(f"  - Vacuum:        {mat_conf['vacuum']} Angstrom")

    log("\n[Slip Systems & Crystallography]")
    log("  ------------------------------------------------------------")
    log("  * Main Slip Path (v1):")
    log(f"      Direction : {info.get('main_label', 'N/A')}")
    log(f"      Type      : {info.get('main_desc', 'N/A')}")
    log(f"      Vector    : {format_vec(v1)}")
    log("  ------------------------------------------------------------")
    log("  * Secondary Slip Path (v2):")
    log(f"      Direction : {info.get('sec_label', 'N/A')}")
    log(f"      Type      : {info.get('sec_desc', 'N/A')}")
    log(f"      Vector    : {format_vec(v2)}")
    log("  ------------------------------------------------------------")

    log("\n[Output Directories]")
    for k, v in dirs.items():
        log(f"  - {k:<12}: {os.path.abspath(v)}")

    log("\n[Task Switches]")
    for k, cfg in tasks.items():
        log(f"  - {k:<10}: {'ON' if cfg.get('enable', False) else 'OFF'}")

    log("=" * 80 + "\n")


# ======================================================================
# =============================== MAIN ==================================
# ======================================================================
if __name__ == "__main__":

    try:
        # output dirs
        run_tag = RUN_NAME or _now_tag()
        ROOT_DIR = os.getcwd()

        OUT_ROOT = _ensure_dir(os.path.join(ROOT_DIR, "outputs", PROJECT_NAME, run_tag))
        DIR_BASE = _ensure_dir(os.path.join(OUT_ROOT, "00_base_model"))
        DIR_BATCH = _ensure_dir(os.path.join(OUT_ROOT, "10_batch"))
        DIR_RELAX = _ensure_dir(os.path.join(DIR_BATCH, "to_relax"))
        DIR_POLY = _ensure_dir(os.path.join(OUT_ROOT, "20_poly_block"))
        DIR_BASIC = _ensure_dir(os.path.join(OUT_ROOT, "30_basic"))

        dirs = {
            "out_root": OUT_ROOT,
            "base_model": DIR_BASE,
            "batch_out": DIR_BATCH,
            "batch_relax": DIR_RELAX,
            "poly_block": DIR_POLY,
            "basic": DIR_BASIC,
        }

        BASIC_GEN_CONFIG["output_dir"] = DIR_BASIC

        # inject batch common fields
        for k in TASKS:
            TASKS[k]["amp"] = MAT_CONFIG["amp"]
            TASKS[k]["element"] = MAT_CONFIG["element"]
            TASKS[k]["struct_type"] = STRUCTURE_TYPE

        # 1) build base slab
        with stage("BUILD BASE SLAB"):
            log("Initializing Model Construction...")
            builder = mc.StructureFactory()

            if STRUCTURE_TYPE == "FCC":
                base_atoms, v_main, v_sec, info = builder.build_fcc(
                    MAT_CONFIG["element"],
                    MAT_CONFIG["lattice_a"],
                    MAT_CONFIG["size"],
                    vacuum=MAT_CONFIG["vacuum"],
                )
            elif STRUCTURE_TYPE == "BCC":
                base_atoms, v_main, v_sec, info = builder.build_bcc(
                    MAT_CONFIG["element"],
                    MAT_CONFIG["lattice_a"],
                    MAT_CONFIG["size"],
                    vacuum=MAT_CONFIG["vacuum"],
                )
            elif STRUCTURE_TYPE == "HCP":
                c_val = MAT_CONFIG["lattice_c"] if MAT_CONFIG["lattice_c"] else 1.633 * MAT_CONFIG["lattice_a"]
                base_atoms, v_main, v_sec, info = builder.build_hcp(
                    MAT_CONFIG["element"],
                    MAT_CONFIG["lattice_a"],
                    c_val,
                    MAT_CONFIG["size"],
                    vacuum=MAT_CONFIG["vacuum"],
                )
            else:
                raise ValueError(f"Unknown STRUCTURE_TYPE: {STRUCTURE_TYPE}")

            mc.AtomUtils.write_vasp_file(
                base_atoms,
                f"Base_{STRUCTURE_TYPE}_{MAT_CONFIG['element']}",
                DIR_BASE
            )

        print_job_info(STRUCTURE_TYPE, MAT_CONFIG, TASKS, dirs, v_main, v_sec, info)

        # 2) batch tasks
        batch_dirs = {"base_dir": DIR_BASE, "out": DIR_BATCH, "relax": DIR_RELAX}
        total_gen = 0

        with stage("BATCH: 1d_main", heartbeat=True):
            if TASKS["1d_main"]["enable"]:
                total_gen += batch.generate_1d_main(base_atoms, v_main, TASKS["1d_main"], batch_dirs)

        with stage("BATCH: 1d_sec", heartbeat=True):
            if TASKS["1d_sec"]["enable"]:
                total_gen += batch.generate_1d_sec(base_atoms, v_sec, TASKS["1d_sec"], batch_dirs)

        with stage("BATCH: 2d_surf", heartbeat=True):
            if TASKS["2d_surf"]["enable"]:
                raw_v1 = base_atoms.get_cell()[0] / MAT_CONFIG["size"][0]
                raw_v2 = base_atoms.get_cell()[1] / MAT_CONFIG["size"][1]
                total_gen += batch.generate_2d_surface(base_atoms, raw_v1, raw_v2, TASKS["2d_surf"], batch_dirs)

        with stage("BATCH: key_pts", heartbeat=True):
            if TASKS["key_pts"]["enable"]:
                total_gen += batch.generate_key_points(base_atoms, v_main, TASKS["key_pts"], batch_dirs)

        with stage("BATCH: strain", heartbeat=True):
            if TASKS["strain"]["enable"]:
                total_gen += batch.generate_strain(base_atoms, v_main, TASKS["strain"], batch_dirs)

        with stage("BATCH: twin", heartbeat=True):
            if TASKS["twin"]["enable"]:
                total_gen += batch.generate_twin(base_atoms, TASKS["twin"], batch_dirs)

        log(f"\n[Batch] Done. Total files: {total_gen}")

        # 3) poly+block
        with stage("POLY+BLOCK", heartbeat=True):
            if POLY_BLOCK_CONFIG.get("enable", False):
                log("[Polycrystal] Generating Voronoi polycrystal + blocks using Atomsk ...")
                res = mpb.run_from_config(POLY_BLOCK_CONFIG, DIR_POLY)
                log("[Polycrystal] Output summary:")
                for k, v in res.items():
                    log(f"  - {k:15s}: {v}")
            else:
                log("[Polycrystal] Skipped (enable=False)")

        # 4) basic
        with stage("BASIC: supercell + remove ONE atom", heartbeat=False):
            if BASIC_GEN_CONFIG.get("enable", False):
                out = basicgen.run_from_config(BASIC_GEN_CONFIG)
                log("[Basic] Output summary:")
                for k, v in out.items():
                    log(f"  - {k:15s}: {v}")
            else:
                log("[Basic] Skipped (enable=False)")

        log("\nAll done.")
        log(f"Outputs saved under: {os.path.abspath(OUT_ROOT)}")

    except KeyboardInterrupt:
        # 你按 Ctrl+C 时，会走这里，不会再出现“像代码错了一样”的 forrtl 信息误导
        log("\n[Interrupted] User pressed Ctrl+C. Exiting safely.")
