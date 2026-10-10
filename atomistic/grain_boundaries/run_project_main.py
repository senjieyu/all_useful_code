# -*- coding: utf-8 -*-
# Run_Project_Main.py
# 
# 【项目主入口】
# 功能说明：
# 1. 负责整个工程的配置管理（材料参数、任务开关等）。
# 2. 调度各个子模块执行任务：
#    - GSFE 生成 (Task_GSFE_Generator)
#    - 均匀切胞 (Task_Grid_Cutter)
#    - 基础晶胞构建 (Task_Basic_Cell_Builder)
#    - 层错堆垛穷举 (Task_Stacking_Enumerator)
# 3. 负责日志记录、目录创建和心跳监测。
#
# 注意事项：
# - 已集成 Project_Config.py 配置文件。
# - 包含 GSFE, GridCutter, BasicGen, StackingEnum 四大任务。

import os
import time
import threading
from datetime import datetime
from contextlib import contextmanager
import numpy as np

# 导入配置
import Project_Config as CFG

# 导入功能模块
import Lib_Structure_Builder as mc
import Task_Grid_Cutter as grid_cutter
import Task_Basic_Cell_Builder as basicgen
import Task_GSFE_Generator as batch
import Task_Stacking_Enumerator as stacking_enum  # 新增

# ======================================================================
# LOGGING & UTILS
# ======================================================================
class DualLogger:
    def __init__(self, filepath):
        self.file = open(filepath, 'w', encoding='utf-8')
        print(f"[Log] Log file: {filepath}")

    def write(self, msg):
        print(msg)
        self.file.write(msg + '\n')
        self.file.flush()

    def close(self): self.file.close()


GLOBAL_LOGGER = None


def log(msg):
    if GLOBAL_LOGGER:
        GLOBAL_LOGGER.write(msg)
    else:
        print(msg)


def log_structure_details(name, atoms, desc, params=None):
    if GLOBAL_LOGGER is None: return
    cell = atoms.get_cell()
    msg = [f"[Struct] {name}", f"  > Desc: {desc}", f"  > Atoms: {len(atoms)}",
           f"  > Cell: {np.diag(cell)}"]
    if params: msg.append(f"  > Params: {params}")
    msg.append("-" * 30)
    log("\n".join(msg))


def _now_tag(): return datetime.now().strftime("%Y%m%d_%H%M%S")


def _ensure_dir(p): os.makedirs(p, exist_ok=True); return p


class Heartbeat:
    def __init__(self, name, heartbeat_sec=10, enable=True):
        self.name, self.sec, self.enable = name, heartbeat_sec, enable
        self._stop, self._thread = threading.Event(), None

    def start(self):
        if self.enable: self._stop.clear(); self._thread = threading.Thread(target=self._run,
                                                                            daemon=True); self._thread.start()
        return self

    def _run(self):
        while not self._stop.wait(self.sec): log(f"    ... {self.name} ...")

    def stop(self):
        if self.enable: self._stop.set(); self._thread.join(1.0)


@contextmanager
def stage(name, heartbeat=False):
    log(f"\n>>> [START] {name}")
    hb = Heartbeat(name, enable=heartbeat).start()
    try:
        yield
    finally:
        hb.stop(); log(f">>> [DONE ] {name}")


def print_job_info(struct, conf, tasks, dirs):
    log("=" * 60 + f"\nJOB DASHBOARD: {struct} {conf['element']}\n" + "=" * 60)
    log(f"Output Dir: {dirs['out_root']}")


# ======================================================================
# MAIN EXECUTION
# ======================================================================
if __name__ == "__main__":
    try:
        ROOT = os.getcwd()
        OUT = _ensure_dir(os.path.join(ROOT, "outputs", CFG.PROJECT_NAME, CFG.RUN_NAME or _now_tag()))
        GLOBAL_LOGGER = DualLogger(os.path.join(OUT, "run.log"))

        # 目录准备
        D_BASE = _ensure_dir(os.path.join(OUT, "00_base_model"))
        D_BATCH = _ensure_dir(os.path.join(OUT, "10_batch"))
        D_RELAX = _ensure_dir(os.path.join(D_BATCH, "relax"))
        D_BASIC = _ensure_dir(os.path.join(OUT, "30_basic"))
        
        # 注入基础生成配置的输出目录
        CFG.BASIC_GEN_CONFIG["output_dir"] = D_BASIC
        
        # 注入 GSFE 任务通用配置
        for k in CFG.GSFE_TASKS: 
            CFG.GSFE_TASKS[k].update({
                "amp": CFG.MAT_CONFIG["amp"], 
                "element": CFG.MAT_CONFIG["element"], 
                "struct_type": CFG.STRUCTURE_TYPE
            })

        # 1. Base Slab (如果需要生成基础 GSFE 结构)
        # 如果所有 GSFE 任务都关闭，也可以跳过此步，但生成基础结构总是好的
        with stage("Build Base"):
            bf = mc.StructureFactory()
            if CFG.STRUCTURE_TYPE == "FCC":
                base, v1, v2, inf = bf.build_fcc(CFG.MAT_CONFIG["element"], CFG.MAT_CONFIG["lattice_a"], CFG.MAT_CONFIG["size"],
                                                 CFG.MAT_CONFIG["vacuum"])
            else:
                # 暂时只支持 FCC
                base, v1, v2, inf = bf.build_fcc(CFG.MAT_CONFIG["element"], CFG.MAT_CONFIG["lattice_a"], CFG.MAT_CONFIG["size"],
                                                 CFG.MAT_CONFIG["vacuum"])

            mc.AtomUtils.write_vasp_file(base, f"Base_{CFG.STRUCTURE_TYPE}", D_BASE)
            log_structure_details("Base", base, "Initial Matrix")

        # 定义目录字典供 batch 调用
        dirs = {"out_root": OUT, "out": D_BATCH, "relax": D_RELAX}
        print_job_info(CFG.STRUCTURE_TYPE, CFG.MAT_CONFIG, CFG.GSFE_TASKS, dirs)

        # 2. GSFE Batch Generation
        with stage("GSFE Batch Gen", heartbeat=True):
            if CFG.GSFE_TASKS["1d_main"]["enable"]:
                batch.generate_1d_main(base, v1, CFG.MAT_CONFIG, dirs, log_structure_details)

            if CFG.GSFE_TASKS["1d_sec"]["enable"]:
                batch.generate_1d_sec(base, v2, CFG.MAT_CONFIG, dirs, log_structure_details)

            if CFG.GSFE_TASKS["2d_surf"]["enable"]:
                batch.generate_2d_surface(base, v1, v2, CFG.MAT_CONFIG, dirs, log_structure_details)
            
            if CFG.GSFE_TASKS["key_pts"]["enable"]:
                batch.generate_key_points(base, v1, CFG.MAT_CONFIG, dirs, log_structure_details)

        # 3. Grid Cutting (均匀切胞)
        if CFG.CUTTING_CONFIG.get("enable"):
            try:
                # Inject element for fallback
                CFG.CUTTING_CONFIG["element"] = CFG.MAT_CONFIG["element"]
                with stage("Grid Cutting", heartbeat=True):
                    grid_cutter.run_from_config(CFG.CUTTING_CONFIG, OUT)
            except Exception as e:
                log(f"[GridCutter] Failed: {e}")

        # 4. Basic Cell Generation
        if CFG.BASIC_GEN_CONFIG.get("enable"):
            try:
                with stage("Basic Gen"):
                    res = basicgen.run_from_config(CFG.BASIC_GEN_CONFIG)
                    log(f"[Basic] Generated: {res['vasp_path']}")
            except Exception as e:
                log(f"[Basic] Failed: {e}")
                
        # 5. Stacking Enumeration (新增)
        if CFG.STACKING_ENUM_CONFIG.get("enable"):
            try:
                # 注入是否生成少于 n_layers 的结构的配置
                # 默认设为 False (只生成指定层数的)，或者 True (全部生成)，用户可在 main.py 修改
                CFG.STACKING_ENUM_CONFIG["generate_shorter_layers"] = True 
                
                with stage("Stacking Enumeration", heartbeat=True):
                    stacking_enum.run_from_config(CFG.STACKING_ENUM_CONFIG, OUT)
            except Exception as e:
                log(f"[StackingEnum] Failed: {e}")

        log(f"Done. Logs: {GLOBAL_LOGGER.file.name}")

    except KeyboardInterrupt:
        log("Interrupted.")
    finally:
        if GLOBAL_LOGGER: GLOBAL_LOGGER.close()
