# -*- coding: utf-8 -*-
# Project_Config.py
#
# 【项目全局配置文件】
# 功能说明：
# 1. 包含所有可修改的参数，方便用户统一管理。
# 2. 控制各个任务模块的开关和参数。

# ==============================================================================
# 1. 通用设置
# ==============================================================================
PROJECT_NAME = "AOTU_run"
RUN_NAME = None # Set to a string to force a specific run name, or None for auto-timestamp
STRUCTURE_TYPE = "FCC"

# 监控设置
MONITOR = {
    "enable_heartbeat": True,
    "heartbeat_sec": 10
}

# ==============================================================================
# 2. 材料与基础模型配置 (用于 GSFE 和 Basic 任务)
# ==============================================================================
MAT_CONFIG = {
    "element": "Cu",
    "lattice_a": 3.615,
    "lattice_c": None,  # HCP 需要
    # size: (x, y, z) in unit cells
    # 注意：FCC 周期性要求 Z 层数为 3 的倍数 (如 6, 9, 12...)
    "size": (4, 4, 6),
    "amp": 0.0,         # 随机微扰幅度
    "vacuum": 0.0,      # 真空层厚度
}

# ==============================================================================
# 3. GSFE (批量生成) 任务配置
# ==============================================================================
# 控制 Task_GSFE_Generator 的行为
GSFE_TASKS = {
    "1d_main": {
        "enable": True,
        "n_points": 21,  # 采样点数
        "n_perturb": 0
    },
    "1d_sec": {
        "enable": True,
        "n_points": 21,
        "n_perturb": 0
    },
    "2d_surf": {
        "enable": True,
        "grid_n": 6,     # 网格密度 (6x6)
        "n_perturb": 0
    },
    "key_pts": {
        "enable": True,
        "fractions": [0.0, 1.0],
        "n_perturb": 0
    }
}

# ==============================================================================
# 4. 均匀切胞 (Grid Cutting) 任务配置
# ==============================================================================
# 控制 Task_Grid_Cutter 的行为
CUTTING_CONFIG = {
    "enable": True,
    
    # 源文件目录 (相对于项目根目录)
    # 将对该目录下所有的 .vasp / POSCAR 文件进行切割
    "source_dir": "twin_complex",
    
    # 切割网格大小 (Angstrom)
    # 类似于切魔方，将整个盒子切成 grid_size 大小的子块
    # 例如: (10, 10, 10) 表示切成 10x10x10 埃的小方块
    "grid_size": (12.0, 12.0, 12.0),
    
    # 原子数限制 (min, max)
    # 只有原子数在此范围内的子块才会被保存
    "atom_limit": (10, 200),
    
    # 输出子目录名
    "output_subdir": "20_grid_cuts",
    
    # 是否添加随机微扰
    "perturb_amp": 0.0
}

# ==============================================================================
# 5. 基础晶胞 (Basic Cell) 构建配置
# ==============================================================================
# 控制 Task_Basic_Cell_Builder 的行为
BASIC_GEN_CONFIG = {
    "enable": True,
    "structure": STRUCTURE_TYPE,
    "element": MAT_CONFIG["element"],
    "lattice_a": MAT_CONFIG["lattice_a"],
    "lattice_c": MAT_CONFIG["lattice_c"],
    "supercell": (3, 3, 3),
    "vacancy_index": 1,        # 移除第几个原子制造空位
    "name_prefix": "Basic",
    "verbose": True
}

# ==============================================================================
# 6. 层错堆垛穷举 (Stacking Enumeration) 配置
# ==============================================================================
# 控制 Task_Stacking_Enumerator 的行为
STACKING_ENUM_CONFIG = {
    "enable": True,
    
    # 堆垛总层数 (必须 >= 2)
    "n_layers": 8,
    
    # 面内超胞尺寸
    "nx": 4,
    "ny": 4,
    
    # 材料参数继承自 MAT_CONFIG
    "element": MAT_CONFIG["element"],
    "lattice_a": MAT_CONFIG["lattice_a"],
    
    "output_subdir": "40_stacking_enum"
}
