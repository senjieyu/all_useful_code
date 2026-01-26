# -*- coding: utf-8 -*-
# setting.py

class Config:
    # ==========================
    # 基础材料参数
    # ==========================
    ELEMENT = "Cu"
    LATTICE_A = 3.615
    
    # ==========================
    # 孪晶生成参数 (Twin Configuration)
    # ==========================
    # 孪晶模式: 
    # 1: 'random' (随机间距)
    # 2: 'equal'  (等间距)
    # 3: 'ratio'  (根据比例随机间距)
    TWIN_MODE = 'equal'
    
    # 晶向设定 (参照用户提供的 Al 孪晶多晶建模步骤)
    # X=[11-2], Y=[111], Z=[-110]
    # Y轴 ([111]) 为堆垛/孪晶生长方向
    TWIN_ORIENT = ["[11-2]", "[111]", "[-110]"]
    
    # 堆垛方向 (必须与 TWIN_ORIENT 中的 [111] 对应)
    # 用户示例中使用 Y 轴进行 duplicate 和 mirror
    STACK_AXIS = "Y"
    
    # 种子在非堆垛方向的尺寸 (单位: Unit Cells)
    # 保持为 1 以避免种子文件过大导致内存溢出
    SEED_DUP_X = 1
    SEED_DUP_Z = 1
    
    # --- 模式 1: Random (随机间距) ---
    # 每一层 (Matrix 或 Twin) 的厚度范围 (单位: Unit Cells)
    RANDOM_MIN_LAYERS = 2
    RANDOM_MAX_LAYERS = 8
    # 生成多少层 (Matrix + Twin + Matrix ...)
    # 注意：如果总长度过长，Atomsk 可能会在 polycrystal 步骤消耗大量内存
    # 建议保持种子总长度在 20~50 个单元左右
    RANDOM_LAYER_COUNT = 6
    
    # --- 模式 2: Equal (等间距) ---
    # 这里的 thickness 是单层厚度 (Matrix 或 Twin 的厚度)
    # 最终种子将是 [Matrix(N) + Twin(N)] 的一个周期
    EQUAL_THICKNESS = 3  # 单位: Unit Cells
    
    # --- 模式 3: Ratio (比例) ---
    RATIO_TWIN_FRACTION = 0.3
    RATIO_TOTAL_LAYERS = 30
    RATIO_MIN_THICKNESS = 2
    
    # ==========================
    # 多晶生成参数 (Polycrystal)
    # ==========================
    # 最终多晶盒子的尺寸 (Angstrom)
    # 用户示例: 600 400 200
    # 为了演示快速运行，这里暂时设小一点，正式跑可改回 [600, 400, 200]
    POLY_BOX = [200.0, 200.0, 200.0]
    
    # 晶粒数量
    POLY_GRAINS = 10
    
    # ==========================
    # 系统控制
    # ==========================
    # 最大原子数限制
    MAX_ATOMS = 5000000
    
    # Atomsk 命令超时时间 (秒)
    ATOMSK_TIMEOUT = 7200
    
    # 输出目录
    OUT_DIR = "./output_project_v2"
