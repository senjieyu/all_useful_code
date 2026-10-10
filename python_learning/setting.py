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
    # 4: 'GNT'    (梯度孪晶结构: 晶粒/晶面间距从上往下增大，取向一致)
    # 5: 'NT'     (均匀孪晶结构: 间距均匀，取向一致)
    # 6: 'all_Equal' (所有晶粒一致间距，取向随机)
    TWIN_MODE = 'GNT'
    
    # 晶向设定 (参照用户提供的 Al 孪晶多晶建模步骤)
    # X=[11-2], Y=[111], Z=[-110]
    # Y轴 ([111]) 为堆垛/孪晶生长方向
    TWIN_ORIENT = ["[11-2]", "[111]", "[-110]"]
    
    # 堆垛方向 (必须与 TWIN_ORIENT 中的 [111] 对应)
    # 用户示例中使用 Y 轴进行 duplicate 和 mirror
    STACK_AXIS = "Y"
    
    # 种子在非堆垛方向的尺寸 (单位: Unit Cells)
    # 保持为 1 以避免种子文件过大导致内存溢出
    # 对于 random 取向 (all_Equal 模式)，建议适当增大 (如 2~4) 以避免极细种子旋转时的计算错误
    # 增大至 4 以确保旋转覆盖，减少晶界真空空洞
    # 如果依然觉得不够完美，可以继续尝试增大SEED_DUP_X / Z （例如到6或8），但请注意这会显著增加内存消耗
    # 降低至 10 以平衡覆盖率和内存 (原20导致 OOM)
    SEED_DUP_X = 10
    SEED_DUP_Z = 10
    
    '''# --- 模式 1: Random (随机间距) ---'''
    # 每一层 (Matrix 或 Twin) 的厚度范围 (单位: Unit Cells)
    RANDOM_MIN_LAYERS = 2
    RANDOM_MAX_LAYERS = 8
    # 生成多少层 (Matrix + Twin + Matrix ...)
    # 注意：如果总长度过长，Atomsk 可能会在 polycrystal 步骤消耗大量内存
    # 建议保持种子总长度在 20~50 个单元左右
    RANDOM_LAYER_COUNT = 6
    
    '''# --- 模式 2: Equal (等间距) ---'''
    # 这里的 thickness 是单层厚度 (Matrix 或 Twin 的厚度)
    # 最终种子将是 [Matrix(N) + Twin(N)] 的一个周期
    EQUAL_THICKNESS = 3  # 单位: Unit Cells
    
    '''# --- 模式 3: Ratio (比例) ---'''
    RATIO_TWIN_FRACTION = 0.3
    RATIO_TOTAL_LAYERS = 30
    RATIO_MIN_THICKNESS = 2

    # --- 模式 4: GNT (梯度孪晶) ---
    # 采用"分层拼接"法 (Stacking Blocks)
    # --------------------------------------------------------------------------------------
    # 使用说明 (GNT Configuration Guide):
    # GNT_ZONES 是一个列表，每个元素代表 Y 轴上的一个区域。程序会根据 ratio 计算每个区域的高度。
    # 每个区域 (Zone) 支持以下参数：
    #
    # 1. 区域划分参数:
    #    - "ratio": float, 高度占比。例如三个区域设为 2, 3, 2，则高度分别为总高的 2/7, 3/7, 2/7。
    #    - "n_grains": int, 该区域内生成的晶粒(种子节点)数量。Atomsk将基于这些点生成Voronoi晶粒。
    #
    # 2. 内部孪晶结构参数 (决定种子类型):
    #    - "mode": string, 种子生成模式，可选 "equal", "random", "ratio", "all_Equal"
    #      (注意: 此处的 mode 仅影响该区域种子内部的层分布，所有晶粒的孪晶面始终平行于 X-Z 平面)
    #
    # 3. 模式详细参数 (根据 mode 选择对应参数):
    #    [Mode: "equal"] (等间距)
    #      - "equal_thickness": int, Matrix 和 Twin 的单层厚度 (单位: 晶胞)。
    #      - 默认 Matrix 和 Twin 厚度相同。
    #
    #    [Mode: "random"] (随机间距)
    #      - "min_layers": int, 最小层厚。
    #      - "max_layers": int, 最大层厚。
    #      - "layer_count": int, 种子中包含的总层数 (Matrix+Twin)。
    #
    #    [Mode: "ratio"] (比例控制)
    #      - "ratio_fraction": float, 孪晶体积占比 (0.0~1.0, 如 0.3)。
    #      - "total_layers": int, 种子总层数。
    #
    #    [Mode: "all_Equal"] (独立控制 Matrix/Twin 厚度)
    #      - "matrix_thickness": int, 基体层厚度。
    #      - "twin_thickness": int, 孪晶层厚度。
    #      - "layer_count": int, 种子总层数。
    # --------------------------------------------------------------------------------------
    
    GNT_ZONES = [
        # --- 区域 1 (底部) ---
        # 示例: 占比 2/7, 20个晶粒, 等间距(10层Matrix/10层Twin)
        {
            "ratio": 2.0,
            "n_grains": 20,
            "mode": "equal",
            "equal_thickness": 2,
            "random_orient_X": True,
            "random_orient_Y": True,
            "random_orient_Z": True,
            "enable_twin": True
        },
        
        # --- 区域 2 (中部) ---
        # 示例: 占比 3/7, 50个晶粒 (更密), 随机间距 (2~5层), 种子共10层
        {
            "ratio": 2.0,
            "n_grains": 20,
            "mode": "random",
            "min_layers": 2,
            "max_layers": 5,
            "layer_count": 4,
            "random_orient_X": True,
            "random_orient_Y": True,
            "random_orient_Z": True,
            "enable_twin": True
        },
        
        # --- 区域 3 (顶部) ---
        # 示例: 占比 2/7, 20个晶粒, 极细等间距(2层Matrix/2层Twin)
        {
            "ratio": 3.0,
            "n_grains": 20,
            "mode": "equal",
            "equal_thickness": 3,
            "random_orient_X": True,
            "random_orient_Y": True,
            "random_orient_Z": True,
            "enable_twin": True
        },
        
        # 如需添加更多区域 (例如 2,3,2,2)，只需在此处继续添加字典即可
        # {
        #     "ratio": 2.0,
        #     "n_grains": 30,
        #     "mode": "all_Equal",
        #     "matrix_thickness": 5,
        #     "twin_thickness": 1,
        #     "layer_count": 20,
        #     "random_orient_X": True,
        #     "random_orient_Y": True,
        #     "random_orient_Z": True,
        #     "enable_twin": True     # True=生成孪晶, False=纯晶粒(无孪晶)
        # }
    ]
    
    '''# --- 模式 5: NT (均匀孪晶) ---'''
    # 晶粒大小均匀 (网格状分布)，但可有扰动
    # 内部孪晶结构可选: 'random', 'equal', 'ratio', 'all_Equal'
    # 晶粒尺寸 [X, Y, Z] (Angstrom)，可分别设置以生成长条形或扁平晶粒
    # 恢复为均匀尺寸以生成规则的多面体/砖块
    NT_GRAIN_SIZE = [40.0, 200.0, 60.0]
    # 晶粒中心位置的随机扰动范围 (+- Angstrom)
    # 设为 0 则生成规则的长方体砖块；设为 >0 则生成不规则的多面体(Voronoi)
    NT_GRAIN_SIZE_VAR = 5.0       
    NT_INTERNAL_MODE = 'all_Equal' # 内部孪晶模式: 'random', 'equal', 'ratio', 'all_Equal'
    
    # 晶粒取向是否随机 (分别控制X/Y/Z轴的旋转)
    NT_GRAIN_ORIENT_RANDOM_X = True
    NT_GRAIN_ORIENT_RANDOM_Y = True
    NT_GRAIN_ORIENT_RANDOM_Z = True
    
    # 是否生成孪晶 (True=生成孪晶, False=生成纯晶粒)
    NT_ENABLE_TWIN = True

    # NT 模式下使用的具体参数将复用上方对应的参数:
    # - 'random'    -> 使用 RANDOM_... 参数
    # - 'equal'     -> 使用 EQUAL_... 参数
    # - 'ratio'     -> 使用 RATIO_... 参数
    # - 'all_Equal' -> 使用 ALL_EQUAL_... 参数

    '''# --- 模式 6: all_Equal (一致间距) ---'''
    # 所有晶粒内部具有相同的孪晶间距参数，但晶粒取向随机
    ALL_EQUAL_MATRIX_THICKNESS = 3  # Matrix 层厚度
    ALL_EQUAL_TWIN_THICKNESS = 3    # Twin 层厚度
    ALL_EQUAL_LAYER_COUNT = 4      # 种子层数 (降低以减少内存，20层足够覆盖大部分 Box 高度)


    
    # ==========================
    # 多晶生成参数 (Polycrystal)
    # ==========================
    # 最终多晶盒子的尺寸 (Angstrom)
    # 用户示例: 600 400 200
    # 为了演示快速运行，这里暂时设小一点，正式跑可改回 [600, 400, 200]
    POLY_BOX = [200.0, 600.0, 200.0]
    
    # 晶粒数量
    POLY_GRAINS = 30
    
    # ==========================
    # 系统控制
    # ==========================
    # 最大原子数限制
    MAX_ATOMS = 10000000
    
    # Atomsk 命令超时时间 (秒)
    ATOMSK_TIMEOUT = 720000
    
    # 输出目录
    OUT_DIR = "./output_project_v2"
