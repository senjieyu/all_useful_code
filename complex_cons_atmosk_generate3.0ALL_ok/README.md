# Complex Twin Polycrystal Generator

这是一个基于 **Atomsk** 的自动化脚本工具，用于生成包含复杂内部孪晶结构的多晶模型。该工具专为模拟高熵合金或金属材料（如 Cu, Al）的微观结构而设计，支持高度自定义的孪晶间距控制，特别是**梯度孪晶 (GNT)** 和 **纳米孪晶 (NT)** 结构。

## 1. 核心功能

*   **自动化多晶生成**：基于 Voronoi 泰森多边形法生成多晶结构。
*   **复杂内部孪晶**：每个晶粒内部不再是单一晶向，而是包含平行排列的孪晶层（Twin Lamellae）。
*   **GNT 无缝拼接**：采用“幽灵节点 (Ghost Node)”策略，实现梯度区域之间的**无缝自然过渡**，晶界可跨越区域边界，消除人工拼接痕迹。
*   **3D 随机取向控制**：支持对 X、Y、Z 轴分别进行随机旋转控制，可生成织构或完全随机取向的多晶。
*   **灵活的间距控制**：支持多种孪晶生成模式：
    1.  **GNT (梯度孪晶)**：沿 Y 轴分区域控制孪晶密度（如：底部粗晶 -> 顶部细晶），区域间无缝连接。
    2.  **NT (均匀孪晶)**：所有晶粒内部孪晶间距均匀。
    3.  **Random (随机)**：孪晶层和基体层的厚度在指定范围内随机变化。
    4.  **Equal (等间距)**：生成规则的、等厚度的孪晶结构。
    5.  **Ratio (比例)**：根据设定的体积百分比（如 30% 孪晶）自动分配层厚度。
    6.  **all_Equal (一致间距)**：所有晶粒内部孪晶间距相同（可分别设置 Matrix/Twin 厚度），但晶粒取向随机。
*   **自动时间戳存档**：每次运行会自动创建一个带时间戳的新文件夹，防止覆盖旧数据。
*   **内存优化**：针对 Atomsk 的内存限制进行了优化，种子构建时仅在非堆垛方向保留最小尺寸，减少内存占用。

## 2. 建模原理 (GNT 模式)

为了生成具有梯度结构且边界自然的多晶，本工具采用了 **Multi-Pass Ghost Node (多重扫描幽灵节点)** 策略：

1.  **全局节点规划**：首先在整个模拟框内生成所有区域的种子点（Nodes）。
2.  **多重扫描 (Multi-Pass)**：
    *   **Pass 1 (Zone 1)**：使用 Zone 1 的孪晶参数（如层厚 10nm）生成**全尺寸**多晶，但只保留属于 Zone 1 区域的原子。其他区域的晶粒作为“幽灵”存在，仅用于维持晶界形状。
    *   **Pass 2 (Zone 2)**：使用 Zone 2 的参数（如层厚 2nm）再次生成全尺寸多晶，提取属于 Zone 2 的原子。
    *   ...以此类推。
3.  **合并 (Merge)**：将各次扫描提取的原子完美拼合。由于所有 Pass 使用同一套 Voronoi 节点，晶界在空间上是连续且锁合的。

## 3. 使用方法

### 3.1 环境要求
*   Python 3.x
*   **Atomsk** (必须已安装并添加到系统 PATH 环境变量中)

### 3.2 运行
直接运行 `main.py`：

```bash
python main.py
```

程序运行结束后，会在当前目录下生成一个新的文件夹（例如 `output_project_v2_GNT_20260123_123000`），其中包含：
*   `poly.lmp`: 最终的 LAMMPS 数据文件。
*   `merged.cfg`: 中间合并文件（保留了 GrainID 信息）。
*   `log_*.txt`: 各个步骤的运行日志。

### 3.3 参数配置 (`setting.py`)

所有参数均在 `setting.py` 中修改。

#### 3.3.1 全局参数

```python
ELEMENT = "Cu"
LATTICE_A = 3.615
POLY_BOX = [200.0, 600.0, 200.0]  # 模拟盒尺寸 [X, Y, Z]
POLY_GRAINS = 30                  # 默认总晶粒数 (非 GNT 模式下使用)
TWIN_MODE = 'GNT'                 # 模式选择: 'GNT', 'NT', 'random', 'equal', 'ratio', 'all_Equal'
```

#### 3.3.2 GNT 模式配置 (Gradient Nano-Twins)

```python
GNT_ZONES = [
    # --- 区域 1 (底部) ---
    {
        "ratio": 2.0,           # 高度占比
        "n_grains": 20,         # 该区域内的晶粒数量
        "mode": "equal",        # 孪晶模式 (equal/random/ratio)
        "equal_thickness": 10,  # 孪晶层厚 (晶胞数)
        "random_orient_X": True, # X轴随机旋转 (Tilted Twin Planes)
        "random_orient_Y": True, # Y轴随机旋转 (Rotation around Stacking Axis)
        "random_orient_Z": True, # Z轴随机旋转 (Tilted Twin Planes)
        "enable_twin": True     # 是否生成孪晶 (False=纯晶粒)
    },
    
    # --- 区域 2 (中部) ---
    {
        "ratio": 3.0,
        "n_grains": 100,        # 更高的晶粒密度
        "mode": "random",       # 随机层厚
        "min_layers": 2,
        "max_layers": 5,
        "layer_count": 10,
        "random_orient_X": True,
        "random_orient_Y": True,
        "random_orient_Z": True,
        "enable_twin": True
    },
    # ... 更多区域
]
```

#### 3.3.3 取向控制详解

*   **`random_orient_Y`**: 绕堆垛轴（孪晶法线）旋转。
    *   `True` + `X/Z=False`: 孪晶面保持水平平行，但在平面内随机旋转（类似竹节结构）。
*   **`random_orient_X` / `random_orient_Z`**: 绕水平轴旋转。
    *   `True`: 孪晶面会发生倾斜。
*   **全 `True`**: 完全随机取向的 3D 多晶，每个晶粒的孪晶面朝向各不相同。

#### 3.3.4 其他模式参数

*   **NT 模式 (Uniform Nano-Twins)**
    *   `NT_GRAIN_SIZE`: 控制均匀晶粒的近似尺寸。
    *   `NT_GRAIN_ORIENT_RANDOM_X/Y/Z`: 全局随机取向开关。
    *   `NT_INTERNAL_MODE`: 内部孪晶模式（复用下方参数）。

*   **Random (随机间距)**
    *   `RANDOM_MIN_LAYERS` / `RANDOM_MAX_LAYERS`: 随机层厚范围。
    *   `RANDOM_LAYER_COUNT`: 种子层数。

*   **Equal (等间距)**
    *   `EQUAL_THICKNESS`: Matrix 和 Twin 的固定层厚。

*   **Ratio (比例)**
    *   `RATIO_TWIN_FRACTION`: 孪晶相的体积占比 (0.0~1.0)。

*   **all_Equal (一致间距)**
    *   `ALL_EQUAL_MATRIX_THICKNESS`: Matrix 层厚。
    *   `ALL_EQUAL_TWIN_THICKNESS`: Twin 层厚。

## 4. 常见问题 (FAQ)

**Q: 为什么 GNT 模式运行很慢？**
A: 为了保证晶界自然过渡，程序实际上计算了 `N` 次完整的全尺寸多晶（N=区域数）。如果你有 3 个区域，计算量大约是普通多晶的 3 倍。这是为了获得高质量模型所付出的必要代价。

**Q: 为什么晶粒看起来大小差不多？**
A: 晶粒的视觉直径与数量的**立方根**成反比。要看到明显的“大小梯度”，建议数量级按指数增加（例如 20 -> 200 -> 2000）。

**Q: 为什么 enable_twin=False 还是有线条？**
A: 那些线条是**晶界 (Grain Boundaries)**，不是孪晶界。如果你启用了随机取向，晶粒之间会有取向差，显示为边界。如果想要完全单晶，需同时设置 `n_grains=1` 或 `random_orient=False`。
