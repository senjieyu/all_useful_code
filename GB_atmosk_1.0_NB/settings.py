# -*- coding: utf-8 -*-
"""
settings.py
只需要改这里即可。build_atomsk.py 会全自动创建文件夹并生成结构。
"""

# Atomsk 可执行文件（若不在 PATH，请填绝对路径）
ATOMSK = "atomsk"

# 材料参数（示例：Cu fcc）
MATERIAL = {
    "lattice": "fcc",
    "a": 3.615,       # Å
    "element": "Cu",
}

# 输出目录规划：最终结构 vs 过程文件
FINAL_DIR = "final_structures"  # 最终 .cfg 汇总目录
WORK_DIR = "work"               # 过程文件根目录，每个 job 一个子文件夹

# merge 后选项
WRAP = True
SHIFT_SYSID2 = (0.0, 0.0, 0.0)  # 只对 sysID=2 的晶粒平移（Å），γ-surface 扫描可改

# ---------------------------
# 选择要自动生成的任务（开关）
# ---------------------------
MAKE_TILT_210_001 = True          # {210}[001] 对称倾转（orient + merge）
MAKE_TWIST_001 = True             # [001] 对称扭转（orient + merge）
MAKE_ARBITRARY_TILT = True        # 任意角倾转（rotate Z + orthogonal-cell + merge）
MAKE_ARBITRARY_TWIST = True       # 任意角扭转（rotate Y + orthogonal-cell + merge）
MAKE_POLYCRYSTAL_2G = True        # polycrystal 两晶粒（你的 polyX 思路）

# ---------------------------
# 尺寸控制：duplicate 倍数
# ---------------------------
DUPLICATE_TILT = (1, 5, 1)
DUPLICATE_TWIST = (1, 5, 1)
DUPLICATE_ARBITRARY_TILT = (1, 1, 1)
DUPLICATE_ARBITRARY_TWIST = (1, 10, 1)

# 任意角度（单晶粒旋转 ±alpha，总失配=2alpha）
ARBITRARY_ANGLE_DEG = 5.1

# ---------------------------
# polycrystal 参数
# ---------------------------
POLY_BOX = (100.0, 100.0, 10.0)  # Å
POLY_ROT_DEG = 26.57             # 绕 Z 轴分别 ±角度
POLY_FILE_NAME = "polyX.txt"
