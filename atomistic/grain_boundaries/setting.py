# GB_code Configuration File
# 配置文件 - 包含所有可调整的参数

import numpy as np

# ==========================================
# 1. CSL 搜索参数 (用于 Step 1: 寻找 Sigma)
# ==========================================
# 旋转轴 (Rotation Axis)
AXIS = [1, 1, 0]

# 搜索最大 Sigma 值
LIMIT_SIGMA = 10

# ==========================================
# 2. 晶界平面搜索参数 (用于 Step 2: 寻找 GB Plane)
# ==========================================
# 选定的 Sigma 值 (从 Step 1 的结果中选择)
SIGMA = 9

# 晶体结构基组 (basis)
# 可选: 'fcc', 'bcc', 'diamond', 'sc'
BASIS = 'fcc'

# 晶界平面指数限制 (Limit for GB plane indices)
# 值越大，搜索到的高指数晶界平面越多
LIMIT_GB_PLANE = 2

# ==========================================
# 3. 晶界原子结构生成参数 (用于 Step 3: 生成 GB)
# ==========================================
# 选定的晶界平面 (GB Plane) - 从 Step 2 的结果中选择
# 必须是列表形式，例如 [1, 1, 2]
GB_PLANE = [2, -2, 1]

# 晶格常数 (Lattice Parameter) - 单位: Angstrom
# 如果 AUTO_LATTICE_PARAMETER = True，此值将被忽略，程序会根据 ATOM_TYPE 自动查找。
LATTICE_PARAMETER = 2.561

# 自动匹配晶格常数 (Auto Lattice Parameter)
# True: 根据 ATOM_TYPE 自动设置 LATTICE_PARAMETER
# False: 使用手动设置的 LATTICE_PARAMETER
AUTO_LATTICE_PARAMETER = True

# 原子类型 (Atom Type)
# 用于 VASP POSCAR 输出，以及自动查找晶格常数
# 示例: 'Al', 'Cu', 'Fe', 'Si', 'Au'
ATOM_TYPE = 'Cu'

# 原子重叠去除距离 (Overlap Distance)
# 距离小于此值 (LatP 的倍数) 的原子将被合并为一个位于中点的原子。0.0 表示不进行操作。
OVERLAP_DISTANCE = 0.5

# 从哪个晶粒移除原子 ('g1' 或 'g2')
WHICH_G = 'g1'

# 刚体平移 (Rigid Body Translation)
# 是否在晶界平面内进行刚体平移搜索低能结构 (True/False)
RIGID_TRANS = False

# 刚体平移网格划分 (仅当 RIGID_TRANS = True 时有效)
# 将 CSL 向量划分为 a * b 的网格
RIGID_A = 10
RIGID_B = 5

# 超胞尺寸 (Dimensions)
# [l1, l2, l3]
# l1: 垂直于晶界平面的方向 (需足够大以避免周期性镜像相互作用)
# l2, l3: 晶界平面内的方向
DIMENSIONS = [6, 6, 6]

# 输出文件类型 (File Type)
# 可选: 'LAMMPS' 或 'VASP'
FILE_TYPE = 'VASP'
