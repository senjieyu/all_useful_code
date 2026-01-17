# -*- coding: utf-8 -*-
"""
自动生成 (111) GSFE / γ-surface / 孪晶 等相关结构的 VASP POSCAR (.vasp) 文件

依赖：
    pip install ase

注意：
1. 请根据体系修改：
   - element
   - a (晶格常数)
   - size_111 (三维尺寸)
2. 真正的弛豫在 DFT 中完成，这里只负责生成初始结构。
"""

import os
import numpy as np

from ase.build import fcc111, stack
from ase.io.vasp import write_vasp

# ==================== 全局参数 ====================

element = "Cu"        # 体系元素，自行修改
a = 3.615             # 晶格常数（单位 Å），自行修改
size_111 = (3, 3, 12) # (nx, ny, nlayers)：可以先用 (3,3,12) 调流程，再放大

output_root = "C:/Users/senji/Desktop/super/GSFE（广义层错曲线）"
os.makedirs(output_root, exist_ok=True)


# 小工具：加随机微扰
def add_random_perturbation(atoms, amp=0.02):
    """给原子坐标添加均匀分布的随机微扰，幅度 amp（单位 Å）"""
    perturbed = atoms.copy()
    disp = (np.random.rand(len(perturbed), 3) - 0.5) * 2.0 * amp
    perturbed.positions += disp
    return perturbed


# 小工具：写 vasp 文件，规范命名
def write_structure(atoms, name):
    """写 VASP POSCAR 文件，名字为 name.vasp"""
    path = os.path.join(output_root, f"{name}.vasp")
    # 不写注释，以便你自己在后处理里加
    write_vasp(path, atoms, direct=True, sort=True, vasp5=True)
    print("Wrote:", path)


# 小工具：构造 (111) 取向的“块体超胞”（无真空，三向 PBC）
def build_111_bulk_like():
    # vacuum=0，z 方向无真空；三向 PBC
    atoms = fcc111(
        symbol=element,
        size=size_111,
        a=a,
        vacuum=0.0,
        orthogonal=False  # 标准六方胞
    )
    atoms.pbc = [True, True, True]
    return atoms


# 小工具：把原子按 z 坐标切成上下两部分
def split_upper_lower(atoms):
    z = atoms.positions[:, 2]
    z_mid = 0.5 * (z.min() + z.max())
    lower = z < z_mid + 1e-5
    upper = ~lower
    return lower, upper, z_mid


# ==================== 1. 1D GSFE 主路径：沿 dir1 (近似 ⟨11-2⟩) ====================

def generate_1d_gsfe_main(n_points=16, n_perturb=4, max_shift_frac=1.0):
    """
    1D GSFE 主路径（最关键）
    - 沿 (111) 面内的 dir1 从 0 → max_shift_frac * |dir1| 的滑移
    - n_points：总采样点数（建议 12–20）
    - n_perturb：每个点微扰数（建议 3–5）
    """
    base = build_111_bulk_like()
    cell = base.get_cell()
    dir1 = cell[0]        # 把 cell[0] 视作 Shockley 方向（近似 ⟨11-2⟩）
    lower_mask, upper_mask, _ = split_upper_lower(base)

    for i in range(n_points):
        frac = max_shift_frac * i / n_points  # 0 ~ max_shift_frac
        shift_vec = frac * dir1
        atoms = base.copy()
        positions = atoms.positions
        positions[upper_mask] += shift_vec
        atoms.set_positions(positions)

        # 基础结构
        name = f"1D_main_step{i:02d}"
        write_structure(atoms, name)

        # 微扰版本
        for k in range(n_perturb):
            pert = add_random_perturbation(atoms, amp=0.02)
            name_p = f"1D_main_step{i:02d}_pert{k:02d}"
            write_structure(pert, name_p)


# ==================== 2. 1D GSFE 次路径：沿 dir2 (近似 ⟨1-10⟩) ====================

def generate_1d_gsfe_secondary(n_points=12, n_perturb=3, max_shift_frac=1.0):
    """
    1D GSFE 次路径（强烈推荐）
    - 沿 (111) 面内的 dir2 从 0 → max_shift_frac * |dir2|
    """
    base = build_111_bulk_like()
    cell = base.get_cell()
    dir2 = cell[1]        # 把 cell[1] 视作另一剪切方向（近似 ⟨1-10⟩）
    lower_mask, upper_mask, _ = split_upper_lower(base)

    for i in range(n_points):
        frac = max_shift_frac * i / n_points
        shift_vec = frac * dir2
        atoms = base.copy()
        positions = atoms.positions
        positions[upper_mask] += shift_vec
        atoms.set_positions(positions)

        name = f"1D_sec_step{i:02d}"
        write_structure(atoms, name)

        for k in range(n_perturb):
            pert = add_random_perturbation(atoms, amp=0.02)
            name_p = f"1D_sec_step{i:02d}_pert{k:02d}"
            write_structure(pert, name_p)


# ==================== 3. 2D γ-surface：dir1 × dir2 网格 ====================

def generate_2d_gamma_surface(grid_n=5, n_perturb=2, max_shift_frac=1.0):
    """
    2D γ-surface（稀疏网格即可）
    - 在 (dir1, dir2) 平面上做 grid_n × grid_n 网格
    - grid_n=5 → 5×5 共 25 个点；grid_n=7 → 49 个点
    """
    base = build_111_bulk_like()
    cell = base.get_cell()
    dir1, dir2 = cell[0], cell[1]
    lower_mask, upper_mask, _ = split_upper_lower(base)

    for ix in range(grid_n):
        for iy in range(grid_n):
            fx = max_shift_frac * ix / grid_n
            fy = max_shift_frac * iy / grid_n
            shift_vec = fx * dir1 + fy * dir2

            atoms = base.copy()
            positions = atoms.positions
            positions[upper_mask] += shift_vec
            atoms.set_positions(positions)

            name = f"2D_gamma_x{ix:02d}_y{iy:02d}"
            write_structure(atoms, name)

            for k in range(n_perturb):
                pert = add_random_perturbation(atoms, amp=0.02)
                name_p = f"2D_gamma_x{ix:02d}_y{iy:02d}_pert{k:02d}"
                write_structure(pert, name_p)


# ==================== 4. 法向放松版本（关键点：0, USF, ISF, 周期） ====================

def generate_normal_relax_keypoints(key_fracs=(0.0, 0.5, 0.8, 1.0),
                                    n_perturb=4):
    """
    法向放松版本（可选但很值）
    - 这里只生成“初始几何 + 标记为 key_*”，真正的 out-of-plane 弛豫放在 DFT 中做。
    - key_fracs：沿主路径（dir1）的关键滑移分数（0：完美；~USF；~ISF；~周期）
    """
    base = build_111_bulk_like()
    cell = base.get_cell()
    dir1 = cell[0]
    lower_mask, upper_mask, _ = split_upper_lower(base)

    for idx, frac in enumerate(key_fracs):
        shift_vec = frac * dir1
        atoms = base.copy()
        positions = atoms.positions
        positions[upper_mask] += shift_vec
        atoms.set_positions(positions)

        # 这里不改 cell，只是标记为 "normal_relax"，方便在 DFT 中允许 z 方向自由弛豫
        name = f"normal_relax_key{idx:02d}_frac{frac:.2f}"
        write_structure(atoms, name)

        for k in range(n_perturb):
            pert = add_random_perturbation(atoms, amp=0.02)
            name_p = f"{name}_pert{k:02d}"
            write_structure(pert, name_p)


# ==================== 5. 应力/体积敏感性：±1–2% 体积或法向应变 ====================

def generate_strain_sensitivity(strain_list=( -0.02, -0.01, 0.01, 0.02 ),
                                key_fracs=(0.0, 0.5, 0.8),
                                n_perturb=3):
    """
    应力/体积敏感性（进阶）
    - 对若干关键滑移点 key_fracs（沿 dir1）在 ±1–2% 应变下重复生成
    - 这里示例只对 cell 等比例缩放（近似各向同性体积/三轴应变）
      你可以改成只对法向 (z) 拉伸收缩。
    """
    base0 = build_111_bulk_like()
    cell0 = base0.get_cell()
    dir1 = cell0[0]
    lower_mask, upper_mask, _ = split_upper_lower(base0)

    for eps in strain_list:
        # 等比例缩放 cell + positions（体积/三轴应变）
        scale = 1.0 + eps
        cell_eps = cell0 * scale

        for f in key_fracs:
            atoms = base0.copy()
            positions = atoms.positions
            shift_vec = f * dir1
            positions[upper_mask] += shift_vec

            # 统一缩放 positions
            positions *= scale
            atoms.set_positions(positions)
            atoms.set_cell(cell_eps, scale_atoms=False)

            name = f"strain_eps{eps:+.3f}_frac{f:.2f}"
            write_structure(atoms, name)

            for k in range(n_perturb):
                pert = add_random_perturbation(atoms, amp=0.02)
                name_p = f"{name}_pert{k:02d}"
                write_structure(pert, name_p)


# ==================== 6. 孪晶锚点：CTB 结构 + 微扰 ====================

def generate_twin_ctb(n_perturb_per=10):
    """
    孪晶锚点（建议必加，但不必很多）
    - 这里用一个比较简单的构造方法：
        1) 生成 (111) slab
        2) 复制一个镜像 slab（关于中间面反转）
        3) 用 stack 把两个 slab 堆叠，形成一个近似 CTB 的结构
    - 具体堆垛顺序/层数可根据体系精细调整，这里给出一个自动化框架。
    """
    # 先做一个 (111) slab，z 方向一半层数
    nx, ny, nlayers = size_111
    slab_layers = max(4, nlayers // 2)

    slab = fcc111(symbol=element,
                  size=(nx, ny, slab_layers),
                  a=a,
                  vacuum=0.0,
                  orthogonal=False)
    slab.pbc = [True, True, True]

    # 构造一个镜像 slab：对 z 取负，再平移回去
    twin = slab.copy()
    pos = twin.positions
    z_mid = 0.5 * (pos[:, 2].min() + pos[:, 2].max())
    pos[:, 2] = 2 * z_mid - pos[:, 2]
    twin.set_positions(pos)

    # 用 stack 把 slab 与镜像 slab 堆叠在 z 方向，界面处形成 CTB 近似
    ct_win = stack(slab, twin, axis=2, distance=0.0)
    ct_win.pbc = [True, True, True]

    # 写 CTB 基础结构
    name_base = "twin_CTB_base"
    write_structure(ct_win, name_base)

    # 微扰版本
    for k in range(n_perturb_per):
        pert = add_random_perturbation(ct_win, amp=0.02)
        name_p = f"{name_base}_pert{k:02d}"
        write_structure(pert, name_p)


# ==================== main：调用所有模块 ====================

if __name__ == "__main__":
    print("Generating 1D GSFE main path ...")
    generate_1d_gsfe_main(n_points=16, n_perturb=4)

    print("Generating 1D GSFE secondary path ...")
    generate_1d_gsfe_secondary(n_points=12, n_perturb=3)

    print("Generating 2D gamma surface ...")
    generate_2d_gamma_surface(grid_n=5, n_perturb=2)

    print("Generating normal-relax keypoints (geometry only) ...")
    generate_normal_relax_keypoints(key_fracs=(0.0, 0.5, 0.8, 1.0), n_perturb=4)

    print("Generating strain-sensitivity structures ...")
    generate_strain_sensitivity(
        strain_list=(-0.02, -0.01, 0.01, 0.02),
        key_fracs=(0.0, 0.5, 0.8),
        n_perturb=3
    )

    print("Generating twin CTB anchors ...")
    generate_twin_ctb(n_perturb_per=10)

    print("All structures written to:", output_root)
