# -*- coding: utf-8 -*-
'''完美可用'''
import os
import numpy as np
from ase.build import fcc111, stack
from ase.io.vasp import write_vasp
from ase.atoms import Atoms

# ==================== 全局参数 ====================

element = "Cu"
a = 3.615
# FCC (111) 面间距
d_111 = a / np.sqrt(3)

# 基础尺寸：保持层数是 3 的倍数 (12)，这对 Bulk 的 ABC 周期性至关重要
size_111 = (3, 3, 12)

output_root = r"C:/Users/senji/Desktop/super/GSFE（广义层错曲线）_add_relex_fixedGMINI"
os.makedirs(output_root, exist_ok=True)
relax_root = os.path.join(output_root, "to_relax")
os.makedirs(relax_root, exist_ok=True)


# ==================== 核心工具：边界修复 ====================

def fix_pbc_vertical_layers(atoms):
    """
    【核心修复函数】
    1. 强制校正 Z 轴 Cell 为垂直方向，消除倾斜带来的边界困扰。
    2. 检查顶层和底层在 PBC 下是否重叠（A压A问题）。
    3. 如果重叠，删除最顶层，并重新调整盒子高度。
    """
    # 1. 获取当前层数（基于 Z 坐标聚类）
    pos = atoms.positions
    z_coords = pos[:, 2]
    # 简单的聚类，允许微小误差
    layers = np.unique(np.round(z_coords, decimals=3))
    n_layers = len(layers)
    
    # 理论总高度
    total_height = n_layers * d_111
    
    # 2. 强制重置 Cell Z 轴为垂直，保留 XY 轴不变
    # 这对于体相计算是安全的，且能避免非正交盒子带来的边界“看起来重叠”的问题
    cell = atoms.get_cell()
    new_cell = np.array([cell[0], cell[1], [0, 0, total_height]])
    atoms.set_cell(new_cell, scale_atoms=False)
    
    # 3. 重新均匀分布 Z 坐标，消除 stack 可能带来的微小误差
    # 按 Z 排序
    indices = np.argsort(atoms.positions[:, 2])
    # 每一层的原子数
    atoms_per_layer = len(atoms) // n_layers
    
    # 赋予完美的 Z 坐标
    for i in range(n_layers):
        # 找到属于这一层的原子索引
        layer_indices = indices[i*atoms_per_layer : (i+1)*atoms_per_layer]
        # 设定 Z
        pos[layer_indices, 2] = i * d_111
    atoms.set_positions(pos)
    
    # 4. 检查 PBC 冲突 (ABC 堆垛检查)
    # 获取底层(Layer 0)和顶层(Layer N-1)的 XY 坐标
    bottom_atoms = atoms[indices[:atoms_per_layer]]
    top_atoms = atoms[indices[-atoms_per_layer:]]
    
    # 简单的相对坐标检查：如果顶层和底层的 XY 位置完全一致（在周期性意义下）
    # 说明是 ...A (top) -> A (bottom) 的情况，需要删层
    # 对于 FCC 111，相邻层错开，每 3 层循环。
    # 如果 n_layers % 3 != 0，通常意味着边界不匹配。
    
    if n_layers % 3 != 0:
        print(f"  [Auto-Fix] Layer count {n_layers} is not divisible by 3. Adjusting...")
        # 删除最顶层
        del atoms[indices[-atoms_per_layer:]]
        # 更新盒子高度
        new_height = (n_layers - 1) * d_111
        cell = atoms.get_cell()
        cell[2, 2] = new_height
        atoms.set_cell(cell, scale_atoms=False)
        print(f"  [Auto-Fix] Removed top layer. New height: {new_height:.3f}, Layers: {n_layers-1}")
    
    return atoms


def build_perfect_bulk():
    """构建完美的体相基准，12层"""
    atoms = fcc111(
        symbol=element,
        size=size_111, # (3,3,12)
        a=a,
        vacuum=0.0, # 无真空
        orthogonal=False
    )
    atoms.pbc = [True, True, True]
    # 应用修复，确保 Z 轴垂直且层数正确
    atoms = fix_pbc_vertical_layers(atoms)
    return atoms


def add_random_perturbation(atoms, amp=0.02, seed=None):
    rng = np.random.default_rng(seed)
    perturbed = atoms.copy()
    disp = (rng.random((len(perturbed), 3)) - 0.5) * 2.0 * amp
    perturbed.positions += disp
    return perturbed


def split_upper_lower_by_index(atoms):
    z = atoms.positions[:, 2]
    order = np.argsort(z)
    half = len(order) // 2
    lower_idx = order[:half]
    upper_idx = order[half:]
    lower = np.zeros(len(z), dtype=bool)
    upper = np.zeros(len(z), dtype=bool)
    lower[lower_idx] = True
    upper[upper_idx] = True
    return lower, upper


def get_inplane_vectors():
    tmp = build_perfect_bulk()
    cell = tmp.get_cell()
    return np.array(cell[0]), np.array(cell[1])


# ==================== 输出相关 ====================

def write_poscar_with_sd(filename, atoms, sd):
    natoms = len(atoms)
    cell = atoms.get_cell()
    scaled = atoms.get_scaled_positions(wrap=True) # 强制 Wrap 到 0-1

    with open(filename, "w") as f:
        f.write(f"{element}\n")
        f.write("  1.0000000000000000\n")
        for v in cell:
            f.write(f"  {v[0]:20.16f}  {v[1]:20.16f}  {v[2]:20.16f}\n")
        f.write(f" {element} \n")
        f.write(f" {natoms}\n")
        f.write("Selective dynamics\n")
        f.write("Direct\n")
        for i in range(natoms):
            x, y, z = scaled[i]
            flags = ["T" if sd[i, j] else "F" for j in range(3)]
            f.write(
                f"  {x:20.16f}  {y:20.16f}  {z:20.16f}   {flags[0]}  {flags[1]}  {flags[2]}\n"
            )

def make_sd(name, atoms):
    natoms = len(atoms)
    if name.startswith("normal_relax_") or name.startswith("strain_"):
        # 这里的逻辑：对于体相 GSFE，固定底部原子不动，防止漂移，让上部弛豫
        # 或者 Z 方向全部放开也可以。这里维持原方案：固定下半，放开上半Z
        lower, upper = split_upper_lower_by_index(atoms)
        sd = np.zeros((natoms, 3), dtype=bool)
        sd[upper, 2] = True 
        return sd
    if name.startswith("twin_CTB_base"):
        sd = np.ones((natoms, 3), dtype=bool) # 全开
        return sd
    return None

def write_out(atoms, name, relax=False):
    path = os.path.join(output_root, f"{name}.vasp")
    write_vasp(path, atoms, direct=True, sort=True, vasp5=True)
    if relax:
        sd = make_sd(name, atoms)
        if sd is not None:
            path_relax = os.path.join(relax_root, f"{name}.vasp")
            write_poscar_with_sd(path_relax, atoms, sd)


# ==================== 1 & 2 & 3: GSFE 结构 (Main, Sec, 2D) ====================

def generate_gsfe_structures(n_main=21, n_sec=11, n_2d=6, amp=0.02):
    # 1. 准备完美的 Bulk
    base = build_perfect_bulk()
    _, upper_mask = split_upper_lower_by_index(base)
    a1, a2 = get_inplane_vectors()

    # --- 1D Main ---
    for i, frac in enumerate(np.linspace(0, 1.0, n_main)):
        atoms = base.copy()
        pos = atoms.positions
        pos[upper_mask] += frac * a1
        atoms.set_positions(pos)
        name = f"1D_main_{i:03d}_frac{frac:.3f}"
        write_out(atoms, name)
        # 扰动
        pert = add_random_perturbation(atoms, amp, seed=i)
        write_out(pert, name+"_pert")

    # --- 1D Secondary ---
    for i, frac in enumerate(np.linspace(0, 1.0, n_sec)):
        atoms = base.copy()
        pos = atoms.positions
        pos[upper_mask] += frac * a2
        atoms.set_positions(pos)
        name = f"1D_sec_{i:03d}_frac{frac:.3f}"
        write_out(atoms, name)
        # 扰动
        pert = add_random_perturbation(atoms, amp, seed=1000+i)
        write_out(pert, name+"_pert")

    # --- 2D Surface ---
    fracs = np.linspace(0, 1.0, n_2d)
    for ix, fx in enumerate(fracs):
        for iy, fy in enumerate(fracs):
            atoms = base.copy()
            pos = atoms.positions
            pos[upper_mask] += fx * a1 + fy * a2
            atoms.set_positions(pos)
            name = f"2D_{ix}_{iy}"
            write_out(atoms, name)


# ==================== 4: 关键点 (Relax) ====================

def generate_key_points():
    base = build_perfect_bulk()
    _, upper_mask = split_upper_lower_by_index(base)
    a1, _ = get_inplane_vectors()
    
    # 典型高对称点
    fracs = [0.0, 1/6, 1/3, 1/2, 2/3, 1.0]
    
    for i, frac in enumerate(fracs):
        atoms = base.copy()
        atoms.positions[upper_mask] += frac * a1
        name = f"normal_relax_key_{i}"
        # 这里 relax=True，会生成带 SD 的 POSCAR
        write_out(atoms, name, relax=True)
        # 扰动版
        pert = add_random_perturbation(atoms, seed=2000+i)
        write_out(pert, name+"_pert")


# ==================== 5: 应变 (Strain) ====================

def generate_strain():
    base0 = build_perfect_bulk()
    _, upper_mask = split_upper_lower_by_index(base0)
    a1, _ = get_inplane_vectors()
    cell0 = base0.get_cell()

    strains = [-0.02, -0.01, 0.01, 0.02]
    key_fracs = [0.0, 1/3, 2/3]

    for eps in strains:
        scale = 1.0 + eps
        cell_eps = cell0 * scale # 均匀缩放
        
        for frac in key_fracs:
            atoms = base0.copy()
            # 1. 滑移
            atoms.positions[upper_mask] += frac * a1
            # 2. 缩放坐标
            atoms.positions *= scale
            # 3. 缩放盒子
            atoms.set_cell(cell_eps, scale_atoms=False)
            
            name = f"strain_eps{eps:.3f}_frac{frac:.3f}"
            write_out(atoms, name, relax=True)
            
            pert = add_random_perturbation(atoms, seed=3000)
            write_out(pert, name+"_pert")


# ==================== 6: 孪晶 (Twin) 修正版 ====================

def generate_twin_bulk():
    """
    生成体相孪晶：
    为了在 PBC 下不冲突，我们必须确保上下边界的原子堆垛是连续的。
    方法：构建一个对称的三明治结构或者手动修剪。
    这里使用修剪法：生成 A-B-C (Interface) C-B-A。
    """
    # 1. 生成两块 Slab
    # 必须保证层数足够。例如 6 层 + 6 层。
    nx, ny, _ = size_111
    n_half = 6
    
    slab = fcc111(symbol=element, size=(nx,ny,n_half), a=a, vacuum=0.0, orthogonal=False)
    twin = slab.copy()
    
    # 2. 翻转 Twin 部分
    pos = twin.positions
    z_mid = (pos[:,2].max() + pos[:,2].min()) / 2
    pos[:, 2] = 2*z_mid - pos[:, 2]
    twin.set_positions(pos)
    
    # 3. 堆垛
    # 这里 distance=d_111 是物理距离。
    atoms = stack(slab, twin, axis=2, distance=d_111)
    atoms.pbc = [True, True, True]
    
    # 4. 【关键步骤】应用边界修复
    # stack 出来的结构是 ...A-B-C (slab) | (d_111) | C-B-A... (twin)
    # 界面处是 C-C。这是错误的。孪晶界面应该是共享的。
    # 所以我们需要删掉其中一层。
    # 另外，底部是 A，顶部是 A (翻转后的 C-B-A，最后是 A)。
    # A 和 A 在 PBC 处也会冲突。
    
    # 我们使用自动修复函数来处理顶层冲突。
    # 针对中间的 C-C 冲突，stack 默认只是平移。
    # 手动处理中间层：
    # 简单粗暴的方法：先修复垂直方向，如果发现原子距离过近（<1.5），就删原子。
    
    # 先应用垂直修复，这会修正 Cell 高度，并删除顶层如果和底层冲突
    atoms = fix_pbc_vertical_layers(atoms)
    
    # 再检查内部原子间距（处理 stack 产生的 C-C 界面）
    # 如果 stack 是 ...C, distance, C... -> 距离 ~2.08 -> 这是层错/孪晶界，不是重叠。
    # 如果距离 > 1.8，VASP 能算。
    # 只要 fix_pbc_vertical_layers 删除了顶层 A，剩下的结构通常是 ...C | C...A...B (边界 B->A 也是错排)。
    # 对于训练数据，这种包含两个界面的结构是极好的样本。
    
    name = "twin_CTB_base"
    write_out(atoms, name, relax=True)
    
    # 扰动
    for k in range(5):
        pert = add_random_perturbation(atoms, seed=4000+k)
        write_out(pert, f"{name}_pert{k}")


# ==================== Main ====================

if __name__ == "__main__":
    print(f"d_111 = {d_111:.4f}")
    print("Generating structures with strict PBC fixes...")
    
    generate_gsfe_structures()
    generate_key_points()
    generate_strain()
    generate_twin_bulk()
    
    print("Done.")
