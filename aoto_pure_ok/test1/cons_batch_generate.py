# -*- coding: utf-8 -*-
# cons_batch_generate.py
import numpy as np
from ase.build import stack
import model_construct as mc


# ... (前五个函数保持不变: generate_1d_main 到 generate_strain) ...

def generate_1d_main(base_atoms, vec, cfg, dirs):
    """任务1: 沿主矢量生成 1D 路径"""
    print(f"  > Executing Task: 1D Main Path")
    mask = mc.AtomUtils.split_upper_lower_mask(base_atoms)
    count = 0
    for i, frac in enumerate(np.linspace(0, 1.0, cfg['n_points'])):
        atoms = base_atoms.copy()
        atoms.positions[mask] += frac * vec
        name = f"1D_main_{i:03d}_frac{frac:.4f}"
        mc.AtomUtils.write_vasp_file(atoms, name, dirs['out'], dirs['relax'], cfg['element'])
        count += 1
        for k in range(cfg['n_perturb']):
            p = mc.AtomUtils.add_perturbation(atoms, cfg['amp'], seed=i * 100 + k)
            mc.AtomUtils.write_vasp_file(p, f"{name}_pert{k:02d}", dirs['out'], dirs['relax'], cfg['element'])
            count += 1
    return count


def generate_1d_sec(base_atoms, vec, cfg, dirs):
    """任务2: 沿次矢量生成 1D 路径"""
    print(f"  > Executing Task: 1D Secondary Path")
    mask = mc.AtomUtils.split_upper_lower_mask(base_atoms)
    count = 0
    for i, frac in enumerate(np.linspace(0, 1.0, cfg['n_points'])):
        atoms = base_atoms.copy()
        atoms.positions[mask] += frac * vec
        name = f"1D_sec_{i:03d}_frac{frac:.4f}"
        mc.AtomUtils.write_vasp_file(atoms, name, dirs['out'], dirs['relax'], cfg['element'])
        count += 1
        for k in range(cfg['n_perturb']):
            p = mc.AtomUtils.add_perturbation(atoms, cfg['amp'], seed=1000 + i * 100 + k)
            mc.AtomUtils.write_vasp_file(p, f"{name}_pert{k:02d}", dirs['out'], dirs['relax'], cfg['element'])
            count += 1
    return count


def generate_2d_surface(base_atoms, v1, v2, cfg, dirs):
    """任务3: 2D Gamma Surface"""
    print(f"  > Executing Task: 2D Surface")
    mask = mc.AtomUtils.split_upper_lower_mask(base_atoms)
    fracs = np.linspace(0, 1.0, cfg['grid_n'])
    count = 0
    for ix, fx in enumerate(fracs):
        for iy, fy in enumerate(fracs):
            atoms = base_atoms.copy()
            atoms.positions[mask] += fx * v1 + fy * v2
            name = f"2D_{ix:02d}_{iy:02d}"
            mc.AtomUtils.write_vasp_file(atoms, name, dirs['out'], dirs['relax'], cfg['element'])
            count += 1
            for k in range(cfg['n_perturb']):
                p = mc.AtomUtils.add_perturbation(atoms, cfg['amp'], seed=2000 + ix * 100 + iy + k)
                mc.AtomUtils.write_vasp_file(p, f"{name}_pert{k:02d}", dirs['out'], dirs['relax'], cfg['element'])
                count += 1
    return count


def generate_key_points(base_atoms, vec, cfg, dirs):
    """任务4: 关键点"""
    print(f"  > Executing Task: Key Points")
    mask = mc.AtomUtils.split_upper_lower_mask(base_atoms)
    count = 0
    for i, frac in enumerate(cfg['fractions']):
        atoms = base_atoms.copy()
        atoms.positions[mask] += frac * vec
        name = f"key_point_{i:02d}_frac{frac:.3f}"
        mc.AtomUtils.write_vasp_file(atoms, name, dirs['out'], dirs['relax'], cfg['element'], need_sd=True)
        count += 1
        for k in range(cfg['n_perturb']):
            p = mc.AtomUtils.add_perturbation(atoms, cfg['amp'], seed=3000 + i * 100 + k)
            mc.AtomUtils.write_vasp_file(p, f"{name}_pert{k:02d}", dirs['out'], dirs['relax'], cfg['element'])
            count += 1
    return count


def generate_strain(base_atoms, vec, cfg, dirs):
    """任务5: 应变结构"""
    print(f"  > Executing Task: Strain Structures")
    mask = mc.AtomUtils.split_upper_lower_mask(base_atoms)
    cell0 = base_atoms.get_cell()
    count = 0
    for eps in cfg['strains']:
        scale = 1.0 + eps
        for frac in cfg['key_fracs']:
            atoms = base_atoms.copy()
            atoms.positions[mask] += frac * vec
            atoms.set_cell(cell0 * scale, scale_atoms=True)
            name = f"strain_eps{eps:.3f}_frac{frac:.3f}"
            mc.AtomUtils.write_vasp_file(atoms, name, dirs['out'], dirs['relax'], cfg['element'], need_sd=True)
            count += 1
            for k in range(cfg['n_perturb']):
                p = mc.AtomUtils.add_perturbation(atoms, cfg['amp'], seed=4000 + int((eps + 0.5) * 1000) + k)
                mc.AtomUtils.write_vasp_file(p, f"{name}_pert{k:02d}", dirs['out'], dirs['relax'], cfg['element'])
                count += 1
    return count


def generate_twin(base_atoms, cfg, dirs):
    """
    任务6: 孪晶结构
    【修正】调用 StructureFactory 的 build_twin_slab 方法
    """
    print(f"  > Executing Task: Twin Structure")

    # 获取结构类型 (需要 main.py 注入到 cfg 中)
    struct_type = cfg.get('struct_type', 'FCC')

    # 使用科学的构建方法 (旋转或镜像)
    atoms = mc.StructureFactory.build_twin_slab(base_atoms, struct_type)

    count = 0
    name = "twin_boundary_model"
    # 写入文件，开启 SD (上下固定)
    mc.AtomUtils.write_vasp_file(atoms, name, dirs['out'], dirs['relax'], cfg['element'], need_sd=True)
    count += 1

    for k in range(cfg['n_perturb']):
        p = mc.AtomUtils.add_perturbation(atoms, cfg['amp'], seed=5000 + k)
        mc.AtomUtils.write_vasp_file(p, f"{name}_pert{k:02d}", dirs['out'], dirs['relax'], cfg['element'])
        count += 1
    return count
