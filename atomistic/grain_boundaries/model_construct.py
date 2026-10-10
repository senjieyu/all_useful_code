# -*- coding: utf-8 -*-
# model_construct.py
import os
import numpy as np
import shutil
import subprocess
import tempfile

from ase.build import fcc111, bcc110, hcp0001, bulk, surface, stack
from ase.io.vasp import write_vasp
from ase.io import read as ase_read
from ase.neighborlist import NeighborList


# ============================================================
#  BCC Σ3 {112} twin Atomsk duplicate 固定为你验证过的规模
#  等价于你脚本里的：-duplicate 1 14 1
#  注意：这里 duplicate 的三个方向分别对应：
#    X=[-110], Y=[112], Z=[22-2]
# ============================================================
BCC_TWIN_DUPLICATE = (1, 14, 1)


class AtomUtils:
    """原子操作通用工具集"""

    @staticmethod
    def get_d_spacing(atoms, axis=2):
        """计算沿特定轴的平均层间距"""
        coords = atoms.positions[:, axis]
        layers = np.unique(np.round(coords, 3))
        if len(layers) > 1:
            return np.mean(np.diff(layers))
        return 2.0

    @staticmethod
    def fix_pbc_vertical(atoms, vacuum=0.0):
        """
        修复垂直方向周期性 (Bulk模式下 vacuum=0)
        严格保证盒子高度 = 层数 * 层间距
        """
        z_coords = atoms.positions[:, 2]
        unique_layers = np.unique(np.round(z_coords, 3))
        n_layers = len(unique_layers)

        if n_layers < 2:
            atoms.pbc = [True, True, True]
            return atoms

        # 1) 均匀化 Z 坐标
        d_spacing = np.mean(np.diff(unique_layers))
        indices = np.argsort(z_coords)
        atoms_per_layer = len(atoms) // n_layers

        pos = atoms.positions
        for i in range(n_layers):
            start = i * atoms_per_layer
            end = min((i + 1) * atoms_per_layer, len(atoms))
            if start >= len(atoms):
                break
            layer_idx = indices[start:end]
            pos[layer_idx, 2] = i * d_spacing
        atoms.set_positions(pos)

        # 2) 设置盒子高度
        total_height = n_layers * d_spacing + vacuum
        cell = atoms.get_cell()
        new_cell = np.array([cell[0], cell[1], [0, 0, total_height]])
        atoms.set_cell(new_cell, scale_atoms=False)
        atoms.pbc = [True, True, True]
        return atoms

    @staticmethod
    def split_upper_lower_mask(atoms):
        """用于滑移时的上下分层"""
        z = atoms.positions[:, 2]
        z_mid = (np.max(z) + np.min(z)) / 2.0
        return z > z_mid

    @staticmethod
    def remove_overlaps(atoms, cutoff=0.8):
        """
        强力去重（快速版）：只要距离小于 cutoff (Å) 就删除
        使用 NeighborList，避免 get_all_distances O(N^2) 卡死
        """
        n = len(atoms)
        if n == 0:
            return atoms

        nl = NeighborList([cutoff] * n, bothways=True, self_interaction=False)
        nl.update(atoms)

        to_delete = set()
        for i in range(n):
            if i in to_delete:
                continue
            idxs, offsets = nl.get_neighbors(i)
            for j, off in zip(idxs, offsets):
                if j <= i or j in to_delete:
                    continue
                rij = atoms.positions[j] + np.dot(off, atoms.cell) - atoms.positions[i]
                if np.linalg.norm(rij) < cutoff:
                    to_delete.add(j)

        if to_delete:
            print(f"    [Auto-Fix] Removing {len(to_delete)} overlapped atoms (dist<{cutoff}A)...")
            del atoms[sorted(to_delete)]
        return atoms

    @staticmethod
    def add_perturbation(atoms, amp, seed):
        if amp <= 0:
            return atoms
        rng = np.random.default_rng(seed)
        pert = atoms.copy()
        disp = (rng.random((len(pert), 3)) - 0.5) * 2.0 * amp
        pert.positions += disp
        return pert

    @staticmethod
    def write_vasp_file(atoms, filename, output_dir, relax_dir=None, element="X", need_sd=False):
        """写入文件，包含 SD 设置 (上下固定，中间放开)"""
        full_path = os.path.join(output_dir, f"{filename}.vasp")
        write_vasp(full_path, atoms, direct=True, sort=True, vasp5=True)

        if need_sd and relax_dir:
            sd_path = os.path.join(relax_dir, f"{filename}.vasp")
            sd_flags = np.ones((len(atoms), 3), dtype=bool)

            z_coords = atoms.positions[:, 2]
            z_min = np.min(z_coords)
            z_max = np.max(z_coords)
            tolerance = 0.5

            sd_flags[z_coords < z_min + tolerance, :] = False
            sd_flags[z_coords > z_max - tolerance, :] = False

            with open(sd_path, "w") as f:
                f.write("Generated generic\n1.0\n")
                cell = atoms.get_cell()
                for v in cell:
                    f.write(f"  {v[0]:12.8f}  {v[1]:12.8f}  {v[2]:12.8f}\n")
                f.write(f" {element}\n {len(atoms)}\nSelective dynamics\nDirect\n")
                scaled = atoms.get_scaled_positions(wrap=True)
                for i in range(len(atoms)):
                    s = scaled[i]
                    flags = " ".join(["T" if x else "F" for x in sd_flags[i]])
                    f.write(f"  {s[0]:12.8f}  {s[1]:12.8f}  {s[2]:12.8f} {flags}\n")


class StructureFactory:

    @staticmethod
    def _get_vectors_by_symmetry(atoms, lattice_type, size):
        cell = atoms.get_cell()
        a1_prim = np.array(cell[0]) / size[0]
        a2_prim = np.array(cell[1]) / size[1]
        info = {}
        if lattice_type == "FCC":
            v_main = a1_prim
            info['main_label'] = "[1-10]"
            info['main_desc'] = "Full Slip"
            v_sec = a2_prim
            info['sec_label'] = "[11-2]"
            info['sec_desc'] = "Partial Slip"
        elif lattice_type == "BCC":
            v_main = a1_prim
            info['main_label'] = "[1-11]"
            info['main_desc'] = "Burgers Vector"
            v_sec = a2_prim
            info['sec_label'] = "[001]"
            info['sec_desc'] = "Secondary Slip"
        elif lattice_type == "HCP":
            v_main = a1_prim
            info['main_label'] = "<11-20>"
            info['main_desc'] = "Basal a-type"
            v_sec =  a2_prim
            info['sec_label'] = "<10-10>"
            info['sec_desc'] = "Prism"
        else:
            v_main = a1_prim
            v_sec = a2_prim
        return v_main, v_sec, info

    @staticmethod
    def build_fcc(element, a, size, vacuum=0.0):
        atoms = fcc111(symbol=element, size=size, a=a, vacuum=0.0, orthogonal=False)
        atoms = AtomUtils.fix_pbc_vertical(atoms, vacuum=vacuum)

        # 注入 main 参数，供 twin 读取（不改 main/batch）
        atoms.info['struct_type'] = "FCC"
        atoms.info['element'] = element
        atoms.info['lattice_a'] = float(a)
        atoms.info['lattice_c'] = None
        atoms.info['size'] = tuple(size)
        atoms.info['vacuum'] = float(vacuum)

        v1, v2, info = StructureFactory._get_vectors_by_symmetry(atoms, "FCC", size)
        return atoms, v1, v2, info

    @staticmethod
    def build_bcc(element, a, size, vacuum=0.0):
        atoms = bcc110(symbol=element, size=size, a=a, vacuum=0.0, orthogonal=False)
        atoms = AtomUtils.fix_pbc_vertical(atoms, vacuum=vacuum)

        # 注入 main 参数，供 twin 读取（不改 main/batch）
        atoms.info['struct_type'] = "BCC"
        atoms.info['element'] = element
        atoms.info['lattice_a'] = float(a)
        atoms.info['lattice_c'] = None
        atoms.info['size'] = tuple(size)
        atoms.info['vacuum'] = float(vacuum)

        v1, v2, info = StructureFactory._get_vectors_by_symmetry(atoms, "BCC", size)
        return atoms, v1, v2, info

    @staticmethod
    def build_hcp(element, a, c, size, vacuum=0.0):
        atoms = hcp0001(symbol=element, size=size, a=a, c=c, vacuum=0.0, orthogonal=False)
        atoms = AtomUtils.fix_pbc_vertical(atoms, vacuum=vacuum)

        # 注入 main 参数，供 twin 读取（不改 main/batch）
        atoms.info['struct_type'] = "HCP"
        atoms.info['element'] = element
        atoms.info['lattice_a'] = float(a)
        atoms.info['lattice_c'] = float(c)
        atoms.info['size'] = tuple(size)
        atoms.info['vacuum'] = float(vacuum)

        v1, v2, info = StructureFactory._get_vectors_by_symmetry(atoms, "HCP", size)
        return atoms, v1, v2, info

    # ============================================================
    #  Atomsk: BCC Σ3 {112} twin（mirror + merge，你验证正确）
    # ============================================================
    @staticmethod
    def _atomsk_build_bcc_sigma3_112_twin_fixed(element, lattice_a):
        """
        调用 Atomsk 生成 BCC Σ3 {112} twin
        duplicate 固定使用 (1,14,1)，对应你确认的原子数量
        """
        atomsk = shutil.which("atomsk") or shutil.which("atomsk.exe")
        if atomsk is None:
            raise RuntimeError("找不到 atomsk/atomsk.exe，请把 Atomsk 加入 PATH。")

        nx, ny, nz = BCC_TWIN_DUPLICATE
        print(f"    [Twin] Atomsk duplicate fixed to: ({nx},{ny},{nz})")

        with tempfile.TemporaryDirectory() as td:
            A = os.path.join(td, "Fe_A.xsf")
            B = os.path.join(td, "Fe_B.xsf")
            T = os.path.join(td, "Fe_twin.xsf")
            POS = os.path.join(td, "POSCAR")

            def run(cmd):
                p = subprocess.run(cmd, cwd=td, text=True, capture_output=True)
                if p.returncode != 0:
                    raise RuntimeError(
                        f"[Atomsk ERROR]\nCMD: {' '.join(cmd)}\n"
                        f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
                    )

            # 1) create + orient + duplicate
            run([
                atomsk,
                "--create", "bcc", str(lattice_a), element,
                "orient", "[-110]", "[112]", "[22-2]",
                "-duplicate", str(nx), str(ny), str(nz),
                A
            ])

            # 2) mirror about Y
            run([atomsk, A, "-mirror", "0", "Y", "-wrap", B])

            # 3) merge grains along Y
            run([atomsk, "--merge", "Y", "2", A, B, T])

            # 4) export POSCAR
            run([atomsk, T, "pos"])

            if not os.path.exists(POS):
                raise RuntimeError("Atomsk 未生成 POSCAR，请检查 Atomsk 输出。")

            atoms = ase_read(POS)
            atoms.wrap()
            return atoms

    # =========================================================================
    #  Twin builder
    # =========================================================================
    @staticmethod
    def build_twin_slab(base_atoms, structure_type):
        """
        生成孪晶结构：
        - BCC：用 Atomsk Σ3 {112}（正确，且 duplicate 固定 1 14 1）
        - FCC/HCP：保留你原来的 stacking（不动）
        """
        elem = base_atoms.get_chemical_symbols()[0]

        # ---------------- BCC：Atomsk 正确孪晶 ----------------
        if structure_type == "BCC":
            print("    [Twin] BCC: Using Atomsk Σ3 {112} (mirror + merge, validated)")

            element = base_atoms.info.get("element", elem)
            lattice_a = base_atoms.info.get("lattice_a", None)
            vacuum = base_atoms.info.get("vacuum", 0.0)

            # 兜底：如果没注入成功，才用估算
            if lattice_a is None:
                dists = base_atoms.get_all_distances(mic=True)
                np.fill_diagonal(dists, 100.0)
                min_dist = np.min(dists)
                lattice_a = 2.0 / np.sqrt(3.0) * min_dist

            atoms = StructureFactory._atomsk_build_bcc_sigma3_112_twin_fixed(
                element=str(element),
                lattice_a=float(lattice_a),
            )

            atoms = AtomUtils.fix_pbc_vertical(atoms, vacuum=float(vacuum))
            atoms = AtomUtils.remove_overlaps(atoms, cutoff=1.0)
            atoms = AtomUtils.fix_pbc_vertical(atoms, vacuum=float(vacuum))
            return atoms

        # ---------------- 下面保持你原来的 stacking twin（FCC/HCP） ----------------
        dists = base_atoms.get_all_distances(mic=True)
        np.fill_diagonal(dists, 100.0)
        min_dist = np.min(dists)

        slab_lower = None

        if structure_type == "FCC":
            print("    [Twin] FCC: Stacking (111) Matrix + Rotated Twin")
            a_est = min_dist * np.sqrt(2)
            slab_lower = fcc111(elem, size=(3, 3, 6), a=a_est, vacuum=0.0, orthogonal=False)

        elif structure_type == "HCP":
            print("    [Twin] HCP: Stacking {10-12} Matrix + Rotated Twin")
            a_est = min_dist
            c_est = 1.633 * a_est
            b_obj = bulk(elem, 'hcp', a=a_est, c=c_est)
            slab_lower = surface(b_obj, (1, 0, 2), layers=8, vacuum=0.0)
            slab_lower.center(vacuum=0.0, axis=2)
            slab_lower = slab_lower.repeat((3, 2, 1))
        else:
            return base_atoms

        slab_lower = AtomUtils.fix_pbc_vertical(slab_lower, vacuum=0.0)
        d_spacing = AtomUtils.get_d_spacing(slab_lower)
        print(f"    [Twin] Layer spacing detected: {d_spacing:.3f} A")

        slab_upper = slab_lower.copy()

        # 180°绕 Z 旋转
        pos = slab_upper.positions
        pos[:, 0] = -pos[:, 0]
        pos[:, 1] = -pos[:, 1]
        slab_upper.set_positions(pos)

        atoms = stack(slab_lower, slab_upper, axis=2, distance=d_spacing)

        atoms = AtomUtils.fix_pbc_vertical(atoms, vacuum=0.0)
        atoms = AtomUtils.remove_overlaps(atoms, cutoff=1.0)
        atoms = AtomUtils.fix_pbc_vertical(atoms, vacuum=0.0)

        return atoms
