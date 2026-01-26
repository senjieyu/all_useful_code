import glob
import os.path
import shutil

from ase.io import read, write
#from pynep.calculate import NEP
from ase.optimize import BFGS
import numpy as np
from ase.constraints import ExpCellFilter
from ase.units import GPa
from mechelastic.core import ElasticProperties
from ase.data import atomic_numbers
from pymlip.core import MTPCalactor,PyConfiguration
from ase.calculators.calculator import Calculator
from typing import Optional, Union, List
from multiprocessing import Pool
from ase import Atoms
import time
# =============================================================================
# 使用说明（把这段复制到脚本开头即可）
# =============================================================================
# 这个脚本做什么：
#   1) 批量读取结构文件：stru/*.vasp
#   2) 用 MTP 机器学习势（.mtp）+ ASE 做结构弛豫（可含晶胞自由度）
#   3) 调用 vaspkit 的 201 功能生成应变结构（strain_*），并对每个应变结构做离子弛豫
#   4) vaspkit 后处理输出 ELASTIC_TENSOR（6×6 Cij）
#   5) 用 mechelastic 从 Cij 计算 B_vrh（体积模量）和 G_vrh（剪切模量）
#   6) 输出 results.csv 和 pre_BG.pkl
#
# -----------------------------------------------------------------------------
# 一、运行前的环境要求
# -----------------------------------------------------------------------------
# 1) Python 包：
#   - ase
#   - numpy
#   - pandas（脚本里会用来写 results.csv）
#   - pymlip（提供 MTPCalactor / PyConfiguration，用于加载 .mtp 势）
#   - mechelastic（ElasticProperties，用于从 Cij 求 B/G）
#
# 2) 外部程序：
#   - vaspkit 必须能在命令行直接运行（脚本使用：echo 201 | vaspkit）
#   - 脚本会在每个应变目录里运行：python opt_ion.py
#
# 3) 文件/目录准备：
#   - 必须存在目录：stru/
#   - 结构文件放在：stru/*.vasp
#     （每个文件名会作为结构标签，例如 stru/xxx.vasp -> 标签 xxx）
#   - 需要一个 MTP 势文件：current_0.mtp（或你自己的 .mtp）
#
# -----------------------------------------------------------------------------
# 二、你需要改哪些参数（脚本末尾 __main__ 区域）
# -----------------------------------------------------------------------------
# 下面这些变量是你最常改的：
#
# (1) 结构路径（默认不改）：
#     structures = glob.glob(os.path.join('stru','*.vasp'))
#
# (2) 并行进程数：
#     num_processes = 1
#     - 想并行就改成 CPU 核数，例如 4/8/16
#
# (3) 体系维度：
#     dim = 3
#     - 3：块体（bulk）
#     - 2：二维材料/薄膜（slab），脚本里会在 relax 阶段设 constant_volume=True
#
# (4) 输出目录：
#     output_base_dir = 'sus2_BG'
#     - 每个结构会在 sus2_BG/<结构名>/ 下生成计算目录
#
# (5) 元素列表（用于 unique_numbers 给 MTP）：
#     ele = ['Cu']
#     - 如果是多元体系，例如 Cu-Zn：ele = ['Cu','Zn']
#     - 注意：这些元素必须与 .mtp 势训练时支持的元素一致
#
# (6) MTP 势文件路径（最重要）：
#     sus2mlip_path = r'/work/.../current_0.mtp'
#     - 这里必须是你机器上的真实路径
#     - 脚本会检查该文件是否存在，不存在会直接报错退出
#
# (7) 应变幅度（7 个值，对应 vaspkit.in 的 strain 行）：
#     strain_str = '-0.015 -0.010 -0.005 0.000 0.005 0.010 0.015'
#     - 必须是 7 个数（脚本写死：number of strain = 7）
#     - 典型小应变：±0.005 或 ±0.01；你可以按需要改
#
# -----------------------------------------------------------------------------
# 三、目录与文件输出说明（脚本运行后会生成什么）
# -----------------------------------------------------------------------------
# 顶层输出目录（output_base_dir，例如 sus2_BG）：
#   sus2_BG/
#     results.csv              # 汇总表（每个结构一行）
#     <结构名1>/
#       POSCAR                 # 初始结构弛豫后的结构（如果弛豫成功才写）
#       ELASTIC_TENSOR         # vaspkit 后处理输出（如果流程走完才有）
#       ...（中间会产生 C*/ strain_* 等目录，最后脚本会 rm -r 清理部分文件）
#     <结构名2>/
#       ...
#
# 另外脚本最后还会保存：
#   pre_BG.pkl                 # Python pickle，保存 results 列表（包含 Cij/B/G）
#
# -----------------------------------------------------------------------------
# 四、运行方式（最常用）
# -----------------------------------------------------------------------------
# 1) 进入脚本所在目录（要求下面这个结构）：
#    .
#    ├── your_script.py
#    └── stru/
#        ├── a.vasp
#        ├── b.vasp
#        └── ...
#
# 2) 修改 __main__ 区域的参数（至少改 MTP 势路径 sus2mlip_path）
#
# 3) 运行：
#    python your_script.py
#
# 4) 查看结果：
#    - 汇总：sus2_BG/results.csv
#    - 单个结构：sus2_BG/<结构名>/ELASTIC_TENSOR
#    - pickle：pre_BG.pkl
#
# -----------------------------------------------------------------------------
# 五、这个脚本的计算流程细节（你想知道它到底在每一步干啥）
# -----------------------------------------------------------------------------
# 对每个结构文件（例如 stru/a.vasp），会做：
#
# Step A：结构弛豫（main_relax -> relax）
#   - ASE + BFGS + ExpCellFilter
#   - fmax = 0.001（力收敛阈值）
#   - max_step = 1000
#   - dim=3 时：constant_volume=False（允许晶胞变化）
#   - dim=2 时：constant_volume=True（固定体积）
#   - 成功后写出：sus2_BG/a/POSCAR
#
# Step B：调用 vaspkit 201 生成应变任务（VPKIT.in + vaspkit）
#   - 写 VPKIT.in（signal=1：预处理）
#   - 运行：echo 201 | vaspkit
#   - vaspkit 会创建类似 Cij/strain_* 的目录结构（取决于 vaspkit版本/设置）
#
# Step C：对每个应变结构做离子弛豫（sub.sh + opt_ion.py）
#   - 生成 sub.sh，循环进入每个 strain_* 目录
#   - 在每个 strain_* 目录里：
#       ln -s ${root_path}/opt_ion.py ./
#       python opt_ion.py > OUTCAR
#   - opt_ion.py 内部：
#       * 读 POSCAR
#       * 用同一个 MTP 势算力
#       * BFGS 弛豫离子（cell 不动：mask 全 False）
#       * 输出 CONTCAR
#
# Step D：vaspkit 后处理得到弹性张量（再次 vaspkit 201）
#   - 写 VPKIT.in（signal=2：后处理）
#   - 运行：echo 201 | vaspkit
#   - 得到：ELASTIC_TENSOR
#
# Step E：读取 Cij 并计算 B、G（read_Cij_BG）
#   - 从 ELASTIC_TENSOR 读出 6×6 Cij（GPa）
#   - mechelastic 计算：
#       B_vrh = elastic_props.K_vrh
#       G_vrh = elastic_props.G_vrh
#
# -----------------------------------------------------------------------------
# 六、常见问题与排错建议（非常实用）
# -----------------------------------------------------------------------------
# 1) 找不到 vaspkit / 运行失败：
#   - 终端里先测试：vaspkit 是否能运行
#   - 再测试：echo 201 | vaspkit 是否会卡住/报错
#
# 2) MTP 势文件不匹配元素：
#   - ele 列表必须和势支持的元素一致
#   - 多元素一定要把所有元素都写进 ele
#
# 3) stru/*.vasp 读不进来：
#   - 确认文件确实是 VASP 格式且 ASE 能读
#   - 可以单独测试：python -c "from ase.io import read; read('stru/a.vasp')"
#
# 4) 结果全是 None：
#   - 表示弛豫失败或 vaspkit/后处理失败
#   - 去对应目录看日志：
#       sus2_BG/<结构名>/OUTCAR（opt_ion.py 输出）
#       以及 vaspkit 产生的文件是否存在
#
# 5) 应变幅度太大导致不稳定：
#   - 先用小应变（例如 ±0.005）验证流程通了再加大
#
# -----------------------------------------------------------------------------
# 七、重要提示（脚本里一个容易误判的地方）
# -----------------------------------------------------------------------------
# main_relax 里判断收敛的代码：
#     if (np.max(atom_force <= f_max)):
# 这在逻辑上容易误判（它对每个分量比较，而不是比较力的模长/绝对值）。
# 更稳的写法通常是：
#     if np.max(np.linalg.norm(atom_force, axis=1)) <= f_max:
# 或：
#     if np.max(np.abs(atom_force)) <= f_max:
# 如果你发现明明没收敛也被当成收敛，建议改这里。
# =============================================================================

class MTPCalculator(Calculator):
    """
    MTP calculator based on ase Calculator
    """

    implemented_properties = ["energy",  "forces", "energies" ,"stress"]

    def __init__(self, potential: str = "p.mtp", mtpcalc: MTPCalactor = None, unique_numbers: List = None ,compute_stress: bool = True, stress_weight: float = 1.0,print_EK: bool = True, **kwargs):
        """

        Args:
            potential (str): xxx.mtp
            compute_stress (bool): whether to calculate the stress
            stress_weight (float): the stress weight.
            **kwargs:
        """
        super().__init__(**kwargs)
        self.potential = potential   ##  xxx.mtp
        self.compute_stress = compute_stress
        self.print_EK = print_EK
        self.stress_weight = stress_weight
        self.mtpcalc = MTPCalactor(self.potential)
        self.unique_numbers=unique_numbers

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties: Optional[list] = None,
        system_changes: Optional[list] = None,
        unique_numbers: Optional[list] = None
    ):
        """
        Args:
            atoms (ase.Atoms): ase Atoms object
            properties (list): list of properties to calculate
            system_changes (list): monitor which properties of atoms were
                changed for new calculation. If not, the previous calculation
                results will be loaded.
        Returns:

        """
        properties = properties or ["energy"]
        system_changes = system_changes or all_changes
        super().calculate(atoms=atoms, properties=properties, system_changes=system_changes)

        # graph = self.potential.graph_converter(atoms)
        # graph_list = graph.as_tf().as_list()
        # results = self.potential.get_efs_tensor(graph_list, include_stresses=self.compute_stress)
        cfg=PyConfiguration.from_ase_atoms(atoms,unique_numbers=self.unique_numbers)
       # print(cfg.types)
        V=atoms.cell.volume
        self.mtpcalc.calc(cfg)
        self.results['energy'] = np.array(cfg.energy)
        self.results['forces'] = cfg.force
        self.results['energies'] = np.array(cfg.site_energys)
        #print(f"PE:{np.array(cfg.energy):.4f} (ev)")

        if self.compute_stress:
            self.results['stress'] = -np.array([cfg.stresses[0,0],cfg.stresses[1,1],cfg.stresses[2,2],cfg.stresses[1,2],cfg.stresses[0,2],cfg.stresses[0,1]])* self.stress_weight/V


def relax(atoms, mlip_calc, pressure=0.0, maxstep=0.1, eps=None, max_step=None, dim=None):
    atoms.set_calculator(mlip_calc)
    mask = [True, True, True, True, True, True]
    if (dim == 2):
        ucf = ExpCellFilter(atoms, scalar_pressure=pressure * GPa, mask=mask, constant_volume=True)
    else:
        ucf = ExpCellFilter(atoms, scalar_pressure=pressure * GPa, mask=mask, constant_volume=False)

    gopt = BFGS(ucf, maxstep=maxstep,logfile=None)
    gopt.run(fmax=eps, steps=max_step)
    return atoms

def main_relax(atoms,calc,dim,process_dir):
    # model = read("POSCAR")
    # NEPcalc = NEP('nep.txt')
    f_max = 0.001
    Relaxed_atoms = relax(atoms=atoms, mlip_calc=calc, eps=f_max, max_step=1000, dim=dim)

    label = False
    atom_force = Relaxed_atoms.get_forces()
    #U_atom = Relaxed_atoms.get_potential_energy()
    if (np.max(atom_force <= f_max)):
        # print("  free  energy   TOTEN  =", U_atom, "eV")
        # print("                 Voluntary context switches:")
        label = True
        write(os.path.join(process_dir,'POSCAR'), Relaxed_atoms, vasp5=True, direct=True)
    return label

def VPKIT_in(signal,dim,write_file_path,strain_str):
    context = f"""{signal}      ! 1 for prep-rocessing, 2 for post-processing
{dim}D                                                             ! 2D for slab, 3D for bulk
7                                                              ! number of strain
{strain_str}                  ! magnitude of strain
"""
    with open(write_file_path,'w') as f:
        f.write(context)

def submit(write_file_path):
    context = '''#!/bin/bash
root_path=`pwd`
for cij in `ls -F | grep /$`
do
  cd ${root_path}/$cij
  for s in strain_*
  do
    cd ${root_path}/$cij/$s
    #echo `pwd`
  ln -s ${root_path}/opt_ion.py  ./
  #ln -s ${root_path}/nep.txt ./
  python opt_ion.py > OUTCAR
    # Add here your vasp_submit_job_script
  done
done
'''
    with open(write_file_path,'w') as f:
        f.write(context)

def opt_ion(sus2mlip_path,u_n,write_file_path):
    context = f'''from ase.io import read, write
#from pynep.calculate import NEP
from ase.optimize import BFGS
import numpy as np
from ase.constraints import ExpCellFilter
from ase.units import GPa
from ase.data import atomic_numbers
from pymlip.core import MTPCalactor,PyConfiguration
from ase.calculators.calculator import Calculator
from typing import Optional, Union, List
from ase import Atoms


class MTPCalculator(Calculator):
    """
    MTP calculator based on ase Calculator
    """

    implemented_properties = ["energy",  "forces", "energies" ,"stress"]

    def __init__(self, potential: str = "p.mtp", mtpcalc: MTPCalactor = None, unique_numbers: List = None ,compute_stress: bool = True, stress_weight: float = 1.0,print_EK: bool = True, **kwargs):
        """

        Args:
            potential (str): xxx.mtp
            compute_stress (bool): whether to calculate the stress
            stress_weight (float): the stress weight.
            **kwargs:
        """
        super().__init__(**kwargs)
        self.potential = potential   ##  xxx.mtp
        self.compute_stress = compute_stress
        self.print_EK = print_EK
        self.stress_weight = stress_weight
        self.mtpcalc = MTPCalactor(self.potential)
        self.unique_numbers=unique_numbers

    def calculate(
        self,
        atoms: Optional[Atoms] = None,
        properties: Optional[list] = None,
        system_changes: Optional[list] = None,
        unique_numbers: Optional[list] = None
    ):
        """
        Args:
            atoms (ase.Atoms): ase Atoms object
            properties (list): list of properties to calculate
            system_changes (list): monitor which properties of atoms were
                changed for new calculation. If not, the previous calculation
                results will be loaded.
        Returns:

        """
        properties = properties or ["energy"]
        system_changes = system_changes or all_changes
        super().calculate(atoms=atoms, properties=properties, system_changes=system_changes)

        # graph = self.potential.graph_converter(atoms)
        # graph_list = graph.as_tf().as_list()
        # results = self.potential.get_efs_tensor(graph_list, include_stresses=self.compute_stress)
        cfg=PyConfiguration.from_ase_atoms(atoms,unique_numbers=self.unique_numbers)
       # print(cfg.types)
        V=atoms.cell.volume
        self.mtpcalc.calc(cfg)
        self.results['energy'] = np.array(cfg.energy)
        self.results['forces'] = cfg.force
        self.results['energies'] = np.array(cfg.site_energys)
        #print(f"PE:{{np.array(cfg.energy):.4f}} (ev)")

        if self.compute_stress:
            self.results['stress'] = -np.array([cfg.stresses[0,0],cfg.stresses[1,1],cfg.stresses[2,2],cfg.stresses[1,2],cfg.stresses[0,2],cfg.stresses[0,1]])* self.stress_weight/V

def relax(atoms, calc, pressure=0.0, maxstep=0.1, eps=None, max_step=None):
    atoms.set_calculator(calc)
    mask = [False, False, False, False, False, False]
    ucf = ExpCellFilter(atoms, scalar_pressure=pressure * GPa, mask=mask, constant_volume=False)
    gopt = BFGS(ucf, maxstep=maxstep)
    gopt.run(fmax=eps, steps=max_step)
    return atoms

if __name__ == '__main__':
    model = read("POSCAR")
    calc = MTPCalculator(potential=r"{sus2mlip_path}", unique_numbers={u_n})
    f_max = 0.001
    Relaxed_atoms = relax(atoms=model, calc=calc, eps=f_max, max_step=1000)
    atom_force = Relaxed_atoms.get_forces()

    U_atom = Relaxed_atoms.get_potential_energy()
    if (np.max(atom_force <= f_max)):
        print("  free  energy   TOTEN  =", U_atom, "eV")
        print("                 Voluntary context switches:")

    write('CONTCAR', Relaxed_atoms, vasp5=True, direct=True)
'''
    with open(write_file_path,'w') as f:
        f.write(context)

def read_Cij_BG(file_path):
    """从 ELASTIC_TENSOR 文件中提取 6x6 弹性常数矩阵 (GPa)"""
    cij = []
    with open(file_path, 'r') as f:
        for line in f:
            # 跳过说明行（如第一行）
            if line.strip() and not line.startswith(('Read', '#')):
                # 提取每行的数字
                row = list(map(float, line.split()))
                if len(row) == 6:  # 确保是6个元素的行
                    cij.append(row)
    cij = np.array(cij)
    elastic_props = ElasticProperties(elastic_tensor=cij)
    B_vrh = elastic_props.K_vrh
    G_vrh = elastic_props.G_vrh
    return cij, B_vrh, G_vrh


def main(atoms,calc,dim,process_dir,sus2mlip_path,u_n,strain_str):

    #1. 弛豫结构
    try:
        label = main_relax(atoms, calc, dim, process_dir)
    except:
        label = False

    os.chdir(process_dir)
    if label:
        try:
            #2.产生应变结果并计算
            VPKIT_in(1, dim, 'VPKIT.in',strain_str)
            os.system('touch INCAR POTCAR KPOINTS')
            os.system('echo 201 | vaspkit > /dev/null 2>&1')
            submit('sub.sh')
            opt_ion(sus2mlip_path, u_n, 'opt_ion.py')
            os.system('bash sub.sh')

            #3. 分析结果
            VPKIT_in(2, dim, 'VPKIT.in',strain_str)
            os.system('echo 201 | vaspkit > /dev/null 2>&1')
            os.system('rm -r C* INCAR POTCAR KPOINTS SYMMETRY VPKIT.in opt_ion.py sub.sh')

            cij, B_vrh, G_vrh = read_Cij_BG('ELASTIC_TENSOR')
            return cij, B_vrh, G_vrh
        except:
            return None, None, None
    else:
        return None, None, None

def process_structure(structure_file, dim, output_base_dir, struct_name, sus2mlip_path,u_n,main_path,strain_str):
    """处理单个结构文件的函数"""
    os.chdir(main_path)
    # 创建输出目录
    output_dir = os.path.join(output_base_dir, struct_name)
    os.makedirs(output_dir, exist_ok=True)

    # 读取结构
    atoms = read(structure_file)

    # 设置计算器
    calc = MTPCalculator(potential=sus2mlip_path, unique_numbers=u_n)

    # 运行主计算
    cij, B_vrh, G_vrh = main(atoms, calc, dim, output_dir, sus2mlip_path, u_n,strain_str)

    return {
        'structure': structure_file,
        'pre_cij': cij,
        'pre_B_vrh': B_vrh,
        'pre_G_vrh': G_vrh,
    }

def parallel_process_structures(structures_info, dim, output_base_dir,sus2mlip_path, u_n, main_path,strain_str,num_processes=None):
    args = [(a[0], dim, output_base_dir, a[1], sus2mlip_path, u_n, main_path,strain_str) for a in structures_info]
    # 使用进程池并行处理
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(process_structure, args)

    import pandas as pd
    # 保存为CSV
    df = pd.DataFrame(results)
    csv_path = os.path.join(output_base_dir, "results.csv")
    df.to_csv(csv_path, index=False)
    print(f"The results are saved to : {csv_path}")

    return results


def add_results_to_dataframe(original_df, results):
    """将计算结果添加到原始DataFrame中"""
    # 创建结果DataFrame
    results_df = pd.DataFrame(results)

    # 重命名列以添加pre_前缀
    results_df = results_df.rename(columns={
        'B_vrh': 'pre_k',
        'G_vrh': 'pre_g',
        'cij': 'pre_cij'
    })

    # 合并到原始DataFrame
    merged_df = original_df.merge(results_df, how='left')

    return merged_df

if __name__ == '__main__':
    s = time.time()

    # name = 'POSCAR'
    # atoms = read(name)
    # ele = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ca',
    #  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Rb',
    #  'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
    #  'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    #  'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Ac', 'Th', 'Pa', 'U', 'Np',
    #  'Pu']
    # #sus2mlip_path = r'/work/home/xill/MD_test/current.mtp'
    # sus2mlip_path = r'/home/jinghuang/database/mp_property/current.mtp'
    # dim = 3
    #
    # u_n = sorted([atomic_numbers[i] for i in ele])
    # output_base_dir = 'ttt'
    # struct_name = '1'
    # tt = process_structure(name, dim,output_base_dir, struct_name, u_n)
    # print(tt)
    #
    # print(f'run time: {(time.time()-s)/60} min')


    # import pandas as pd
    # import csv
    #
    # data = pd.read_csv('materials_failed.csv')
    # ids = data['material_id'].values
    #
    # dir = r'/work/home/xill/mp_property/failed/structures'
    #
    # structures_info = []
    # for a in ids:
    #     structures_info.append((os.path.join(dir,f'{a}.vasp'),a))

    structures = glob.glob(os.path.join('stru','*.vasp'))
    structures_info = [(stru, str(os.path.basename(stru)).replace('.vasp','')) for stru in structures]

    num_processes=1
    dim = 3
    output_base_dir = 'sus2_BG'
    # ele = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ca',
    #  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Rb',
    #  'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
    #  'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    #  'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Ac', 'Th', 'Pa', 'U', 'Np',
    #  'Pu']
    ele = ['Cu']
    #sus2mlip_path = r'/work/home/xill/MD_test/current.mtp'
    #sus2mlip_path = r'/work/home/xill/mp_property/current.mtp'
    sus2mlip_path = r'/work/phy-luoyf/WORK/Cu_vasp_self_constrct/cu_sample/current_0.mtp'
    if not os.path.exists(sus2mlip_path):
        raise ValueError(f'{sus2mlip_path} is not exist! Please check!')
    u_n = sorted([atomic_numbers[i] for i in ele])
    main_path = os.getcwd()
    # strain_str = '-0.005 -0.003 -0.001 0.000 0.001 0.003 0.005'
    strain_str = '-0.015 -0.010 -0.005 0.000 0.005 0.010 0.015'

    results = parallel_process_structures(structures_info, dim, output_base_dir, sus2mlip_path, u_n, main_path,strain_str,num_processes=num_processes)

    import pickle

    # 指定保存的文件路径
    file_path = 'pre_BG.pkl'

    # 使用 pickle 将列表数据保存为 pkl 文件
    with open(file_path, 'wb') as f:
        pickle.dump(results, f)

    with open(file_path, 'rb') as f:
        loaded_data = pickle.load(f)
    #print(loaded_data)
    #print(results)
    # # 将结果添加到原始数据
    # final_data = add_results_to_dataframe(data,results)
    #
    # # 保存结果
    # output_csv = 'test_with_predictions.csv'
    # final_data.to_csv(output_csv, index=False)

    print(f'run time: {(time.time() - s) / 60} min')

