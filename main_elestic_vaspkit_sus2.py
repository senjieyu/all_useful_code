import glob
import os.path
import shutil

from ase.io import read, write
#from pynep.calculate import NEP
from ase.optimize import BFGS
import numpy as np
from ase.filters import ExpCellFilter
from ase.units import GPa
from mechelastic.core import ElasticProperties
from ase.data import atomic_numbers
from pymlip.core import MTPCalactor,PyConfiguration
from ase.calculators.calculator import Calculator
from typing import Optional, Union, List
from multiprocessing import Pool
from ase import Atoms
import time

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


def relax(atoms, calc, pressure=0.0, maxstep=0.1, eps=None, max_step=None, dim=None):
    atoms.calc = calc
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
    Relaxed_atoms = relax(atoms=atoms, calc=calc, eps=f_max, max_step=1000, dim=dim)

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
from ase.filters import ExpCellFilter
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
    atoms.calc = calc
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

