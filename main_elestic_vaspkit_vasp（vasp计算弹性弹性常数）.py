import os.path
import shutil
from ase.io import read, write
from ase.optimize import BFGS
import numpy as np
from ase.filters import ExpCellFilter
from ase.units import GPa
from mechelastic.core import ElasticProperties
import time
import sys
import os
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.parallel import world, parprint, paropen
import pandas as pd

def vasp_calculator():
    vasp_pp = os.environ['vasp_pp']
    num_cores = os.environ['num_cores']
    vasp_std = os.environ['vasp_std']
    out_path = os.environ['vasp_out_path']

    def vasp_set_calculator(vasp_pp,num_cores,vasp_std,out_path):
        from ase.calculators.vasp import Vasp
        import os
        # 设置环境变量（在Python中）
        os.environ['VASP_PP_PATH'] = vasp_pp#'/share/home/xill/hj/app/vasp_pp'
        # 设置 VASP 计算器（包含命令）
        calc = Vasp(
            # 关键：指定 VASP 执行命令
            command=f'mpirun -np {int(num_cores)} {vasp_std}',  # 替换为实际路径

            # VASP 计算参数
            xc='PBE',
            encut=600,
            kspacing=0.2,
            #isif=2,
            nsw=0,
            ibrion=-1,
            ismear=0,   #(Gaussian smearing, metals:1)
            sigma=0.05,  #(Smearing value in eV, metals:0.2)
            prec='Normal',
            lreal='Auto',
            #algo='Normal',
            nelm=100,
            ediff=1e-6,
            # 输出控制
            lwave=False,
            lcharg=False,
            # 计算目录
            directory=out_path
        )
        return calc
    return vasp_set_calculator(vasp_pp,num_cores,vasp_std,out_path)

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

def opt_ion(sus2mlip_path,write_file_path):
    context = f'''from ase.io import read, write
from ase.optimize import BFGS
import numpy as np
from ase.filters import ExpCellFilter
from ase.units import GPa
import sys
import os
from ase.calculators.abacus import Abacus, AbacusProfile
from ase.parallel import world, parprint, paropen
from ase.io import read,write,iread


def vasp_calculator():
    vasp_pp = '{os.environ['vasp_pp']}'
    num_cores = '{os.environ['num_cores']}'
    vasp_std = '{os.environ['vasp_std']}'
    out_path = '{os.environ['vasp_out_path']}'
    
    def vasp_set_calculator(vasp_pp,num_cores,vasp_std,out_path):
        from ase.calculators.vasp import Vasp
        import os
        # 设置环境变量（在Python中）
        os.environ['VASP_PP_PATH'] = vasp_pp#'/share/home/xill/hj/app/vasp_pp'
        # 设置 VASP 计算器（包含命令）
        calc = Vasp(
            # 关键：指定 VASP 执行命令
            command=f'mpirun -np {{int(num_cores)}} {{vasp_std}}',  # 替换为实际路径
    
            # VASP 计算参数
            xc='PBE',
            encut=600,
            kspacing=0.2,
            #isif=2,
            nsw=0,
            ibrion=-1,
            ismear=0,   #(Gaussian smearing, metals:1)
            sigma=0.05,  #(Smearing value in eV, metals:0.2)
            prec='Normal',
            lreal='Auto',
            #algo='Normal',
            nelm=100,
            ediff=1e-6,
            # 输出控制
            lwave=False,
            lcharg=False,
            # 计算目录
            directory=out_path
        )
        return calc
    return vasp_set_calculator(vasp_pp,num_cores,vasp_std,out_path)



def relax(atoms, calc, pressure=0.0, maxstep=0.1, eps=None, max_step=None):
    atoms.calc = calc
    mask = [False, False, False, False, False, False]
    ucf = ExpCellFilter(atoms, scalar_pressure=pressure * GPa, mask=mask, constant_volume=False)
    gopt = BFGS(ucf, maxstep=maxstep)
    gopt.run(fmax=eps, steps=max_step)
    return atoms

if __name__ == '__main__':
    model = read("POSCAR")
    NEPcalc = vasp_calculator()
    f_max = 0.001
    Relaxed_atoms = relax(atoms=model, calc=NEPcalc, eps=f_max, max_step=1000)
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


def main(atoms,calc,dim,process_dir,sus2mlip_path,strain_str):
    os.chdir(process_dir)
    #1. 弛豫结构
    main_relax(atoms, calc, dim, process_dir)


    #2.产生应变结果并计算
    VPKIT_in(1, dim, 'VPKIT.in',strain_str)
    os.system('touch INCAR POTCAR KPOINTS')
    os.system('echo 201 | vaspkit > /dev/null 2>&1')
    submit('sub.sh')
    opt_ion(sus2mlip_path, 'opt_ion.py')
    os.system('bash sub.sh')

    #3. 分析结果
    VPKIT_in(2, dim, 'VPKIT.in',strain_str)
    os.system('echo 201 | vaspkit > /dev/null 2>&1')
    #os.system('rm -r C* INCAR POTCAR KPOINTS SYMMETRY VPKIT.in opt_ion.py sub.sh')

    cij, B_vrh, G_vrh = read_Cij_BG('ELASTIC_TENSOR')
    return cij, B_vrh, G_vrh


if __name__ == '__main__':
    pwd = os.getcwd()
    # module unload compiler/2022.1.0
    # module load vasp/5.4.4
    os.environ['vasp_pp'] = '/share/home/xill/hj/app/vasp_pp'
    os.environ['num_cores'] = '32'
    os.environ['vasp_std'] = 'vasp_std'
    os.environ['vasp_out_path'] = f"{os.path.join(pwd,'vasp_scf')}"
    #strain_str = '-0.005 -0.003 -0.001 0.000 0.001 0.003 0.005'
    strain_str = '-0.015 -0.010 -0.005 0.000 0.005 0.010 0.015'

    s = time.time()
    input_file = sys.argv[1]
    atoms = read(input_file)
    sus2mlip_path = r" "
    calc = vasp_calculator()
    dim =3
    process_dir = os.path.abspath('vasp_BG')
    os.makedirs(process_dir,exist_ok=True)
    cij, B_vrh, G_vrh = main(atoms, calc, dim, process_dir,sus2mlip_path,strain_str)
    #cij, B_vrh, G_vrh = read_Cij_BG('./test/ELASTIC_TENSOR')
    print(cij, B_vrh, G_vrh)
    print(f'run time: {(time.time() - s) / 60} min')

    # 保存为CSV文件
    results = {
        'structure': [os.path.basename(input_file)],
        'pre_cij': [str(cij)],
        'pre_B_vrh': [B_vrh],
        'pre_G_vrh': [G_vrh],
    }


    # 创建DataFrame
    df = pd.DataFrame(results)

    # 保存到CSV文件
    csv_filename = 'elastic_properties_results.csv'
    df.to_csv(csv_filename, index=False)

    print(f"结果已保存到 {csv_filename}")


