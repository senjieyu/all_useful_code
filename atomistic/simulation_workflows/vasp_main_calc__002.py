from ase.io import iread, write
import os
import shutil
#from .gen_calc_file import poscar2STRU,INPUT
from tqdm import trange
import yaml

'''获取目录的最大深度'''
def get_directory_depth(directory):
    max_depth = 0
    for root, _, _ in os.walk(directory):
        depth = root[len(directory):].count(os.sep)  # 计算当前目录相对于基础目录的层级
        max_depth = max(max_depth, depth)
    return max_depth

'''获取npy和tpye的路径'''
def search_directories_in_directory(directory,depth):
    npy_path = []
    for root, dirs, files in os.walk(directory):
        temp_depth = root[len(directory):].count(os.sep)
        if temp_depth == depth:
            npy_path.append(root)
    return npy_path

def check_and_modify_calc_dir(pwd,xyz_list,calc_dir_num):
    parameter_yaml = os.path.join(pwd,'parameter.yaml')
    end_yaml = os.path.join(pwd,'end.yaml')
    xyz_num = len(xyz_list)

    if xyz_num < calc_dir_num:
        new_calc_dir_num = xyz_num
        with open(parameter_yaml, 'r') as file:
            data_1 = yaml.safe_load(file)
        dft = data_1['dft']
        dft['calc_dir_num'] = new_calc_dir_num
        with open(parameter_yaml, 'w') as file:
            yaml.safe_dump(data_1, file, default_flow_style=False)

        with open(end_yaml, 'r') as file:
            data_2 = yaml.safe_load(file)
        dft = data_2['dft']
        dft['calc_dir_num'] = new_calc_dir_num
        with open(end_yaml, 'w') as file:
            yaml.safe_dump(data_2, file, default_flow_style=False)
        return new_calc_dir_num
    else:
        return calc_dir_num

def mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

####数据集划分#######
def mkdir_Placefile(dir,num_posacr,calc_dir_num,xyz):
    quotient = num_posacr // calc_dir_num
    remainder = num_posacr % calc_dir_num
    temp = [quotient] * calc_dir_num
    divide = [a + 1 if index < remainder else a for index, a in enumerate(temp)]
    range_list = [[sum(divide[:a + 1]) - divide[a] + 1, sum(divide[:a + 1])] for a in range(len(divide))]
    for i in trange(1,calc_dir_num+1):
        sub_dir = os.path.join(dir, 'dir_' + str(i))
        mkdir(sub_dir)
        index_poscar = [a for a in range(range_list[i - 1][0], range_list[i - 1][1] + 1)]

        for ii in index_poscar:
            sub_sub_dir = os.path.join(sub_dir, str(ii))
            mkdir(sub_sub_dir)
            b = xyz[ii - 1]
            sorted_indices = sorted(range(len(b)), key=lambda i: b[i].symbol)
            sorted_atoms = b.__class__()
            sorted_atoms.cell = b.cell
            sorted_atoms.pbc = b.pbc
            for i in sorted_indices:
                sorted_atoms.append(b[i])

            write(os.path.join(sub_sub_dir,'POSCAR'), sorted_atoms, format='vasp')
            os.chdir(sub_sub_dir)
            #poscar2STRU(dir)
            #INPUT(dir,sub_sub_dir)
            #os.remove('ase_sort.dat')

def main_calc(atom_list,calc_dir_num,server,scf_cores,scf_queue,scf_ptile):
    #file = [f for f in os.listdir() if f.endswith('.xyz')][0]
    # calc_dir_num=13 #计算文件数目
    #num = 1  #1:scf 2:cell-relax
    #server = 2 #1:qiming 2:wulixi

    pwd = os.getcwd()
    dir = os.path.join(pwd,'filter')
    mkdir(dir)

    num_posacr = len(atom_list)

    mkdir_Placefile(dir,num_posacr,calc_dir_num,atom_list)

    os.chdir(pwd)
    jobs_script = os.path.join(pwd,'jobs_script')
    mkdir(jobs_script)

    quotient = num_posacr // calc_dir_num
    remainder = num_posacr % calc_dir_num
    temp = [quotient] * calc_dir_num
    divide = [a + 1 if index < remainder else a for index, a in enumerate(temp)]
    range_list = [[sum(divide[:a + 1]) - divide[a] + 1, sum(divide[:a + 1])] for a in range(len(divide))]

    ####创建任务提交脚本#####
    for i in range(1,calc_dir_num+1):
        range1,range2 = range_list[i - 1][0], range_list[i - 1][1]
        #print(range1,range2)
        content_1 = ''
        if server == 'qiming':
            content_1 = f'''#!/bin/bash
#BSUB -J a_job_{str(i)}
#BSUB -q {scf_queue}
#BSUB -n {scf_cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={scf_ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std"
'''
        elif server == 'taiyi':
            content_1 = f'''#!/bin/bash
#BSUB -J a_job_{str(i)}
#BSUB -q {scf_queue}
#BSUB -n {scf_cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={scf_ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
module load vasp/5.4.4
#-------------intelmpi+ifort------------------------------------------
source /share/intel/2018u4/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
source /share/intel/2018u4/impi/2018.4.274/intel64/bin/mpivars.sh
COMMAND_std="mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std"
'''
        elif server == 'physics':
            content_1 = f'''#!/bin/bash
#SBATCH --job-name  a_job_{str(i)}
#SBATCH --partition {scf_queue}
##SBATCH --nodelist c0[01-40]
#SBATCH --ntasks {scf_cores}  #number of core
##SBATCH --qos=840cpu
#####

module unload compiler/2022.1.0
module load vasp/5.4.4
COMMAND_std="mpirun vasp_std"

#start_time=$(date +%s.%N)
#$COMMAND_std 
#end_time=$(date +%s.%N)
#runtime=$(echo "$end_time - $start_time" | bc)
'''
        content_2 = f'''dir_1="filter"
dir_2="{'dir_' + str(i)}"
cd ../$dir_1/$dir_2

path=$(pwd)

for item in {{{range1}..{range2}}}; do
    cd $path/$item
    start_time=$(date +%s.%N)
    touch __start__
    $COMMAND_std > logout 2>&1
    touch __ok__
    
    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc)
    cd $path
    cd ..
    echo "job_{str(i)}_$item total_runtime:$runtime s" >> time.txt
done'''
        with open(os.path.join(jobs_script,f'bsub_{i}.lsf'), 'w') as file:
            file.write(content_1+content_2)
        os.chdir(jobs_script)
        #os.system(f'dos2unix bsub_{i}.lsf')

    os.chdir(pwd)
    start_calc = ''

    if server =='physics':
        start_calc = f'''import os
    
start = 1
end = {calc_dir_num}

pwd = os.getcwd()
jobs_script = os.path.join(pwd,'jobs_script')
os.chdir(jobs_script)
for i in range(start,end+1):
    #print(f'python bsub_{{i}}.lsf')
    if not os.path.exists(f'bsub_{{i}}.lsf'):
        print('errors')
    os.system(f'sbatch bsub_{{i}}.lsf')'''

    elif server =='qiming' or server =='taiyi':
        start_calc = f'''import os
start = 1
end = {calc_dir_num}

pwd = os.getcwd()
jobs_script = os.path.join(pwd,'jobs_script')
os.chdir(jobs_script)
for i in range(start,end+1):
    #print(f'python bsub_{{i}}.lsf')
    if not os.path.exists(f'bsub_{{i}}.lsf'):
        print('errors')
    os.system(f'bsub<bsub_{{i}}.lsf')'''

    with open('start_calc.py', 'w') as f:
        f.write(start_calc)
    return num_posacr

def gen_cal_file(INCAR,dirs):
    depth = get_directory_depth(dirs)
    calc_list = search_directories_in_directory(dirs, depth)
    for a in calc_list:
        # 删除目标目录中可能已存在的INCAR（无论是文件还是链接）
        incar_path = os.path.join(a, 'INCAR')
        if os.path.exists(incar_path):
            os.remove(incar_path)  # 删除现有文件/链接
        # 创建符号链接（源INCAR -> 目标目录/INCAR）
        os.symlink(os.path.abspath(INCAR), incar_path)
        #shutil.copy(INCAR,a)
        os.chdir(a)
        os.system(f'(echo 103)|vaspkit')

if __name__ == '__main__':
    pwd = os.getcwd()
    INCAR = os.path.join(pwd,'INCAR')
    dirs = os.path.join(pwd,'filter')

    t = [os.path.join(pwd,a) for a in os.listdir(pwd) if a.endswith('.xyz')][0]
    atom_list = list(iread(t))
    #atom_list = list(iread('vac.xyz'))
    xyz_num = len(atom_list)
    calc_dir_num = 7
    # server = 'taiyi'
    # scf_cores = 40
    # scf_queue = 'short'
    # scf_ptile = 40

    server = 'qiming'
    scf_cores = 40
    scf_queue = '33'
    scf_ptile = 40


    if xyz_num < calc_dir_num:
        new_calc_dir_num = xyz_num
    else:
        new_calc_dir_num = calc_dir_num

    main_calc(atom_list, new_calc_dir_num, server, scf_cores, scf_queue, scf_ptile)
    gen_cal_file(INCAR, dirs)