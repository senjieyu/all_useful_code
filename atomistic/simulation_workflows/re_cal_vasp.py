import pickle
import os
import shutil


def load_from_pkl(filename):
    """
    从 pkl 文件读取数据

    Args:
        filename: 文件名

    Returns:
        读取的数据
    """
    with open(filename, 'rb') as file:  # 'rb' 表示以二进制读取模式打开
        data = pickle.load(file)
    return data

def gen_file(calc_list,INCAR,BSUB):
    for dir_path in calc_list:
        # 创建INCAR符号链接
        incar = os.path.join(dir_path, 'INCAR')
        if os.path.exists(incar) or os.path.islink(incar):
            os.remove(incar)  # 删除已存在的文件或链接
        os.symlink(os.path.abspath(INCAR), incar)

        # 创建BSUB符号链接
        if BSUB is not None:
            bsub = os.path.join(dir_path, 'bsub.lsf')
            if os.path.exists(bsub) or os.path.islink(bsub):
                os.remove(bsub)  # 删除已存在的文件或链接
            os.symlink(os.path.abspath(BSUB), bsub)
        #shutil.copy(INCAR,a)
        os.chdir(dir_path)
        os.system(f'(echo 103)|vaspkit')

def content(server,cores,queue,ptile,dir_name,dir_num,num_list):
    tt = ' '.join(f"{num}" for num in num_list)
    if server == 'physics':
        s = f"""#!/bin/bash
#SBATCH --job-name  a_{dir_name}
#SBATCH --partition {queue}
##SBATCH --nodelist c0[01-40]
#SBATCH --ntasks {cores}  #number of core
##SBATCH --qos=840cpu
#####

module unload compiler/2022.1.0
module load vasp/5.4.4
COMMAND_std="mpirun vasp_std"

#start_time=$(date +%s.%N)
#$COMMAND_std > logout 2>&1
#end_time=$(date +%s.%N)
#runtime=$(echo "$end_time - $start_time" | bc)

users=({tt})
dir_1="filter"
dir_2="{dir_name}"
cd ../$dir_1/$dir_2

path=$(pwd)

for item in "${{users[@]}}"; do
    cd $path/$item
    start_time=$(date +%s.%N)
    touch __start__
    $COMMAND_std > logout 2>&1
    touch __ok__

    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc)
    cd $path
    cd ..
    echo "re_job_{dir_num}_$item total_runtime:$runtime s" >> time.txt
done
    """
    elif server =='qiming':
        s = f"""#!/bin/bash
#BSUB -J  a_{dir_name}
#BSUB -q {queue}
#BSUB -n {cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std"

users=({tt})
dir_1="filter"
dir_2="{dir_name}"
cd ../$dir_1/$dir_2

path=$(pwd)

for item in "${{users[@]}}"; do
    cd $path/$item
    start_time=$(date +%s.%N)
    touch __start__
    $COMMAND_std > logout 2>&1
    touch __ok__

    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc)
    cd $path
    cd ..
    echo "re_job_{dir_num}_$item total_runtime:$runtime s" >> time.txt
done
    """
    else:
        raise ValueError(f'sever is equal to physics or qiming!')
    return s

def sample_bsub(server,cores,queue,ptile,INCAR_path,cal_bool):
    dir_list = load_from_pkl('failed.pkl')
    print(len(dir_list))
    if server == 'qiming':
        s = f"""#!/bin/bash
#BSUB -J vasp_scf
#BSUB -q {queue}
#BSUB -n {cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std"

start_time=$(date +%s.%N)
touch __start__
$COMMAND_std > logout 2>&1
touch __ok__
end_time=$(date +%s.%N)
runtime=$(echo "$end_time - $start_time" | bc)
dir_id=$(pwd | awk -F'/' '{{split($(NF-1),a,"_"); print a[2] "_" $NF}}')
cd ../../
echo "re_job_$dir_id total_runtime:$runtime s" >> time.txt
"""
    elif server == 'physics':
        s = f"""#!/bin/bash
#SBATCH --job-name  vasp_scf
#SBATCH --partition {queue}
##SBATCH --nodelist c0[01-40]
#SBATCH --ntasks {cores}  #number of core
##SBATCH --qos=840cpu
#####

module unload compiler/2022.1.0
module load vasp/5.4.4
COMMAND_std="mpirun vasp_std"

start_time=$(date +%s.%N)
touch __start__
$COMMAND_std > logout 2>&1
touch __ok__
end_time=$(date +%s.%N)
runtime=$(echo "$end_time - $start_time" | bc)
dir_id=$(pwd | awk -F'/' '{{split($(NF-1),a,"_"); print a[2] "_" $NF}}')
cd ../../
echo "re_job_$dir_id total_runtime:$runtime s" >> time.txt
"""
    else:
        raise ValueError(f'sever is equal to physics or qiming!')

    with open('bsub.lsf', 'w') as file:
        file.write(s)
    pwd = os.getcwd()
    #INCAR, BSUB = 'INCAR', 'bsub_1.lsf'
    BSUB = os.path.join(pwd,'bsub.lsf')
    gen_file(dir_list, INCAR_path, BSUB)
    if cal_bool:
        for a in dir_list:
            os.chdir(a)
            if server == 'physics':
                os.system('sbatch bsub.lsf')
            elif server == 'qiming':
                os.system('bsub < bsub.lsf')
    else:
        print(f'submit jobs')

def complex_bsub(server,cores,queue,ptile,INCAR_path,bsub_script=None):
    pwd = os.getcwd()
    dir_list = load_from_pkl('failed.pkl')
    result = {}
    for path in dir_list:
        # 分割路径
        parts = path.split('/')
        # 获取目录名（倒数第二个部分）
        dir_name = parts[-2]
        # 获取数字（最后一个部分）
        number = int(parts[-1])

        # 添加到字典
        if dir_name not in result:
            result[dir_name] = []
        result[dir_name].append(number)
    #INCAR, BSUB = os.path.join(pwd,'INCAR'), os.path.join(pwd,'bsub_1.lsf')
    gen_file(dir_list, INCAR_path, bsub_script)
    #print(result)
    os.chdir(pwd)
    os.makedirs('re_cal_jobs_script',exist_ok=True)
    for dir_name,num_list in result.items():
        dir_num = dir_name.split('_')[1]
        s = content(server,cores,queue,ptile,dir_name, dir_num, num_list)
        with open(os.path.join('re_cal_jobs_script',f're_bsub_{dir_num}.lsf'), 'w') as file:
            file.write(s)


if __name__ == '__main__':
    # server = 'physics'
    # cores = 56
    # queue = '256G56c'
    # ptile = 56

    server = 'qiming'
    cores = 48
    queue = '1t88c'
    ptile = 96

    import sys
    mode = sys.argv[1]

    INCAR_path = os.path.join(os.getcwd(), 'INCAR')
    if mode == 'collective':
        complex_bsub(server,cores,queue,ptile,INCAR_path,bsub_script=None)
        dir_list = load_from_pkl('failed.pkl')
        print(len(dir_list))
    elif mode == 'single':
        cal_bool = False
        pwd = os.getcwd()
        sample_bsub(server,cores,queue,ptile,INCAR_path,cal_bool)
        os.chdir(pwd)
        dir_list = load_from_pkl('failed.pkl')
        print(len(dir_list))
    else:
        raise ValueError(f'mode is equal to collective or single')
    # from ase.io import read
    # a = read(c).get_potential_energy()