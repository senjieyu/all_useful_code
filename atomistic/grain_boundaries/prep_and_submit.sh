#!/bin/bash
# 这是一个提交前准备并提交作业的脚本。请在集群上运行。
# 确保安装了ase： pip install ase --user

echo "准备并提交 VASP 计算..."
cd VASP
for i in {0..9}; do
    mkdir -p $i
    cp ../${i}.vasp $i/POSCAR
    cp POTCAR 弛豫INCAR INCAR单点计算 sub.sh $i/
    cd $i
    sed -i "s/#BSUB -J .*/#BSUB -J vasp_${i}/" sub.sh
    sed -i '/calc_energy/d' sub.sh # 清理旧的遗留指令
    bsub < sub.sh
    cd ..
done
cd ..

echo "准备并提交 LAMMPS EAM 计算..."
cd LAMMPS_EAM
for i in {0..9}; do
    mkdir -p $i
    # 转换VASP POSCAR为LAMMPS DATA文件
    python3 vasp2lammps.py ../${i}.vasp $i/data.lammps
    cp in.eam sub.sh $i/
    cd $i
    sed -i "s/#BSUB -J .*/#BSUB -J eam_${i}/" sub.sh
    bsub < sub.sh
    cd ..
done
cd ..

echo "准备并提交 LAMMPS SUS2 计算..."
cd LAMMPS_SUS2
for i in {0..9}; do
    mkdir -p $i
    python3 vasp2lammps.py ../${i}.vasp $i/data.lammps
    cp in.sus2 sub.sh $i/
    cd $i
    sed -i "s/#BSUB -J .*/#BSUB -J sus2_${i}/" sub.sh
    bsub < sub.sh
    cd ..
done
cd ..

echo "所有计算作业已提交到集群！"
echo "计算完成后，请使用以下命令统计晶界能 (请替换对应的单原子体能量 E_bulk):"
echo "python3 calc_gb.py --e_eam -3.54 --e_sus2 -3.54 --e_vasp -3.72"
