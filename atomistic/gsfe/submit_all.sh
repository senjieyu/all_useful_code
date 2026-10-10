#!/bin/bash

# 延迟时间（秒），防止瞬时提交过多任务给调度器造成压力
DELAY=0.5

echo "==========================================="
echo "  Batch Submission Script for GSFE curves"
echo "==========================================="

# 提交 VASP 任务
echo -e "\n---> Submitting VASP jobs..."
for d in vasp/calc_*/; do
    if [ -d "$d" ]; then
        echo "Submitting job in $d"
        cd "$d"
        bsub < sub.sh
        cd - > /dev/null
        sleep $DELAY
    fi
done

# 提交 LAMMPS EAM 任务 (单节点)
echo -e "\n---> Submitting LAMMPS EAM job..."
cd lammps_eam
bsub < sub.sh
cd - > /dev/null
sleep $DELAY

# 提交 LAMMPS SUS2 任务 (单节点)
echo -e "\n---> Submitting LAMMPS SUS2 job..."
cd lammps_sus2
bsub < sub.sh
cd - > /dev/null
sleep $DELAY

echo -e "\n==========================================="
echo "All jobs have been submitted successfully!"
echo "Use 'bjobs' to check your job status."
