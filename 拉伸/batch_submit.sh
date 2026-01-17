#!/bin/bash

# ================= 配置区域 =================
# LAMMPS 可执行文件路径
LAMMPS_EXE="/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp"

# 基础输入文件名
BASE_INPUT="in_cu_bulk_var.in"
# 势函数文件名
POT_FILE="Cu_mishin1.eam.alloy"
# 结构文件夹路径
STRU_DIR="stru"
# ===========================================

# 检查必要文件是否存在
if [[ ! -f "$BASE_INPUT" ]]; then
    echo "Error: Base input file $BASE_INPUT not found!"
    exit 1
fi

if [[ ! -f "$POT_FILE" ]]; then
    echo "Error: Potential file $POT_FILE not found!"
    exit 1
fi

# 遍历 stru 文件夹下的所有 .data 文件
for data_path in "${STRU_DIR}"/*.data; do
    # 检查是否有文件 (处理空文件夹情况)
    [ -e "$data_path" ] || continue

    # 获取文件名 (例如: 10.data)
    filename=$(basename "$data_path")
    # 获取任务名 (去掉后缀, 例如: 10)
    job_name="${filename%.*}"
    
    # 创建工作目录 (例如: task_10)
    work_dir="task_${job_name}"
    
    echo "Preparing job for structure: $filename -> Dir: $work_dir"
    
    # 创建目录并清理旧数据
    if [ -d "$work_dir" ]; then
        echo "  Warning: Directory $work_dir exists. Overwriting files inside."
    else
        mkdir -p "$work_dir"
    fi

    # 复制文件到工作目录
    cp "$BASE_INPUT" "$work_dir/"
    cp "$POT_FILE" "$work_dir/"
    cp "$data_path" "$work_dir/"

    # 生成该任务的提交脚本 run.sh
    cat > "$work_dir/run.sh" <<EOF
#!/bin/bash
#BSUB -J Cu_${job_name}
#BSUB -q 1t88c
#BSUB -n 96
#BSUB -R "span[ptile=96]"
#BSUB -e %J.err
#BSUB -o %J.out

# ---- environment ----
hostfile="\$LSB_DJOB_HOSTFILE"
NP=\$(wc -l < "\$hostfile")
cd "\$LS_SUBCWD"

module load mpi/2021.6.0 compiler/2022.1.0 mkl/2022.2.0

# 运行 LAMMPS
# 使用 -var struct_file 指定当前目录下的数据文件
mpirun -machinefile "\$hostfile" -np "\$NP" "$LAMMPS_EXE" -in "$BASE_INPUT" -var struct_file "$filename"

echo "Done."
EOF

    # 提交作业
    cd "$work_dir"
    # 为了防止提交太快导致调度器压力过大，稍微sleep一下（可选）
    bsub < run.sh
    
    # 返回上一级目录
    cd ..
    
    echo "  Submitted job in $work_dir"
    echo "----------------------------------------"

done

echo "All jobs submitted."

