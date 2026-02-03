#!/bin/bash

# ================= 配置区域 =================
# LAMMPS 可执行文件路径
LAMMPS_EXE="/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp"

# 输入文件定义
RELAX_INPUT="relax.in"
TENSILE_INPUT="tensile.in"
ANALYZE_SCRIPT="analyze.py"
POT_FILE="Cu_mishin1.eam.alloy"
STRU_DIR="stru"
# ===========================================

# 检查必要文件是否存在
for f in "$RELAX_INPUT" "$TENSILE_INPUT" "$ANALYZE_SCRIPT" "$POT_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: File $f not found!"
        exit 1
    fi
done

# 创建全局结果保存目录
if [ ! -d "stru_relaxed" ]; then mkdir -p stru_relaxed; fi
if [ ! -d "results_X" ]; then mkdir -p results_X; fi
if [ ! -d "results_Y" ]; then mkdir -p results_Y; fi
if [ ! -d "results_Z" ]; then mkdir -p results_Z; fi

# 初始化汇总文件
SUMMARY_FILE="$(pwd)/yield_points_summary.txt"
echo "# Yield Point Analysis Results - $(date)" > "$SUMMARY_FILE"

# 遍历 stru 文件夹下的所有 .data 文件
for data_path in "${STRU_DIR}"/*.data; do
    # 检查是否有文件
    [ -e "$data_path" ] || continue

    # 获取文件名和任务名
    filename=$(basename "$data_path")
    job_name="${filename%.*}"
    work_dir="task_${job_name}"
    
    echo "Preparing job for structure: $filename -> Dir: $work_dir"
    
    # 重新创建工作目录
    if [ -d "$work_dir" ]; then
        echo "  Cleaning existing directory..."
        rm -rf "$work_dir"
    fi
    mkdir -p "$work_dir"

    # 复制文件到工作目录
    cp "$RELAX_INPUT" "$work_dir/"
    cp "$TENSILE_INPUT" "$work_dir/"
    cp "$ANALYZE_SCRIPT" "$work_dir/"
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
if [ -z "\$hostfile" ]; then
    # Fallback if not running in LSF
    NP=4
    CMD="mpirun -np \$NP"
else
    NP=\$(wc -l < "\$hostfile")
    CMD="mpirun -machinefile \$hostfile -np \$NP"
fi
cd "\$LS_SUBCWD"

module load mpi/2021.6.0 compiler/2022.1.0 mkl/2022.2.0

echo "=========================================="
echo "Job Start: \$(date)"
echo "Structure: $job_name"
echo "=========================================="

# ----------------------------------------
# 1. 全弛豫 (Relaxation)
# ----------------------------------------
echo "[Step 1] Running Relaxation (ISIF=3 equivalent)..."
\$CMD "$LAMMPS_EXE" -in "$RELAX_INPUT" -var struct_file "$filename" > log.relax

if [ ! -f relaxed.data ]; then
    echo "Error: Relaxation failed. relaxed.data not found."
    exit 1
fi

# 保存弛豫后的结构到全局目录
cp relaxed.data ../stru_relaxed/${job_name}_relaxed.data
echo "Relaxed structure saved to stru_relaxed/${job_name}_relaxed.data"

# ----------------------------------------
# 2. 拉伸 (Tensile) - X, Y, Z
# ----------------------------------------

# 定义拉伸函数
run_tensile() {
    local dir=\$1
    echo "[Step 2-\$dir] Running Tensile Test along \$dir direction..."
    
    # 创建子目录以防文件冲突
    mkdir -p \$dir
    cd \$dir
    
    # 链接必要文件
    ln -s ../relaxed.data .
    ln -s ../$POT_FILE .
    
    # 运行 LAMMPS
    \$CMD "$LAMMPS_EXE" -in ../"$TENSILE_INPUT" -var struct_file relaxed.data -var tensile_dir \$dir > log.tensile_\$dir
    
    # 运行 Python 分析
    # 输出: png图, 结果txt
    python ../"$ANALYZE_SCRIPT" stress_strain_\$dir.dat stress_strain_\$dir.png yield_info.txt \$dir
    
    # 整理结果
    # 复制到全局目录
    GLOBAL_RES_DIR="../../results_\$(echo \$dir | tr '[:lower:]' '[:upper:]')"
    cp stress_strain_\$dir.png "\$GLOBAL_RES_DIR/${job_name}_stress_strain.png"
    cp stress_strain_\$dir.dat "\$GLOBAL_RES_DIR/${job_name}_stress_strain.dat"
    cp final_structure_\$dir.data "\$GLOBAL_RES_DIR/${job_name}_final.data"
    
    # 汇总屈服点信息
    echo "Structure: ${job_name} | Direction: \$dir" >> "$SUMMARY_FILE"
    cat yield_info.txt >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
    
    cd ..
}

# 依次执行三个方向
run_tensile x
run_tensile y
run_tensile z

echo "=========================================="
echo "Job Done: \$(date)"
echo "=========================================="
EOF

    # 提交作业
    cd "$work_dir"
    bsub < run.sh
    cd ..
    
    echo "  Submitted job in $work_dir"
    echo "----------------------------------------"

done

echo "All jobs submitted."
