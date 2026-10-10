#!/bin/bash

set -euo pipefail

LAMMPS_EXE="/work/phy-luoyf/apps/lammps2025/lammps-22Jul2025/build/lmp"

RELAX_INPUT="relax.in"
SHEAR_INPUT="tensile.in"
ANALYZE_SCRIPT="analyze.py"
POT_FILE="Cu_mishin1.eam.alloy"
STRU_DIR="stru"
TARGET_STRUCT="POSCAR_Cu111_twin_singleplane_d1_7_d2_12_d3_5_1100.vasp_relaxed.data"

for f in "$RELAX_INPUT" "$SHEAR_INPUT" "$ANALYZE_SCRIPT" "$POT_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: missing required file: $f"
        exit 1
    fi
done

rm -rf stru_relaxed results_XY results_XZ results_YZ task_*
mkdir -p stru_relaxed results_XY results_XZ results_YZ

SUMMARY_FILE="$(pwd)/shear_summary.txt"
echo "# Shear Analysis Results - $(date)" > "$SUMMARY_FILE"

if [[ ! -f "${STRU_DIR}/${TARGET_STRUCT}" ]]; then
    echo "Error: target structure not found: ${STRU_DIR}/${TARGET_STRUCT}"
    exit 1
fi

for data_path in "${STRU_DIR}/${TARGET_STRUCT}"; do
    [ -e "$data_path" ] || continue

    filename=$(basename "$data_path")
    job_name="${filename%.*}"
    work_dir="task_${job_name}"

    echo "Preparing shear jobs for structure: $filename -> $work_dir"

    rm -rf "$work_dir"
    mkdir -p "$work_dir"

    cp "$RELAX_INPUT" "$work_dir/"
    cp "$SHEAR_INPUT" "$work_dir/"
    cp "$ANALYZE_SCRIPT" "$work_dir/"
    cp "$POT_FILE" "$work_dir/"
    cp "$data_path" "$work_dir/"

    cat > "$work_dir/run_relax.sh" <<EOF
#!/bin/bash
#BSUB -J Cu_${job_name}_relax
#BSUB -q v3-64
#BSUB -n 26
#BSUB -R "span[ptile=24]"
#BSUB -e relax_%J.err
#BSUB -o relax_%J.out

set -euo pipefail

hostfile="\$LSB_DJOB_HOSTFILE"
if [ -z "\$hostfile" ]; then
    NP=4
    CMD="mpirun -np \$NP"
else
    NP=\$(wc -l < "\$hostfile")
    CMD="mpirun -machinefile \$hostfile -np \$NP"
fi

cd "\$LS_SUBCWD"
module load mpi/2021.6.0 compiler/2022.1.0 mkl/2022.2.0

echo "=========================================="
echo "Relax Job Start: \$(date)"
echo "Structure: $job_name"
echo "=========================================="

\$CMD "$LAMMPS_EXE" -in "$RELAX_INPUT" -var struct_file "$filename" > log.relax

if [ ! -f relaxed.data ]; then
    echo "Error: relaxed.data not found after relaxation."
    exit 1
fi

cp relaxed.data ../stru_relaxed/${job_name}_relaxed.data

echo "=========================================="
echo "Relax Job Done: \$(date)"
echo "=========================================="
EOF

    for mode in xy xz yz; do
        upper_mode=$(echo "$mode" | tr '[:lower:]' '[:upper:]')
        cat > "$work_dir/run_shear_${mode}.sh" <<EOF
#!/bin/bash
#BSUB -J Cu_${job_name}_${mode}
#BSUB -q v3-64
#BSUB -n 26
#BSUB -R "span[ptile=24]"
#BSUB -e shear_${mode}_%J.err
#BSUB -o shear_${mode}_%J.out

set -euo pipefail

hostfile="\$LSB_DJOB_HOSTFILE"
if [ -z "\$hostfile" ]; then
    NP=4
    CMD="mpirun -np \$NP"
else
    NP=\$(wc -l < "\$hostfile")
    CMD="mpirun -machinefile \$hostfile -np \$NP"
fi

cd "\$LS_SUBCWD"
module load mpi/2021.6.0 compiler/2022.1.0 mkl/2022.2.0

echo "=========================================="
echo "Shear Job Start: \$(date)"
echo "Structure: $job_name | Mode: $mode"
echo "=========================================="

mkdir -p $mode
cd $mode

ln -sf ../relaxed.data .
ln -sf ../$POT_FILE .

\$CMD "$LAMMPS_EXE" -in ../"$SHEAR_INPUT" -var struct_file relaxed.data -var shear_mode $mode > log.shear_$mode

python ../"$ANALYZE_SCRIPT" stress_strain_$mode.dat stress_strain_$mode.png yield_info.txt $mode

GLOBAL_RES_DIR="../../results_${upper_mode}"
cp stress_strain_$mode.png "\$GLOBAL_RES_DIR/${job_name}_stress_strain.png"
cp stress_strain_$mode.dat "\$GLOBAL_RES_DIR/${job_name}_stress_strain.dat"
cp final_structure_$mode.data "\$GLOBAL_RES_DIR/${job_name}_final.data"

echo "Structure: ${job_name} | Mode: $mode" >> "$SUMMARY_FILE"
cat yield_info.txt >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

echo "=========================================="
echo "Shear Job Done: \$(date)"
echo "=========================================="
EOF
    done

    cd "$work_dir"
    relax_res=$(bsub < run_relax.sh)
    echo "  $relax_res"
    relax_job_id=$(echo "$relax_res" | awk -F'[<>]' '{print $2}')

    if [ -n "$relax_job_id" ]; then
        for mode in xy xz yz; do
            bsub -w "done($relax_job_id)" < "run_shear_${mode}.sh"
        done
    else
        echo "  Warning: could not parse relax job ID, submitting shear jobs without dependency."
        for mode in xy xz yz; do
            bsub < "run_shear_${mode}.sh"
        done
    fi

    cd ..
    echo "----------------------------------------"
done

echo "All shear jobs submitted."
