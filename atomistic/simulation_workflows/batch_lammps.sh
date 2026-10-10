#!/bin/bash

# Function to generate and submit LAMMPS EAM job
submit_eam() {
    cat << 'EOF' > sub_eam_all.sh
#!/bin/bash
#BSUB -J Cu_stress_eam_all
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD

# Define LAMMPS command
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP lmp"
# If needed, add path: export PATH=/work/phy-huangj/app/lammps/bin:$PATH

echo "Starting serial LAMMPS EAM calculations for all strains..."

# Find all strain directories under work_*
for d in work_*/strain_*; do
    if [ -d "$d" ]; then
        echo "Processing $d..."
        
        cd "$d"
        # Check if input file exists
        if [ -f "in.eam" ]; then
            # Check if data file exists
            if [ -f "data.lammps" ]; then
                $COMMAND_std -in in.eam > lammps_eam.out
            elif [ -f "data.lammps_init" ]; then
                echo "Warning: data.lammps not found, using data.lammps_init..."
                cp data.lammps_init data.lammps
                $COMMAND_std -in in.eam > lammps_eam.out
            else
                echo "Warning: data.lammps not found in $d. Skipping."
            fi
        else
            echo "Warning: in.eam not found in $d. Skipping."
        fi
        cd - > /dev/null
    fi
done

echo "All LAMMPS EAM calculations finished."
EOF

    echo "Submitting EAM job..."
    bsub < sub_eam_all.sh
}

# Function to generate and submit LAMMPS SUS2 job
submit_sus2() {
    cat << 'EOF' > sub_sus2_all.sh
#!/bin/bash
#BSUB -J Cu_stress_sus2_all
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD

# Define LAMMPS command for SUS2
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP /work/phy-luoyf/apps/lammps_sus2/interface-lammps-mlip-2/lmp_intel_cpu_intelmpi"

echo "Starting serial LAMMPS SUS2 calculations for all strains..."

for d in work_*/strain_*; do
    if [ -d "$d" ]; then
        echo "Processing $d..."
        
        cd "$d"
        if [ -f "in.sus2" ]; then
            if [ -f "data.lammps" ]; then
                $COMMAND_std -in in.sus2 > lammps_sus2.out
            elif [ -f "data.lammps_init" ]; then
                echo "Warning: data.lammps not found, using data.lammps_init..."
                cp data.lammps_init data.lammps
                $COMMAND_std -in in.sus2 > lammps_sus2.out
            else
                echo "Warning: data.lammps not found in $d. Skipping."
            fi
        else
            echo "Warning: in.sus2 not found in $d. Skipping."
        fi
        cd - > /dev/null
    fi
done

echo "All LAMMPS SUS2 calculations finished."
EOF

    echo "Submitting SUS2 job..."
    bsub < sub_sus2_all.sh
}

# Submit both
submit_eam
submit_sus2
