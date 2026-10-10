#!/bin/bash
#BSUB -J eam_relax_sp
#BSUB -q 33
#BSUB -n 40
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile=40]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
# export PATH=/work/phy-huangj/app/lammps/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP lmp"

echo "Starting LAMMPS EAM relaxation and single point calculation..."

echo "Calculating energy..."
# If there's a calc energy script for lammps later

echo 'Starting serial calculations...'
for d in calc_*; do
    if [ -d "$d" ]; then
        echo "Running $d"
        cd $d
        $COMMAND_std -in in.eam > lammps.out
        cd ..
    fi
done
