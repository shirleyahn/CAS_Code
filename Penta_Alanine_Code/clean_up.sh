#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --partition=mc
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gres=gpu:8
#SBATCH --job-name=clean_up
#SBATCH --output=clean_up.out
#SBATCH --error=clean_up.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sahn1@stanford.edu

###
# TODO: set main directory, walker directory, and molecular dynamics simulation program directory
# and number of nodes requested, number of cores per node, number of gpu's requested.
# num_sim_walkers should be equal to num_nodes * num_cpu, since one walker will run with one node.
###

export MAIN_DIRECTORY=/scratch/users/sahn1/Penta_Alanine_5
export WALKER_DIRECTORY=/scratch/users/sahn1/Penta_Alanine_5/CAS
export GROMACS=/home/sahn1/gromacs/4.6.4/bin

cd $WALKER_DIRECTORY
# post-process data after simulations are done. this needs to be modified depending on how the output files are named
# and what collective variables are being collected from the simulations.
for i in `seq 594 4819`;
do
    cd walker$i
    echo "post-processing walker$i"
    while [ ! -f "run.gro" ];  # in case a simulation failed to run
    do
        sed -n 1,1'p' < ../../nodefilelist.txt > nodefile0
        ssh $(cat nodefile0) bash -c "`
        cd ${WALKER_DIRECTORY}

        cd walker$i

        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sahn1/

        ${GROMACS}/grompp -f ../../prod.mdp -p ../../topol.top -c minim.gro -o
run.tpr -maxwarn 1

        ${GROMACS}/mdrun -v -deffnm run
        `" &
    done
    wait
    rm -rf nodefile*
    rm -rf mdout.mdp
    mv run.gro minim.gro
    mv run.tpr minim.tpr
    echo 0 | ${GROMACS}/trjconv -s minim.tpr -f run.xtc -b 50.0 -o run.xtc
    ${GROMACS}/g_angle -f run.xtc -n ../../dihedrals.ndx -ov -all -type dihedral
    awk -v f=1 -v t=2 'END {for(i=1;i<=NF;i++)if(i>=f&&i<=t)continue;else
printf("%s%s",$i,(i!=NF)?OFS:ORS)}' angaver.xvg > coordinates.out
    rm -rf ang*
    if [ -f "traj.xtc" ];  # concatenating trajectories after the first step
    then
             ${GROMACS}/trjcat -f traj.xtc run.xtc -o out.xtc -cat
    fi
    mv run.xtc traj.xtc
    mv out.xtc traj.xtc
    rm -rf run*
    rm -rf \#*
    cd ..
done
exit
