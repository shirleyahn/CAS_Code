#!/bin/bash

###
# TODO: set main directory, walker directory, and molecular dynamics simulation program directory
# and number of nodes requested, number of cores per node, number of gpu's requested.
# num_sim_walkers should be equal to num_nodes * num_cpu, since one walker will run with one node.
###

export MAIN_DIRECTORY=/scratch/users/sahn1/Triazine
export WALKER_DIRECTORY=/scratch/users/sahn1/Triazine/CAS
export GROMACS=/home/sahn1/gromacs/4.6.4/bin

# get sequence of walker indices from sh_input.txt
cd $MAIN_DIRECTORY
output=$(cat bash_script_input_file.txt)
first_walker=$(echo $output | awk '{$NF=""}1')
last_walker=$(echo $output | cut -d'_' -f 2)

# assume walkers run in separate directories, walker0, walker1, ...
cd $WALKER_DIRECTORY
# post-process data after simulations are done. this needs to be modified depending on how the output files are named
# and what collective variables are being collected from the simulations.
for i in `seq $first_walker $last_walker`;
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

        ${GROMACS}/grompp -f ../../prod.mdp -p ../../topol_ions.top -n ../../index.ndx -c minim.gro -o run.tpr -maxwarn 1

        ${GROMACS}/mdrun -v -deffnm run
        `" &
    done
    wait
    rm -rf nodefile*
    rm -rf mdout.mdp
    mv run.gro minim.gro
    echo 0 | ${GROMACS}/trjconv -s ../../init.gro -f run.xtc -b 4.0 -o run.xtc
    #${GROMACS}/g_angle -f run.xtc -n ../../dihedrals.ndx -ov -all -type dihedral
    #awk -v f=1 -v t=2 'END {for(i=1;i<=NF;i++)if(i>=f&&i<=t)continue;else printf("%s%s",$i,(i!=NF)?OFS:ORS)}' angaver.xvg > coordinates.out
    echo 24 27 | ${GROMACS}/g_mindist -s ../../init.tpr -f run.xtc -n ../../index.ndx -nopbc
    awk 'END {print $2*10}' mindist.xvg > coordinates.out
    echo 25 28 | ${GROMACS}/g_mindist -s ../../init.gro -f run.xtc -n ../../index.ndx -nopbc
    awk 'END {print $2*10}' mindist.xvg >> coordinates.out
    echo 26 29 | ${GROMACS}/g_mindist -s ../../init.gro -f run.xtc -n ../../index.ndx -nopbc
    awk 'END {print $2*10}' mindist.xvg >> coordinates.out
    awk '{printf $0" ";next;}' coordinates.out > new_coordinates.out
    mv new_coordinates.out coordinates.out
    #rm -rf ang*
    rm -rf mindist*
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
