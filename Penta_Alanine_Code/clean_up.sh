#!/bin/bash

export MAIN_DIRECTORY=/scratch/users/sahn1/Penta_Alanine  # TODO
export WALKER_DIRECTORY=/scratch/users/sahn1/Penta_Alanine/CAS  # TODO
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sahn1/  # TODO
export GROMACS=/home/sahn1/gromacs/4.6.4/bin  # TODO

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
    rm -rf nodefile*
    rm -rf mdout.mdp
    mv run.gro minim.gro
    mv run.tpr minim.tpr
    echo 0 | ${GROMACS}/trjconv -s minim.tpr -f run.xtc -b 50.0 -o run.xtc
    ${GROMACS}/g_angle -f run.xtc -n ../../dihedrals.ndx -ov -all -type dihedral
    awk -v f=1 -v t=2 'END {for(i=1;i<=NF;i++)if(i>=f&&i<=t)continue;else printf("%s%s",$i,(i!=NF)?OFS:ORS)}' angaver.xvg > coordinates.out
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
