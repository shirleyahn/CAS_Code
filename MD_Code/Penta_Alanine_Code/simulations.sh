#!/bin/bash

WALKER_DIRECTORY=$(echo $PWD)/CAS
LIBRARY_PATH=/home/sahn1/  #TODO: edit library path
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIBRARY_PATH}
export GROMACS=/home/sahn1/gromacs/4.6.4/bin  #TODO: edit MD program path

num_nodes=1  # TODO: set number of nodes requested
num_cpu=16  # TODO: set number of cores per node
#num_gpu=4
num_sim_walkers=$((num_nodes*num_cpu))

# get sequence of walker indices from sh_input.txt
output=$(cat bash_script_input_file.txt)
first_walker=$(echo ${output} | awk '{$NF=""}1')
last_walker=$(echo ${output} | cut -d'_' -f 2)

# assume walkers run in separate directories, walker0, walker1, ...
# TODO: edit MD program commands to start the simulations
cd ${WALKER_DIRECTORY}
counter=0
pin=0
for i in `seq ${first_walker} ${last_walker}`;
do
	if [[ ${counter} -lt ${num_sim_walkers} ]];
	then
		cd walker${i}
                cpu_pin=${i}%${num_cpu}
                #let gpu_pin=${i}%${num_gpu}
                # write hostfile for i-th job to use
                lstart=${counter}+1
                sed -n ${lstart},${lstart}'p' < ../../nodefilelist.txt > nodefile${counter}
                ssh $(cat nodefile${counter}) bash -c "`
                cd ${WALKER_DIRECTORY}

                cd walker${i}

                export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIBRARY_PATH}

                ${GROMACS}/grompp -f ../../prod.mdp -p ../../topol.top -c minim.gro -o run.tpr -maxwarn 1

                ${GROMACS}/mdrun -ntmpi 1 -ntomp 1 -pin on -pinstride 1 -pinoffset $cpu_pin -deffnm run
                `" &
                #${GROMACS}/mdrun -gpu_id $gpu_pin -ntmpi 1 -ntomp 1 -pin on -pinstride 1 -pinoffset $cpu_pin -deffnm run
                echo "running walker$i"
                cd ..
                let counter=counter+1
	fi
  	if [[ ${counter} == ${num_sim_walkers} ]];  # limit to num_sim_walkers concurrent walkers
  	then
      		let counter=0 # reset counter
      		wait
  	fi
done
wait

# post-process data after simulations are done. this needs to be modified depending on how the output files are named
# and what collective variables are being collected from the simulations.
# TODO: edit MD program commands to post-process data
for i in `seq 0 ${last_walker}`;
do
    cd walker${i}
    echo "post-processing walker$i"
    while [[ ! -f "run.gro" ]];  # in case a simulation failed to run
    do
        sed -n 1,1'p' < ../../nodefilelist.txt > nodefile0
        ssh $(cat nodefile0) bash -c "`
        cd ${WALKER_DIRECTORY}

        cd walker${i}

        export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LIBRARY_PATH}

        ${GROMACS}/grompp -f ../../prod.mdp -p ../../topol.top -c minim.gro -o run.tpr -maxwarn 1

        ${GROMACS}/mdrun -deffnm run
        `" &
    done
    wait
    # TODO: EDIT STARTS
    rm -rf nodefile*
    rm -rf mdout.mdp
    mv run.gro minim.gro
    mv run.tpr minim.tpr
    echo 0 | ${GROMACS}/trjconv -s minim.tpr -f run.xtc -b 50.0 -o run.xtc
    ${GROMACS}/g_angle -f run.xtc -n ../../dihedrals.ndx -ov -all -type dihedral
    awk -v f=1 -v t=2 'END {for(i=1;i<=NF;i++)if(i>=f&&i<=t)continue;else printf("%s%s",$i,(i!=NF)?OFS:ORS)}' angaver.xvg > coordinates.out
    rm -rf ang*
    if [[ -f "traj.xtc" ]];  # concatenating trajectories after the first step
    then
             ${GROMACS}/trjcat -f traj.xtc run.xtc -o out.xtc -cat
    fi
    mv run.xtc traj.xtc
    mv out.xtc traj.xtc
    rm -rf run*
    rm -rf \#*
    # TODO: EDIT ENDS
    cd ..
done
exit
