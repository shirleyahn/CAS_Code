#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --partition=mc
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gres=gpu:8
#SBATCH --job-name=9_40_10ps
#SBATCH --output=9_40_10ps.out
#SBATCH --error=9_40_10ps.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sahn1@stanford.edu

###
# Bash shell script is specific for SLURM and GROMACS. This can be easily modified for other cluster systems
# (e.g. PBS) and other molecular dynamics simulation programs (e.g. LAMMPS).
###

echo The master node of this job is `hostname`
echo The working directory is `echo $WORKDIR`
echo This job runs on the following nodes:
echo `scontrol show hostname $SLURM_JOB_NODELIST`

export MAIN_DIRECTORY=/scratch/users/sahn1/WE_Triazine  # TODO: set main directory for WE simulation
num_nodes=1  # TODO: set number of nodes requested
num_cpu=16  # TODO: set number of cores per node

cd $MAIN_DIRECTORY
scontrol show hostname $SLURM_JOB_NODELIST > initial_nodefilelist.txt
for i in `seq 1 $num_nodes`;
do 
    for j in `seq 1 $num_cpu`;
    do
        awk NR==$i initial_nodefilelist.txt >> nodefilelist.txt
    done
done
python we_main.py
exit
