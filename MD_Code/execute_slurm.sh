#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --partition=owners
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gres=gpu:8
#SBATCH --job-name=pi_hbonds
#SBATCH --output=pi_hbonds.out
#SBATCH --error=pi_hbonds.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sahn1@stanford.edu
#SBATCH --no-requeue

### ----------------------------------------------------------------------------------------------------------
### Bash shell script is specific for SLURM and GROMACS. This can be easily modified for other cluster systems
### (e.g. PBS) and other molecular dynamics simulation programs (e.g. LAMMPS).
### ----------------------------------------------------------------------------------------------------------

export MAIN_DIRECTORY=/scratch/users/sahn1/Triazine  # TODO
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sahn1/  # TODO
export GROMACS=/home/sahn1/gromacs/4.6.4/bin  # TODO
num_nodes=1  # TODO: set number of nodes requested
num_cpu=16  # TODO: set number of cores per node

echo The master node of this job is `hostname`
echo The working directory is `echo $MAIN_DIRECTORY`
echo This job runs on the following nodes:
echo `scontrol show hostname $SLURM_JOB_NODELIST`

cd $MAIN_DIRECTORY

scontrol show hostname $SLURM_JOB_NODELIST > initial_nodefilelist.txt
rm -rf nodefilelist.txt
for i in `seq 1 $num_nodes`;
do
    for j in `seq 1 $num_cpu`;
    do
        awk NR==$i initial_nodefilelist.txt >> nodefilelist.txt
    done
done
python main.py
exit
