#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --partition=owners
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=2D_Toy_Model
#SBATCH --output=2D_Toy_Model.out
#SBATCH --error=2D_Toy_Model.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sahn1@stanford.edu

### ----------------------------------------------------------------------------------------------------------
### Bash shell script is specific for SLURM and GROMACS. This can be easily modified for other cluster systems
### (e.g. PBS) and other molecular dynamics simulation programs (e.g. LAMMPS).
### ----------------------------------------------------------------------------------------------------------

export MAIN_DIRECTORY=/scratch/users/sahn1/2D_Toy_Model  # TODO: set main directory for CAS simulation
num_nodes=1  # TODO: set number of nodes requested
num_cpu=1  # TODO: set number of cores per node

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