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

num_nodes=1  # TODO: set number of nodes requested
num_cpu=16  # TODO: set number of cores per node

echo The master node of this job is `hostname`
echo This job runs on the following nodes:
echo `scontrol show hostname $SLURM_JOB_NODELIST`

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

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
