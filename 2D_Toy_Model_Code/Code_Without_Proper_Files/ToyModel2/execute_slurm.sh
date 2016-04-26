#!/bin/bash

#SBATCH --time=168:00:00
#SBATCH --partition=mc
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=toymodel
#SBATCH --output=toymodel.out
#SBATCH --error=toymodel.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sahn1@stanford.edu

### ----------------------------------------------------------------------------------------------------------
### Bash shell script is specific for SLURM. This can be easily modified for other cluster systems (e.g. PBS).
### ----------------------------------------------------------------------------------------------------------

echo The master node of this job is `hostname`
echo This job runs on the following nodes:
echo `scontrol show hostname $SLURM_JOB_NODELIST`

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

num_nodes=1  # TODO: set number of nodes requested
num_cpu=1  # TODO: set number of cores per node

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
