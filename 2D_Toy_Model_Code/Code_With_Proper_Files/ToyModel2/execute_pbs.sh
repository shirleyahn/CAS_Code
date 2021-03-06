#!/bin/bash

#PBS -N toymodel
#PBS -e toymodel.err
#PBS -o toymodel.out
#PBS -k oe
#PBS -m aeb
#PBS -M sahn1@stanford.edu
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -q default
#PBS -V

### ----------------------------------------------------------------------------------------------------------
### Bash shell script is specific for PBS. This can be easily modified for other cluster systems (e.g. SLURM).
### ----------------------------------------------------------------------------------------------------------

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sahn1/  #TODO: edit library path

echo The master node of this job is `hostname`
echo This job runs on the following nodes:
echo `cat $PBS_NODEFILE`

echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH

cd $PBS_O_WORKDIR
python main.py
exit
