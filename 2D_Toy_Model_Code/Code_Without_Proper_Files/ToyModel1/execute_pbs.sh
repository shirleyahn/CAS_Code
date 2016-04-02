#!/bin/bash

#PBS -N ToyModel
#PBS -e ToyModel.err
#PBS -o ToyModel.out
#PBS -m aeb
#PBS -M sahn1@stanford.edu
#PBS -l nodes=1:ppn=1
#PBS -l walltime=168:00:00
#PBS -q long
#PBS -V

### --------------------------------------------------------------------------------------------------------
### Bash shell script is specific for PBS and GROMACS. This can be easily modified for other cluster systems
### (e.g. SLURM) and other molecular dynamics simulation programs (e.g. LAMMPS).
### --------------------------------------------------------------------------------------------------------

export MAIN_DIRECTORY=/home/sahn1/2D_Toy_Model  # TODO: set main directory
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/sahn1/

echo The master node of this job is `hostname`
echo The working directory is `echo $MAIN_DIRECTORY`
echo This job runs on the following nodes:
echo `cat $PBS_NODEFILE`

cd $MAIN_DIRECTORY

cat $PBS_NODEFILE > nodefilelist.txt
python main.py
exit
