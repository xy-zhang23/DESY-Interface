#!/bin/bash

# one process for each physical core:
#SBATCH -c 2

# run on the haswell partition (pax11):
#SBATCH -p haswell
##SBATCH -p sandybridge

# runtime 
#SBATCH --time=47:59:49

# number of tasks/cores
#SBATCH --ntasks=16

# Job name:
#SBATCH --job-name=test
#SBATCH --error=test-%N-%j.err
#SBATCH --output=test-%N-%j.out

# copy environment variables from submit environment
#SBATCH --get-user-env

# send mail on all occasions:
#SBATCH --mail-type=ALL

#export PATH=$PATH:$(pwd)
#module add intel openmpi

#. /afs/ifh.de/group/pitz/data/lixiangk/work/apps/anaconda2/bin/activate
module add python/3
module add gnu8 openmpi3

which python

mpirun python injector_optimization_demo.py
#mpirun -np 101 python simple_mpi.py
