#!/bin/bash
#SBATCH -J svr_GDLISTTEMPLATE
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A highmem
#SBATCH --mem-per-cpu 100G
#SBATCH -t 4:00:00     
#SBATCH -o %x_%A.out     
#SBATCH -e %x_%A.err
#module load intel
#module load netcdf-fortran/4.5.3

./run.solver.1th.sh ${GDLISTTEMPLATE}
