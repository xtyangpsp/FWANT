#!/bin/bash
#SBATCH -J svr
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A highmem
#SBATCH --mem-per-cpu 284800
#SBATCH -t 6:00:00     
#SBATCH -o %x_%A.out     
#SBATCH -e %x_%A.err
#module load intel
#module load netcdf-fortran/4.5.3
for base in layer40km
do
	./run.all.sh ${base}
done
