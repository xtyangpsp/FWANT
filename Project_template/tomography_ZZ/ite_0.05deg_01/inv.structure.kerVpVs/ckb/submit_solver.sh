#!/bin/bash
#SBATCH -J svrckb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p cpu
#SBATCH -A xtyang
#SBATCH --mem-per-cpu 51G
#SBATCH -t 5-00:00     
#SBATCH -o %x_%A.out     
#SBATCH -e %x_%A.err
#module load intel
#module load netcdf-fortran/4.5.3

./run_series_tests.sh
