#!/bin/bash
#SBATCH -J svr
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A xtyang
#SBATCH --mem-per-cpu 19500
#SBATCH -t 4:00:00
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err
#module load intel
#module load netcdf-fortran/4.5.3

ROOT=/depot/xtyang/data/projects/noahtrupin
DATA=$ROOT/LSQR_test/example
SOLVER=$ROOT/FWANT/Codes/InvCkb/inv.LSQR/solver_trupin

time valgrind $SOLVER << EOF
$DATA/FreqAll.Z.model/inv_Gd_list_abs
2
5e-3 15
1 0
1 0
3 0
EOF
