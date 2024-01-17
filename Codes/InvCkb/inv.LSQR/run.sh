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
SOLVER=$ROOT/FWANT/Codes/InvCkb/inv.LSQR/solver_yang
HOSTFILE=hostfile
NPROC=1
BINDIR=bin
TARGET=$BINDIR/solver

OMP_NUM_THREADS=8 mpirun --hostfile ./$HOSTFILE -np $NPROC ./$TARGET << EOF
1
$DATA/FreqAll.Z.model/inv_Gd_list_abs
$DATA/block.46x79x43.1x1x1.1x1x1.smooth1th.dat
2
1 0
1 0
2 0
8 8
5e-3 15
EOF
