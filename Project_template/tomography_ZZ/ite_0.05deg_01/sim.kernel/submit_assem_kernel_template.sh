#!/bin/bash
#SBATCH -J KAJOBTEMPLATE
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -q normal
#SBATCH -A xtyang
#SBATCH -p cpu
#SBATCH --mem-per-cpu 1990
#SBATCH -t 1:00:00     
#SBATCH -o %x_%A.out     
#SBATCH -e %x_%A.err

# 1. Environment: Force a clean, Intel-specific stack
module --force purge
module load intel/2024.1 impi/2021.12 hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1

invblk='../inv.structure.kerVpVs/block.33x35x22.1x1x1.1x1x1.nc'

# This correction is not needed when both the source and sgts are 
# calculated in the same way or properly scaled
fct_ker="1 1 1 1"

# for long period waves, kernel magnitude is lower 
# Kmin for f3 is on the order of 1.0e-16 (see sim.kernel)
# from running ../inv.structure.kerVpVs/inv_make_block_stride_xygridcent.m, 
# Vm at the free surface (half thickness) ~= 2.2e14; at depth Vm ~ 9e14;
# Kmin*Vm(surf) = 2e-2
thres_ker="-3.7e-5 -3.7e-5 -3.7e-5 -3.7e-5"  # [velocity displacement P-wave S-wave]

outdir='../inv.structure.kerVpVs/G_spool/'

if [ ! -d $outdir ]; then
   mkdir -p $outdir
fi

./bin/SI_ker2assm << EOF
STATION_CONF_LIST_TEMPLATE
$invblk
$outdir
$fct_ker
$thres_ker
EOF

date
echo "finished kernel assemble"
