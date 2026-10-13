#!/bin/bash
#SBATCH -J KNJOBTEMPLATE 
#SBATCH -n 105
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH -A xtyang
#SBATCH -N 1
#SBATCH --mem-per-cpu 1990
#SBATCH -t 1:00:00     
#SBATCH -o %x_%A.out     
#SBATCH -e %x_%A.err

# 1. Environment: Force a clean, Intel-specific stack
module --force purge
module load intel/2024.1 impi/2021.12 hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1

# 2. Performance & Stability Flags
# Set fabric log level to only show critical errors
export FI_LOG_LEVEL=error
export FI_PROVIDER=psm3
export SLURM_CPU_BIND=none
ulimit -s unlimited

# 3. Execution: Using srun for the best Slurm-Intel integration
echo "Job started at: $(date) on host $(hostname)"
echo "Allocated tasks: $SLURM_NTASKS"

# 4. run the program
time mpirun -n ${SLURM_NTASKS} ./bin/SI_ker_sta_mpi << EOF
STATION_CONF_LIST_TEMPLATE
EOF

