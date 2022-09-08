#!/bin/bash
#SBATCH -J Csim
#SBATCH -n 32
#SBATCH -A xtyang
#SBATCH --mem-per-cpu 2000
#SBATCH -t 0-5:00     
#SBATCH -o %x_%A.out     
#SBATCH -e %x_%A.err

module load intel
module load netcdf-fortran/4.5.3

PBS_PWD="`pwd`";
THIS_HOST="`hostname`";
MPICH_ROOT="/apps/spack/bell/apps/openmpi/3.1.4-intel-19.0.5-ndc76hl";
MPIRUN_BIN="${MPICH_ROOT}/bin/mpirun";
MPIEXEC_BIN="${MPICH_ROOT}/bin/mpiexec";
FNM_BIN="./bin/seis3d_wave_mpi";

print_job_info()
{
	printf "Torque Job ID: %s\n" "${SLURM_JOB_ID}";
	printf "\nRunning on host %s @ %s\n" "${THIS_HOST}" "`date`";
	printf "\nStarting directory was %s\n" "${PBS_PWD}";
	printf "Working directory is %s\n" "${SLURM_SUBMIT_DIR}";
	printf "The PWD is %s\n" "`pwd`";
	printf "\nThis job runs on the following processors:\n\n\t";
	printf "%s " `echo ${SLURM_JOB_NODELIST} | sort`;
	printf "\n\n";
	printf "This job has allocated %s nodes/processors.\n" "$SLURM_NTASKS";
}

run_mpiexec()
{
	MPIEXEC_CMD="${MPIEXEC_BIN} ${FNM_BIN}";
    	printf "begin simulation, please go to bed ...\n";
	printf "%s\n\n" "${MPIEXEC_CMD}";
	time ${MPIEXEC_CMD};
	#sleep 10;
}

main()
{
	run_mpiexec;
}

time main;

