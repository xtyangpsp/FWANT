#!/bin/bash
#BSUB -n 120
#BSUB -R rusage[mem=2048] 
#BSUB -W 1200 
#BSUB -J "AKsim0.05v1.wave" 
#BSUB -o outfile.%J 
#BSUB -e errorfile.%J
#BSUB -q condo_uma_haiying_gao

module load mpich-intel/3.0.4
module load perl/5.18.1

PBS_PWD="`pwd`";
cd "${PBS_O_WORKDIR}";
THIS_HOST="`hostname`";
MPICH_ROOT="/share/pkg/mpich-intel/3.0.4";
MPIRUN_BIN="${MPICH_ROOT}/bin/mpirun";
##MPIEXEC_BIN="/share/pkg/mpich/3.0.4/bin/mpiexec";
MPIEXEC_BIN="${MPICH_ROOT}/bin/mpiexec";
FNM_BIN="./bin/seis3d_wave_mpi";

#NPROCS="`wc -l < ${PBS_NODEFILE} | tr -d '[:blank:]'`";

P4_GLOBMEMSIZE=$(( 1024 * 1024 * 1024 * 7 ))

hr()
{
	perl -e 'print "\n" . "-"x70 . "\n\n"';
}

print_job_info()
{
	hr;
	printf "Torque Job ID: %s\n" "${PBS_JOBID}";
	printf "\nRunning on host %s @ %s\n" "${THIS_HOST}" "`date`";
	printf "\nStarting directory was %s\n" "${PBS_PWD}";
	printf "Working directory is %s\n" "${PBS_O_WORKDIR}";
	printf "The PWD is %s\n" "`pwd`";
	printf "\nThis job runs on the following processors:\n\n\t";
	printf "%s " `cat ${PBS_NODEFILE} | sort`;
	printf "\n\n";
	printf "This job has allocated %s nodes/processors.\n" "${NPROCS}";
	hr;
}

clean_ipcs()
{
	for NODE in `cat ${PBS_NODEFILE} | sort -u`; do
		ssh ${NODE} ${CLEANIPCS_BIN};
	done;
}

run_mpirun()
{
	clean_ipcs;
	MPI_CMD="${MPIRUN_BIN} -nolocal -np ${NPROCS} -machinefile ${PBS_NODEFILE} ${FNM_BIN}";
    printf "begin simulation, please go to bed ...\n";
	printf "%s\n\n" "${MPI_CMD}";
	time ${MPI_CMD};
	#sleep 10;
}

run_mpiexec()
{
	#clean_ipcs;
	#MPIEXEC_CMD="${MPIEXEC_BIN} --comm p4 -mpich-p4-no-shmem ${FNM_BIN}";
	MPIEXEC_CMD="${MPIEXEC_BIN} ${FNM_BIN}";
    printf "begin simulation, please go to bed ...\n";
	printf "%s\n\n" "${MPIEXEC_CMD}";
	time ${MPIEXEC_CMD};
	#sleep 10;
}

main()
{
	#print_job_info;
	#run_mpirun;
	run_mpiexec;
	hr;
}

time main;

