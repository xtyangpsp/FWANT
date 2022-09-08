#!/bin/bash

#BSUB -n 120
#BSUB -R rusage[mem=2048] 
###BSUB -R "rusage[mem=2048] select[ib]"
###BSUB -R "rusage[mem=2048] select[ib] same[model] span[block=16]"
###BSUB -R "rusage[mem=2048] select[model=Opteron848]"
#BSUB -W 5 
#BSUB -J "spher.media" 
#BSUB -o outfile.%J 
#BSUB -e errorfile.%J
#BSUB -q condo_uma_haiying_gao
### BSUB -q long

module load mpich-intel/3.0.4
module load perl/5.18.1

PBS_PWD="`pwd`";
cd "${PBS_O_WORKDIR}";
THIS_HOST="`hostname`";
MPICH_ROOT="/share/pkg/mpich-intel/3.0.4";
MPIRUN_BIN="${MPICH_ROOT}/bin/mpirun";
MPIEXEC_BIN="${MPICH_ROOT}/bin/mpiexec";
FNM_BIN="./bin/seis3d_media_mpi";

P4_GLOBMEMSIZE=$(( 1024 * 1024 * 1024 * 7 ))

hr()
{
	perl -e 'print "\n" . "-"x70 . "\n\n"';
}

run_mpiexec()
{
	MPIEXEC_CMD="${MPIEXEC_BIN} ${FNM_BIN}";
    printf "begin simulation, please go to bed ...\n";
	printf "%s\n\n" "${MPIEXEC_CMD}";
	time ${MPIEXEC_CMD};
}

main()
{
	run_mpiexec;
	hr;
}

time main;

