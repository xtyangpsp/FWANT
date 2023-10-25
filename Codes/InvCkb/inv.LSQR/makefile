#******************************************************************************#
#*                         makefile of LSQR solver                            *#
#******************************************************************************#

FC=mpifort #ifort #ifort #/opt/intel/composer_xe_2013_sp1.2.139/bin/intel64/ifort
FFLAG=-O3 -check all #-static

#FC=pgf90
#FFLAG=-O3 -Bstatic

all:
	#${FC} -o solver_damp_smooth ${FFLAG} mod_LSQR.f90 solver_damp_smooth.f90
	#${FC} -o solver_damp_smo_loc ${FFLAG} mod_LSQR.f90 solver_damp_smo_loc.f90
	#${FC} -o solver_damp_smo_loc_yang ${FFLAG} mod_LSQR.f90 solver_damp_smo_loc_yang.f90
	${FC} -o solver_trupin ${FFLAG} trupin.f90

# vim:ft=make:ts=4:sw=4:nu:et:ai:
