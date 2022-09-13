#******************************************************************************#
#*                         makefile of LSQR solver                            *#
#******************************************************************************#

FC=ifort
FFLAG=-O3
all:
	${FC} -o inv_ckb_integral ${FFLAG} inv_ckb_integral.f90
	${FC} -o inv_ckb_synthetic ${FFLAG} inv_ckb_synthetic.f90
	${FC} -o inv_make_ckb ${FFLAG} inv_make_ckb.f90

# vim:ft=make:ts=4:sw=4:nu:et:ai:
