#******************************************************************************#
#*                         makefile of LSQR solver                            *#
#******************************************************************************#

FC=ifort
all:
	${FC} -o inv_ckb_integral -O3  -static inv_ckb_integral.f90
	${FC} -o inv_ckb_synthetic -O3  -static inv_ckb_synthetic.f90
	${FC} -o inv_make_ckb -O3  -static inv_make_ckb.f90

# vim:ft=make:ts=4:sw=4:nu:et:ai:
