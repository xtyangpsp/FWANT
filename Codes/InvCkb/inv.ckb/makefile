#******************************************************************************#
#*                         makefile of resolution test                        *#
#******************************************************************************#

FC=ifort
FFLAG=-O3
BIN_DIR=../bin/

all:
	# Compile binaries
	${FC} -o inv_ckb_integral ${FFLAG} inv_ckb_integral.f90
	${FC} -o inv_ckb_synthetic ${FFLAG} inv_ckb_synthetic.f90
	${FC} -o inv_make_ckb ${FFLAG} inv_make_ckb.f90
	
	# Move to bin directory
	mkdir -p ${BIN_DIR}
	cp inv_ckb_integral inv_ckb_synthetic inv_make_ckb ${BIN_DIR}

clean:
	rm -f inv_ckb_integral inv_ckb_synthetic inv_make_ckb
	rm -f ${BIN_DIR}inv_ckb_integral ${BIN_DIR}inv_ckb_synthetic ${BIN_DIR}inv_make_ckb

