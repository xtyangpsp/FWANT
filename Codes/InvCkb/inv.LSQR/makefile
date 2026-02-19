#******************************************************************************#
#*                         makefile of LSQR solver                            *#
#******************************************************************************#

FC=ifort #
FFLAG=-O3 #-static
APP_NAME=solver_damp_smo_loc_yang
BIN_DIR=../bin/

all:
	${FC} -o ${APP_NAME} ${FFLAG} mod_LSQR.f90 solver_damp_smo_loc_yang.f90
	mkdir -p ${BIN_DIR}
	cp ${APP_NAME} ${BIN_DIR}

clean:
	rm -f ${APP_NAME} ${BIN_DIR}${APP_NAME}

