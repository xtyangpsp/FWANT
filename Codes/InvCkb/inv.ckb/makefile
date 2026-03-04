######################################
### ckb codes
##Default is in SP (single precision, for kernal reading).
## To compile in DP (double precision), run: make DP=1
####################################

FC=ifort
BIN_DIR=../bin/

# Base Flags with Pre-processor support
FFLAG_BASE=-O3 -xHost -heap-arrays -assume byterecl -fpp

# Precision Toggle Logic
ifeq ($(DP), 1)
    FFLAG = $(FFLAG_BASE) -DUSE_DP
    NAME_SYN = inv_ckb_synthetic_DP
    NAME_INT = inv_ckb_integral_DP
else
    FFLAG = $(FFLAG_BASE)
    NAME_SYN = inv_ckb_synthetic_SP
    NAME_INT = inv_ckb_integral_SP
endif

all:
	# Compile binaries with the selected precision
	${FC} -o ${NAME_INT} ${FFLAG} inv_ckb_integral.f90
	${FC} -o ${NAME_SYN} ${FFLAG} inv_ckb_synthetic.f90
	${FC} -o inv_make_ckb ${FFLAG_BASE} inv_make_ckb.f90

	# Move to bin directory
	mkdir -p ${BIN_DIR}
	cp ${NAME_INT} ${NAME_SYN} inv_make_ckb ${BIN_DIR}

# Clean target to remove local and installed binaries
clean:
	rm -f ${NAME_SYN} ${NAME_INT}
	rm -f inv_make_ckb
	rm -f ${BIN_DIR}${NAME_SYN} ${BIN_DIR}${NAME_INT}
	rm -f ${BIN_DIR}inv_make_ckb

