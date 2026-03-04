#!/bin/bash

# --- Configuration ---
MAKEFILE_PATH="Makefile_gautschi"  # Change this if the makefile is in another directory

# Define the build modes: "ModeName|MakeArgs"
MODES=(
    "SinglePrecision|"
    "DoublePrecision|DP=1"
)

#load required modules first.
module --force purge
ml intel/2024.1 impi/2021.12 hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1

hr="----------------------------------------------------------------------"

# --- Main Build Loop ---
start_time=$SECONDS

echo "Starting synchronized build for SP and DP versions..."
echo $hr

for entry in "${MODES[@]}"; do
    # Split the entry into Name and Arguments
    IFS="|" read -r MODE_NAME ARGS <<< "${entry}"
    
    echo "BUILDING: ${MODE_NAME}"
    echo "COMMAND:  make -f ${MAKEFILE_PATH} ${ARGS}"

    # 1. Clean is mandatory to reset the #ifdef USE_DP state
    make -f "${MAKEFILE_PATH}" cleanall > /dev/null 2>&1
    
    # 2. Compile and log
    make -f "${MAKEFILE_PATH}" ${ARGS} 2>&1 | tee "make.${MODE_NAME}.log"
    
    # Check the exit status of the 'make' command (not the tee)
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "SUCCESS: ${MODE_NAME} binaries are ready."
    else
        echo "FAILURE: ${MODE_NAME} build failed. Check make.${MODE_NAME}.log"
        exit 1
    fi
    echo $hr
done

# --- Final Summary ---
duration=$(( SECONDS - start_time ))
echo "Build Process Complete."
echo "Total Time: $((duration / 60))m $((duration % 60))s"

