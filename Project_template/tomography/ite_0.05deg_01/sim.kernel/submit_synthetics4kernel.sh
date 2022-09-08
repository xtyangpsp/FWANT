#!/bin/bash
#SBATCH -J synKN        #job name to remember
#SBATCH -N 1    #number of CPU cores you request for the job
#SBATCH -A xtyang  #queue to submit the job, our lab queue.
#SBATCH --mem-per-cpu 2000      #requested memory per CPU
#SBATCH -t 10:00:00                      #requested time day-hour:minute
#SBATCH -o %x.out  #path and name to save the output file.
#SBATCH -e %x.err       #path to save the error file.

echo "submitting matlab job"         

# Load module, and set up environment for Matlab to run
module load matlab

unset DISPLAY

# -nodisplay:        run MATLAB in text mode; X11 server not needed
# -singleCompThread: turn off implicit parallelism
# -r:                read MATLAB program; use MATLAB JIT Accelerator
# Run Matlab, with the above options and specifying our .m file
matlab -nodisplay -singleCompThread -r kern_synthetic4kernel

