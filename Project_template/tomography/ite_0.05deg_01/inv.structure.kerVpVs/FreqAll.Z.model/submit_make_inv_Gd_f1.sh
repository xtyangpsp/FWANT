#!/bin/bash
#SBATCH -J Gd.f1
#SBATCH -n 1
#SBATCH -A xtyang
#SBATCH -t 5:00:00     
#SBATCH -o %x_%A.out     
#SBATCH -e %x_%A.err

./make_inv_Gd_list.sh f1
