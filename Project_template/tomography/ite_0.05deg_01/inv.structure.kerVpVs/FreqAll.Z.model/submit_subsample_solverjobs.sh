#!/bin/bash
mastersubmitfile='submit_solver_bylist.sh'
for f in `ls -1 subsamples/inv_Gd*`
do
	##
	echo 'Submitting '$f
	sed -e 's/GDLISTTEMPLATE/"${f}"/g' | sbatch
done

