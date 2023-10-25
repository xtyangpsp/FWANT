#!/bin/bash

masterfile='submit_solver_bylist.sh'
subsetdirlist='subset.list'
for i in `cat ${subsetdirlist}` #`ls -d */`
do
	echo $i
	iclean=${i} #`echo $i | sed -e 's/\///g'`
	#if [ ${iclean} == "G_spool" ]; then continue; fi
	cd ${iclean}
	trycount=`ls -1 result*/try* 2>/dev/null | wc -c`
	if [ ${trycount} -eq 0 ]
	then
		echo 'submitting '${iclean}
		#sbatch -A standby -J svr_${iclean} ${masterfile}
		sbatch -A xtyang -J svr_${iclean} ${masterfile}
		sleep 1
	fi
	cd ..
done
