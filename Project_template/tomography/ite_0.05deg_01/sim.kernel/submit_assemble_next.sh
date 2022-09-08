#!/bin/bash
if [ $# -lt 1 ]
	then
	echo 'USAGE: submit_assemble_next.sh unfinished_station_list'
	exit 1
fi
stalist=$1 #'unfinishedkernels_1.txt'

for sta in `awk '{print $1}' $stalist`
do
	#echo 'cleaning files'$sta
	#./clean_kernels.sh $sta
	
	echo $sta
	#sbatch -A xtyang -t 4:00:00 ${sta}.assem
	sbatch -A standby -t 3:00:00 ${sta}.assem
	sleep 0.05
done

