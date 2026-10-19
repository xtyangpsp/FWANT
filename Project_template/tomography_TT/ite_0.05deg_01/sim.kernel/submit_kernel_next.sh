#!/bin/bash
if [ $# -lt 1 ]
	then
	echo 'USAGE: submit_kernels_next.sh unfinished_station_list'
	exit 1
fi
stalist=$1 #'unfinishedkernels_1.txt'
#startline=1
#endline=68
echo 'start, end' ${startline} ${endline}
#for sta in `awk '{ if ( NR >= '${startline}' && NR <= '${endline}' ) print $1}' $stalist`
for sta in `awk '{ print $1 }' $stalist`
do
	echo 'submitting '$sta
#	sbatch ${sta}.submit
	sbatch -q standby -t 4:00:00 ${sta}.submit
	sleep 0.1
done

