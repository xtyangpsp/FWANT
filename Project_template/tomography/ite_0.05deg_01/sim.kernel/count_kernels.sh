#!/bin/bash
if [ $# -lt 1 ]
	then
	echo 'USAGE: count_kernels.sh station'
	exit 1
fi
sta=$1
for dir in `cat ${sta}_conf | grep T1T2.P2 | awk '{print $5}'`
do 
	n=`/bin/ls ${dir}/*.nc 2>/dev/null| wc -l`
	echo ${dir} $n
done
