#!/bin/bash
if [ $# -lt 1 ]
	then
	echo 'USAGE: clean_kernels.sh station'
	exit 1
fi
sta=$1
for dir in `cat ${sta}_conf | grep P2 | awk '{print $5}'`
do 
	/bin/rm ${dir}/*.nc 
	echo 'cleaned :' ${dir}; 
done
