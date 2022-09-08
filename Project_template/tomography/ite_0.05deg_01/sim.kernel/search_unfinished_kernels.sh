#!/bin/bash

if [ $# -lt 1 ]
	then
	echo 'USAGE: search_unfinished_kernels.sh stationlist'
	exit 1
fi

stafile=$1 #'station_conf_list'
nfile_min=32

for d in `awk '{print $1}' ${stafile}`
do
	if [ -d $d ]
	then

		dir=`cat ${d}_conf | grep P2 | tail -n 1 | awk '{print $5}'`
		n=`/bin/ls ${dir}/*.nc 2>/dev/null | wc -l`
		#echo ${dir} $n
		
		if [ $n -lt ${nfile_min} ]
		then
			echo $d $dir $n
		fi
	fi
done
