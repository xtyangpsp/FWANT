#!/bin/bash

stafile='station_conf_list'
nfile_min=32

for d in `awk '{print $1}' ${stafile}`
do
	if [ -d $d ]
	then
	
		c=`./count_kernels.sh $d | awk '{ if ( $2 < '${nfile_min}' ) print}' | wc -c`
		if [ $c -lt 1 ]
		then
			echo $d
		fi
	fi
done
