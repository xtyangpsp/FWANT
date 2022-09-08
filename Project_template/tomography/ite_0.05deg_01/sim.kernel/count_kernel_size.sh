#!/bin/bash
stafile=$1

for d in `awk '{print $1}' ${stafile} `
do 
	cd $d
	c=`ls -d ??.* | wc -c`
	echo $d $c
	cd .. 
done

