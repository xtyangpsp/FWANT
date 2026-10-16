#!/bin/bash
if [ $# -lt 2 ]
        then
                echo 'USAGE: ./exclude_measurements.sh meansurementfile excludelistfile'
                exit 1
fi
datafile=$1
listfile=$2
#nl=1
while read line
do
	#echo ${line}
	#echo ${nl}
	tline=`echo ${line} | sed -e 's/\// /g' | sed 's/_/ /g' `
	src=`echo ${tline} | awk '{print $4}'`
	rcv=`echo ${tline} | awk '{print $5}'`
	ftag=`echo ${tline} | awk '{print $12 }'`
	#echo ${src} ${rcv} ${ftag}
	c=0
	c=`grep ${src} ${listfile} | grep ${rcv} | grep -c ${ftag}`
	
	if [ $c == 0 ]
		then
			echo ${line} | awk '{print}'
	fi
	#(( nl = nl + 1 ))
done < ${datafile}

