#!/bin/bash
if [ $# -lt 1 ]
	then
	echo 'USAGE: redo_unfinished_kernels.sh stationlist'
	exit 1
fi
stalist=$1
nmin=32
for s in `awk '{print $1}' $stalist`
do
	#./XO.MG09/XO.LD20/BHZ/7/T1T2.P2 32
	ev=`./count_kernels.sh $s | awk '{ if ( $2 < '$nmin' ) print $1 }' | head -n 1 | sed -e 's/\// /g' | awk '{print $3}'`
	echo $s $ev

	lnum=`awk '{ if ( $1 == "'${ev}'" ) print NR}' ${s}_conf`
	cp ${s}_conf ${s}_conf_bkp
	
	awk '{ if ( NR >= '${lnum}' ) print}' ${s}_conf_bkp > ${s}_conf

	#submit jobs
	echo 'submitting jobs for: '$s
	sbatch -A standby -t 4:00:00 ${s}.submit
done

#change name back
for f in `ls *_bkp`; do nf=`echo $f | sed -e 's/_bkp//g'`; echo $nf;mv $f $nf; done

