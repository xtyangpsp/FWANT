#!/bin/bash
if [ $# -lt 2 ]
	then
		echo 'USAGE: ./measure_cutoff_flat.sh inputresultfile cutoff'
		exit 1
fi
inputfile=$1
cutoff=$2

#awk  -F'[-,]' ' function abs(v) {return v < 0 ? -v : v} \

awk '{ if ( $2<='${cutoff}' && $2>=-1.0*'${cutoff}' ) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${inputfile}
