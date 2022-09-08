#!/bin/bash
if [ $# -lt 2 ]
	then
		echo 'USAGE: ./measure_cutoff_flat.sh inputresultfile dt_min [cc_min] [snr_min]'
		exit 1
fi
inputfile=$1
dt_cutoff=$2
cc_cutoff=0.0
snr_cutoff=0.0
if [ $# -eq 3 ]
	then
		cc_cutoff=$3
		
fi

if [ $# -eq 4 ]
	then
		cc_cutoff=$3
		snr_cutoff=$4
fi

awk '{ if ( $2<='${dt_cutoff}' && $2>=-1.0*'${dt_cutoff}' && $9 >= '${cc_cutoff}' && $10 >= '${snr_cutoff}' ) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${inputfile}
