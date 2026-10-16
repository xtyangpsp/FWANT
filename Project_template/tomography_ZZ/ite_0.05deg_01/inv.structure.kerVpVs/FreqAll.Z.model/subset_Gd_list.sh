#!/bin/bash
if [ $# -lt 1 ]
	then
	echo 'USAGE: subset_Gd_list.sh dtmax'
	exit 1
fi

dtmax=$1
masterfile=inv_Gd_list

awk '{ if ( $3 <='${dtmax}' && $3 >= -1.0*'${dtmax}' ) print }' ${masterfile} > ${masterfile}_dt${dtmax}




