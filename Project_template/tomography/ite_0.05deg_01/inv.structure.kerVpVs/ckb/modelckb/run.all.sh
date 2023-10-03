#!/bin/bash
if [ $# -lt 1 ]
	then
		echo './run.all.sh infilebase'
		exit 1
fi
# edit the following script before excuting run.all.csh
infilebase=$1
./inv.ckb.joint.sh ${infilebase}'.txt' 
./run.ckb.syn.sh ${infilebase}'.txt_joint'
./run.solver.1th.sh 'inv_Gd_list_'${infilebase}'.txt_joint.ori' $infilebase 
