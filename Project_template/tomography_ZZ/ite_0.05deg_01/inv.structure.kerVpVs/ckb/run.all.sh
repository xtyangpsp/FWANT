#!/bin/bash
if [ $# -lt 6 ]
	then
		echo './run.all.sh xc yc zc xwd ywd zwd [dv]'
		exit 1
fi
if [ $# -eq 7 ]
then
        dvmax=$7
else
        dvmax=0.1 #default 0.1
fi

xc=$1
yc=$2
zc=$3
xwd=$4
ywd=$5
zwd=$6

# edit the following script before excuting run.all.csh

./inv.ckb.make.sh ${xc} ${yc} ${zc} ${xwd} ${ywd} ${zwd} ${dvmax}
./inv.ckb.joint.sh ${xwd} ${ywd} ${zwd} ${dvmax}
./run.ckb.syn.sh ${xwd} ${ywd} ${zwd} ${dvmax}
./run.solver.1th.sh ${xwd} ${ywd} ${zwd} ${dvmax}
