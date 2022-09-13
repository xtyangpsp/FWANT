#!/bin/bash
if [ $# -lt 1 ]
        then
                echo 'USAGE: ./reformatNEIC_EQinfo.sh rawNEICfile'
                echo '** rawNEICfile: has to be in csv format when downloaded from NEIC.'
                echo 'By Xiaotao Yang @ UMass Amherst, 2018'
                exit 1
fi

rawfile=$1

sed -e 's/\,/ /g' ${rawfile} | awk '{print $1,$2,$3,$4,$5,$6}' | sed -e 's/T/ /g' | sed -e 's/Z/ /g'
