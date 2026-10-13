#!/bin/bash

for c in `ls *_conf | sed -e 's/_conf//g'`
do
	echo $c ${c}_conf
done
