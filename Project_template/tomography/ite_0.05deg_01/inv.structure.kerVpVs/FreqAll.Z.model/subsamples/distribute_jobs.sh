#!/bin/bash

#setup needed folders and files.
ln -s ../../G_spool .
ln -s ../../block* .

masterfile='submit_solver_bylist.sh'
gdlist='inv_Gd_list.txt'
subsetlist='subset.list'

for i in `cat ${subsetlist}` #`ls -d 1/`
do
	echo $i
	cd $i
	#sed -e 's/\..G_spool/G_spool/g' ${gdlist} > inv_Gd_list
	#mkdir G_spool
	#for listtemp in `awk '{print $1" " $2}' ${gdlist}`
	#do
	#	cp ${listtemp} G_spool/
	#done
	ln -s ../${masterfile} .
	cd ..
done
