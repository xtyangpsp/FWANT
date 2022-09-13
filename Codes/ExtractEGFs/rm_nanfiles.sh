#!/bin/bash
workdir="xcorr"
listfile="xcor_nan_files.txt"
xcorr_nan_dir="/Users/xiaotaoyang/Work/FWANT/Alaska/xcorr_nan_files"
mkdir ${xcorr_nan_dir}
cd $workdir
for netsta in `awk '{print $1}' ../$listfile| sed 's/\./ /g' | awk '{print $4"."$5"."$6"."$7}'`
do
# 	echo $netsta
	mv -v $netsta/xc.2007.016.*.SAC ${xcorr_nan_dir}
done

cd ..
