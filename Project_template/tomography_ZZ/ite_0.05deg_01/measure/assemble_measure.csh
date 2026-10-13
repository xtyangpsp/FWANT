#! /bin/csh
# script to assemble the phase measurements from individual station pairs

set wkdir = `pwd`
echo $wkdir

set phasedir = $wkdir/measure_stnpair
set assembled = $wkdir/ite_0.05deg_01_measure_result_craton.dat

cat /dev/null > $assembled

cd $phasedir

foreach stnpair ( `ls *.dat` )
cat $stnpair | grep f1 >> $assembled
cat $stnpair | grep f2 >> $assembled
cat $stnpair | grep f3 >> $assembled
cat $stnpair | grep f4 >> $assembled
cat $stnpair | grep f5 >> $assembled
cat $stnpair | grep f6 >> $assembled
cat $stnpair | grep f7 >> $assembled
cat $stnpair | grep f8 >> $assembled
end

