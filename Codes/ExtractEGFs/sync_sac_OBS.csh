#!/bin/csh
##########################################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO SYNCHRONIZE SAC FILES SO TRACES START AT
# EXACTLY THE SAME TIME
#
#
# Y.Shen, 02-28-2011
# Haiying Gao, 10-11-2016
# Cong Li, 10-24-2016
#     adding line 57, 58 to fill 0 so that the length of data will be equal.
##########################################################################################
#           SET WORK DIRECTORY
#################################################
set area = Alaska
set wkdir = /Users/xiaotaoyang/Work/FWANT/$area

set sacdir = $wkdir/daily_FTN
cd $sacdir

set sacmacro = sync.m

foreach day ( `/bin/ls -d 2* ` )
cd $day
echo $day

set jday = `echo $day | awk -F. '{print $2}'`

echo "echo off" > $sacmacro

#foreach sacfile ( `/bin/ls ftn.*.$station.*.SAC` ) # UW.*.EHZ 2005-2011
foreach sacfile ( `/bin/ls ftn.*.SAC` )
#echo "echo off" > $sacmacro
echo "cuterr fillz" >> $sacmacro
echo "cut 0 86400" >> $sacmacro
echo "r $sacfile" >> $sacmacro
echo "setbb oldb &1,b&" >> $sacmacro
echo "setbb tmpjday &1,nzjday&" >> $sacmacro # GMT time; The time in some station is the day before the jday, fast several seconds, sync to jday, 0 hour, 0 min, 0 seconds;
echo "ch nzjday $jday" >> $sacmacro
echo "setbb truejday &1,nzjday&" >> $sacmacro
echo "evaluate to jdaydiff %tmpjday - %truejday" >> $sacmacro
echo "setbb tmphour &1,nzhour&" >> $sacmacro
echo "setbb tmpmin &1,nzmin&" >> $sacmacro
echo "setbb tmpsec &1,nzsec&" >> $sacmacro
echo "setbb tmpmsec &1,nzmsec&" >> $sacmacro
echo "evaluate to tmp %tmpmsec * 0.001" >> $sacmacro
echo "evaluate to tmp2 %tmphour * 3600" >> $sacmacro
echo "evaluate to tmp3 %tmpmin * 60" >> $sacmacro
echo "evaluate to tmp4 %jdaydiff * 86400" >> $sacmacro
echo "evaluate to diff %tmpsec + %tmp + %tmp2 + %tmp3 + %tmp4" >> $sacmacro
echo "evaluate to newb %oldb + %diff" >> $sacmacro
echo "ch b %newb" >> $sacmacro
echo "ch nzhour 0" >> $sacmacro
echo "ch nzmin 0" >> $sacmacro
echo "ch nzsec 0" >> $sacmacro
echo "ch nzmsec 0" >> $sacmacro
echo "interpolate d 0.5 b 0" >> $sacmacro #
echo "cuterr fillz" >> $sacmacro # if
echo "cut 0 86400" >> $sacmacro
echo "w $sacfile" >> $sacmacro
end # foreach sacfile
echo "quit" >> $sacmacro
sac <<END
m $sacmacro
END
#end #foreach sacfile
/bin/rm $sacmacro
cd ..
end #foreach day
#end

