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

# echo "echo off" > $sacmacro

#foreach sacfile ( `/bin/ls ftn.*.$station.*.SAC` ) # UW.*.EHZ 2005-2011
foreach sacfile ( `/bin/ls ftn.*.SAC` )
#echo "echo off" > $sacmacro
sac <<END
echo off
cuterr fillz
cut 0 86400
r $sacfile
setbb oldb &1,b&
setbb tmpjday &1,nzjday&
ch nzjday $jday
setbb truejday &1,nzjday&
evaluate to jdaydiff %tmpjday - %truejday
setbb tmphour &1,nzhour&
setbb tmpmin &1,nzmin&
setbb tmpsec &1,nzsec&
setbb tmpmsec &1,nzmsec&
evaluate to tmp %tmpmsec * 0.001
evaluate to tmp2 %tmphour * 3600
evaluate to tmp3 %tmpmin * 60
evaluate to tmp4 %jdaydiff * 86400
evaluate to diff %tmpsec + %tmp + %tmp2 + %tmp3 + %tmp4
evaluate to newb %oldb + %diff
ch b %newb
ch nzhour 0
ch nzmin 0
ch nzsec 0
ch nzmsec 0
interpolate d 0.5 b 0
w $sacfile
cuterr fillz
cut 0 86400
r $sacfile
w $sacfile
quit
END
end #foreach sacfile
# /bin/rm $sacmacro
cd ..
end #foreach day
#end

