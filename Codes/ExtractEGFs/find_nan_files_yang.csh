#!/bin/csh
#####################################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE, CROSS CORRELATION FILES ARE OCCASIONALLY
# CORRUPTED, PERPHAPS DUE TO MATLAB DUMP OR OTHER UNKNOWN PROBLEM. THIS SCRIPT FINDS AND
# REMOVES THOSE CORRUPTED FILES.
#
#
# Y.Shen, 02-28-2011
# Haiying Gao, 10-11-2016
# Cong Li, 08/07/2017
#####################################################################################
#    SET WORK DIRECTORY
#####################################################################################
set area = Alaska
set wkdir = /Users/xiaotaoyang/Work/FWANT/$area
set sacdir = $wkdir/xcorr
cd $sacdir

foreach network ( "`cat /Users/xiaotaoyang/Work/FWANT/$area/alaskanetcodes.txt`" )

foreach day (`/bin/ls -d $network* `)

cd $day
echo $day
set flag = 0

foreach sacfile ( `/bin/ls -f xc.*.SAC ` )
set dep = `echo "r $sacfile; lh DEPMEN; quit" | sac | grep "DEPMEN =" | awk '{print $3}'`
set depmin = `echo "r $sacfile; lh DEPMIN; quit" | sac | grep "DEPMIN =" | awk '{print $3}'`

if ( $dep == nan ) then
echo $sacfile $dep
#/bin/rm $sacfile
endif
if ( $depmin == 0.000000e+00 ) then
echo $sacfile $depmin
#/bin/rm $sacfile
endif

end #foreach sacfile
cd ..

end #foreach day
end # foreach network

