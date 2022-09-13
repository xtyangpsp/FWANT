#!/bin/csh
###############################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO COUNDT THE NUMBER OF
# SAC FILES IN EACH DAY, REMOVE THE DAYS THAT HAVE LESS THAN TWO STATIONS
#
# Y.Shen, 02-28-2011
# Haiying Gao, 10-11-2016
# Cong Li, 10-24-2016
################################################################
#    SET WORK DIRECTORY
###############################################
set area = Alaska
set wkdir = /Users/xiaotaoyang/Work/FWANT/$area
cd $wkdir/daily_FTN
if ( ! ( -e tmp ) ) mkdir tmp
foreach jday ( `/bin/ls -d 2*` )
cd $jday
@ nfile = 0
foreach sacfile ( `/bin/ls -f ftn.*.SAC` )
@ nfile = $nfile + 1
end
echo $jday $nfile
cd ..
if ( $nfile < 2 ) then
echo  " --------- less than 2 ------- " $jday $nfile
mv $jday tmp
#/bin/rm -r $jday
endif
end

