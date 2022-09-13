#!/bin/csh
####################################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE, COUNT THE NUMBER OF SAC FILES. IF THE NUMBER OF FILES IS LESS THAN A CRITERION (E.G., LESS THAN 300 DAYS.
# HERE, I SET 90 DAYS. THE DIRECTORY IS MOVE TO A TEMPT FOLDER
#
# Modified History:
#    Y.Shen, 02-28-2011
#    Haiying Gao, 10-11-2016
#    Cong Li, 08/10/2017
#	Xiaotao Yang, 11/12/2017
#			1. read netsta list directly to avoid the ls argument too long.
#####################################################################################
set area = Alaska
set wkdir = /Users/xiaotaoyang/Work/FWANT/$area/xcorr
cd $wkdir

if ( ! ( -e tmp ) ) mkdir tmp
foreach network ( "`cat /Users/xiaotaoyang/Work/FWANT/$area/alaskanetstalist.txt`" )

foreach stnpair ( `/bin/ls -d $network*` )

cd $stnpair
# @ nfile = 0
# foreach sacfile ( `/bin/ls -f *.SAC` )
# @ nfile = $nfile + 1
# end
set nfile = 0
set nfile = `/bin/ls -1 *.SAC | wc -l`
echo $stnpair $nfile
cd ..
if ( $nfile <= 90 ) then
echo  " -------------- less than 90 ----------- " $stnpair $nfile
/bin/mv $stnpair tmp/
#/bin/rm -r $stnpair
endif

end
end



