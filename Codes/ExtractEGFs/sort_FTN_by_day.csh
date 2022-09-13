#!/bin/csh
#########################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO SORT SAC FTN FILES BY DAY
#
#
# Y.Shen, 02-28-2011
# Haiying Gao, 10-11-2016
# Cong Li, 10-24-2016
########################################################################
#          SET WORK DIRECTORY
###########################################
set area = Alaska
set wkdir = /Users/xiaotaoyang/Work/FWANT/$area
set sacdir = $wkdir/sac_FTN
set dailyFTNdir = $wkdir/daily_FTN
if ( ! ( -e $dailyFTNdir ) ) mkdir $dailyFTNdir
cd $sacdir
foreach network ("`awk '{if ( NR >0 ) print}' $wkdir/alaskanetcodes.txt`") #use cat instead of awk, this is a temporary change.
echo $network
set stnlst = $wkdir/data_reqs/$network".txt"
echo $stnlst
foreach station ( `cat $stnlst | awk '{print $1"."$2}'` )
echo $station
if ( -d $station ) then #only cd to existing stations.
	cd $station
	foreach sacfile ( `/bin/ls ftn.*.SAC` )
	set day = `echo $sacfile | awk '{print substr($1,5,8)}'`
	set dailydir = $dailyFTNdir/$day
	if ( ! (-e $dailydir) ) mkdir $dailydir
	
	echo $sacfile $dailydir
	if (! (-e $dailydir/$sacfile) ) then
		cp $sacfile $dailydir
	else
		echo "***" $sacfile exists in $dailydir
	endif
	end #foreach sacfile
	cd ..
endif

end #foreach station
end #foreach network
