#!/bin/csh
############################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO COUNT THE NUMBER OF SAC FILES
# IN ./sac_data
#
#  Note: change the nfile limit $min_nfile
#
# Y.Shen, 03-09-2011
#############################################################################
# set area = Virginia
set wkdir = '/Users/xiaotaoyang/Work/FWANT/Alaska'
set sacdir = $wkdir/sac_data
#set tempdir = $sacdir/temp
# stations to be used should have more than the minimum number of sac files
set min_nfile = 90
cd $sacdir
if ((-e ../tempdel)) /bin/rm -r ../tempdel
@ nstn = 0
foreach network ( "`cat /Users/xiaotaoyang/Work/FWANT/Alaska/alaskanetcodes.txt`" )
echo $network
#set network = 7D
#set stnlst = ~/FWT/DataRequest/finergrid/SeattleBasinStations/$network"_EN".txt
#set stnlst = ~/FWT/DataRequest/finergrid/SeattleBasinStations/StationList_sb.txt
#foreach stn (`cat $stnlst | awk '{print $2"."$1}'`)
#foreach stn (CN.HOPB CN.MGB CN.PFB CN.SNB CN.VGZ CN.YOUB)
foreach stn ( `/bin/ls -d $network*` )
#echo $stn
#foreach comp (LHZ)
if ( !(-e ../tempdel)) mkdir ../tempdel
cd $stn
@ nstn = $nstn + 1

@ nfile = 0
#foreach sacfile ( `/bin/ls -f *.$comp.SAC` )
foreach sacfile ( `/bin/ls -f *.SAC` )
@ nfile = $nfile + 1
end

cd ..
if ( $nfile <= $min_nfile ) then
echo  " --------- less than " $min_nfile " ------- " $stn $nfile
#/bin/rm -r $stn
mv $stn ../tempdel
endif

#end
end

echo $nstn
end
