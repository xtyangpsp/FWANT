#!/bin/bash
###########################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO EXTRACT DAILY SAC FILES FROM A SEED VOLUME AND INTERPOLATE SAC
# FILES TO A FIXED SPS
#
#  Usage:
#
#  Y.Shen, 03-09-2011
#  Haiying,10-01-2016
#  Cong,   05-01-2017  Add : check out the reduplicated files
# Xiaotao Yang, 9/20/2017: changed to read seed files named by net.sta.*
###########################################################
#   SET AND MAKE DIRECTORY
###########################################################
#set wkdir1 = /Volumes/DataForFWT/Cascadia
#wkdir='/Users/xtyang/Work/Research/Postdoctoral/UMass_PD/Alaska'
wkdir='/Users/xiaotaoyang/Work/FWANT/Alaska'
#set network = $1

netstachanfile='netstachan_afterCN.LIB.txt' #format: NET.STA.CHAN
sacdir=$wkdir'/sac_data'
seedfolder=$wkdir'/RequestedData/stationseedfiles'
respdir=$wkdir'/inst_resp' # Instrument response directory

if [ ! -d $sacdir ]; then mkdir $sacdir; fi
if [ ! -d $respdir ]; then mkdir $respdir;fi
if [ ! -d $seedfolder ]; then echo 'seed dir [ '$seedfoler' ]not found!'; exit 1;fi
cd $seedfolder
/bin/rm -f *.SAC # not sure this is necessary here. keep this from the original script.

for netstachan in `cat $wkdir/$netstachanfile`
do
	#clear existing response files and output log file.
	/bin/rm -f  SAC_PZs*
	/bin/rm -f  rdseed.output
	
	echo $netstachan
	network=`echo $netstachan | sed 's/\./ /g' | awk '{print $1}'`
	stn=`echo $netstachan | sed 's/\./ /g' | awk '{print $2}'`
	comp=`echo $netstachan | sed 's/\./ /g' | awk '{print $3}'`
	
	station=$network"."$stn
	
	stndir=$sacdir'/'$station #mkdir for net.sta sac files.
	if [ ! -d $stndir ];then mkdir $stndir; fi
	nseed=0;
	nseed=`ls -1 *$station'.'*.seed | wc -l`
	if [ $nseed -gt 1 ]
	then
		echo 'more than one '$station' file found! Please clean the seed directory: '$seedfolder
		exit 1
	fi
	
	if [ $nseed -eq 0 ]
	then
		echo 'no seed file for '$station' ! Please check the station file: '$wkdir/$netstachanfile
		echo 'no seed file for '$station $comp' ! Please check the station     file: '$wkdir/$netstachanfile \
		|  mail -s 'rdseed2sac stoped' xiaotaoyang@geo.umass.edu
		exit 1
	fi
	
	seedfile=`ls -f -1 *$station'.'*.seed | awk '{if ( NR==1 ) print}'`
	echo $seedfile
	
	##################################################################
	#  EXTRACT DAILY DATA FROM SEED FILE WITH RDSEED SORFTWARE, AND CONVERT THEM INTO
	#  SAC FORMAT. THE COMMAND OF RDSEED CAN BE REFERENCED IN
	#  http://ds.iris.edu/ds/nodes/dmc/manuals/rdseed/#user-prompt-mode OR TYPE rdseed -h or -u
	##################################################################	
	#2007 2008 2009 2011 2012 2013 2014 2015
	#check data time range, only use start and end years and days to save loop time.
# 	startyear=`rdseed -tvf 1 $seedfile | awk '{ if ( NR ==7 ) print}' | \
# 		sed 's/,/ /g' | awk '{print $6}'`
# 	startday=`rdseed -tvf 1 $seedfile | awk '{ if ( NR ==7 ) print}' | \
# 		sed 's/,/ /g' | awk '{print $7}'`
# 	
# 	#the above fails when there is a local ID colomn, do the following then (shit one column to the right.
# # 	echo $startyear
# # exit 0
# 	if [ $startyear == 'E' ] || [ $startyear == 'D' ] || [ $startyear == 'M' ] \
# 			|| [ $startyear == 'Q' ] || [ $startyear == 'R' ]
# 	then
# 		startyear=`rdseed -tvf 1 $seedfile | awk '{ if ( NR ==7 ) print}' | \
# 	                 sed 's/,/ /g' | awk '{print $7}'`
#         startday=`rdseed -tvf 1 $seedfile | awk '{ if ( NR ==7 ) print}' | \
#         	         sed 's/,/ /g' | awk '{print $8}'`
# 	fi
	((startyear=2007))
	((startday=1))
	endyear=`rdseed -tvf 1 $seedfile | tail -n 1 | sed 's/,/ /g' | awk '{print $6}'`
	endday=`rdseed -tvf 1 $seedfile | tail -n 1 | sed 's/,/ /g' | awk '{print $7}'`
	#the above fails when there is a local ID colomn, do the following then (shit one column to the right.
	if [ $endyear == 'E' ] || [ $endyear == 'D' ] || [ $endyear == 'M' ] \
			|| [ $endyear == 'Q' ] || [ $endyear == 'R' ]
	then
		endyear=`rdseed -tvf 1 $seedfile | tail -n 1 | sed 's/,/ /g' | awk '{print $7}'`
		endday=`rdseed -tvf 1 $seedfile | tail -n 1 | sed 's/,/ /g' | awk '{print $8}'`
	fi
	
	((year=$startyear))
	while [ $year -le $endyear ] 
	#for year in 2007 2008 2009 2011 2012 2013 2014 2015 2016 2017 #good for Alaska data
	do
		#the following way to check if the year is leap year
		((x4=$year%4));((x100=$year%100));((x400=$year%400))
		if [ $x4 -eq 0 ] && ( [ $x100 -gt 0 ] || [ $x400 -eq 0 ] )
# 		if [ $year == 2008 ] || [ $year == 2004 ]|| [ $year == 2000 ]||[ $year == 2012 ]||[ $year == 2016 ]
		then
			jdaymax=366
		else
			jdaymax=365
		fi

		hour=00
		min=00
		sec=00.000
		time1=$hour":"$min":"$sec
		time2=$hour":"$min":"30.000 #?? to avoid taper edge effects.
		##echo $jdaymax $time1 $time2
		
		if [ $year -eq $startyear ]
		then 
			((dd1=$startday))
		else
			((dd1=1))
		fi
		if [ $year -eq $endyear ]
		then
			((jdaymax=$endday))
		fi
		((dd2=$dd1 + 1))
		while [ $dd1 -le $jdaymax ]
		do
			##echo $dd1 $dd2
			jday1=$dd1
			jday2=$dd2
			if [ $dd1 -lt 10 ]; then jday1="00"$dd1; fi
			if [ $dd1 -ge 10 ] && [ $dd1 -lt 100 ]; then jday1="0"$dd1; fi
			if [ $dd2 -lt 10 ]; then jday2="00"$dd2; fi
			if [ $dd2 -ge 10 ] && [ $dd2 -lt 100 ]; then jday2="0"$dd2; fi

			##echo $jday1 $jday2

			if [ $dd1 == $jdaymax ]; then jday2=$jday1; time2="23:59:59.999"; fi
			targetsacname=$year.$jday1.$station.$comp.SAC
			if [ ! -e $stndir/$targetsacname ]
			then				
#
# parameter of rdseed for user prompt mode
# 1----  Input  File (/dev/nrst0) or 'Quit' to Exit
# 2----  Output File (stdout)
# 3----  Volume #  [(1)-N] Selection of the number of volume; for seed file, this parameter is always 1
# 4----  Options [acCsSpRtde] -d or -e
# 5----  Summary file (None)
# 6----  Station List (ALL)
# 7----  Channel List (ALL)
# 8----  Network List (ALL)
# 9----  Loc Ids (ALL ["--" for spaces])
# 10----  Output Format [(1=SAC), 2=AH, 3=CSS, 4=mini seed, 5=seed, 6=SAC ASCII, 7=SEGY, 8=Simple ASCII(SLIST), 9=Simple ASCII(TSPAIR)]
# 11---  Output file names include endtime? [Y/(N)]
# 12---  Output poles & zeroes ? [Y/(N)]
# 13---  Check Reversal [(0=No), 1=Dip.Azimuth, 2=Gain, 3=Both]
# 14---  Select Data Type [(E=Everything), D=Data of Undetermined State, M=Merged data, R=Raw waveform Data, Q=QC'd data]
# 15---  Start Time(s) YYYY,DDD,HH:MM:SS.FFFF
# 16---  End Time(s)   YYYY,DDD,HH:MM:SS.FFFF
# 17---  Sample Buffer Length [20000000]
# 18---  Extract Responses [Y/(N)]
# 19---  Output data format will be sac.binary.	
# echo $seedfile	
echo $station  $comp $year,$jday1,$time1
# echo $year,$jday2,$time2
rdseed <<!>> rdseed.output
$seedfile


d

$stn
$comp


1
N
Y
3
E
$year,$jday1,$time1
$year,$jday2,$time2
20000000
N
quit
!

# exit 0
#end of calling rdseed
					
###########################################################################
# INTERPOLATE UNIFORM SAC-FORMAT DATA TO A FIXED SPS, REMOVE INSTRUMENT RESP.
# REMOVE THE MEAN VALUE, TREND AND FILTER
#
###########################################################################

# interpolate first improves run time, since "transfer" appears to be most time consuming
# inspection of the waveforms shows no difference
#
# The format of sac file converted from seed using rdseed is
# yyyy.dd.hh.mm.ss.ffff.nn.ssss.LL.CCC.Q.SAC
# yyyy is year
# ddd  is day
# hh.mm.ss is hour,minute and second
# ffff is milisecond
# nn is network
# SSSS is station
# LL is the location ID (http://ds.iris.edu/ds/newsletter/vol1/no1/specification-of-seismograms-the-location-identifier/)
# CCC is component
# Q is quality control or local ID
			nsac=0
			nsac=`ls -f -1 *.$station.*.SAC | wc -l`
 			echo $nsac
			
			if [ $nsac -eq 0 ]
			then 
				((dd1=$dd1 + 1))
		                ((dd2=$dd1 + 1))
				continue
			fi 
			if [ $nsac -eq 1 ]
			then #case for only have one local ID(instrument)
				sacfile1=$sacfile
			fi
			if [ $nsac -ge 1 ]
			then
			#remove "10" if there is "00" stream;case for multiple local ID
				streamflag=0
				
				for sacfile in `ls -f -1 *.$station.*.SAC`
				do
# 					echo $sacfile
					stream=`echo $sacfile | awk -F. '{print $9}'`
					# if use ls, should change to $9 instead of $10
					if [ "'$stream'" == "00" ]; then streamflag=1; fi
				done #end loop for local IDs
				
				if [ $streamflag -eq 1 ];
				then
					/bin/rm *.$station.10.$comp.*.SAC 'SAC_PZs_'*'_'$comp'_10'* 
				# ignore the secondary stream
				fi
				
				(( nsac2=0 ))
				# read the data stream 00 and other (but not 10)
				for sacfile in `ls -f *.$station'.'$stream'.'$comp.*'.SAC' ` 
				do
					((nsac2=$nsac2 + 1))
					if [ $nsac2 == 1 ]; then sacfile1=$sacfile; fi
					if [ $nsac2 == 2 ]; then sacfile2=$sacfile; fi
					if [ $nsac2 == 3 ]; then sacfile3=$sacfile; fi
					if [ $nsac2 == 4 ]; then sacfile4=$sacfile; fi
					if [ $nsac2 == 5 ]; then sacfile5=$sacfile; fi
					if [ $nsac2 == 6 ]; then sacfile6=$sacfile; fi
					if [ $nsac2 == 7 ]; then sacfile7=$sacfile; fi
					if [ $nsac2 == 8 ]; then sacfile8=$sacfile; fi
					if [ $nsac2 == 9 ]; then break; fi #this basically will cause skipping days with > 8 fregments.
				done
				
				##this basically will cause skipping days with > 8 fregments.
				if [ $nsac2 -gt 8 ]
                        	then 
                                	((dd1=$dd1 + 1))
                                 	((dd2=$dd1 + 1))
					/bin/rm *.$station'.'$stream'.'$comp.*.SAC
                                 	/bin/rm SAC_PZs*
					continue
                         	fi
# 				echo $stream
				# construct sac macro
# 				echo $network
				respid="SAC_PZs_"$network"_"$stn"_"$comp"_"$stream
				respfile=`ls -f $respid*`
# 				echo '*.'$station'.'$stream'.'$comp'.*.SAC'
 				#echo $nsac2
				echo "echo off errors warnings output commands macros processed"> interpsac.m
				if [ $nsac2 == 1 ]; then
					echo "r $sacfile1" >> interpsac.m
				fi
				if [ $nsac2 == 2 ]; then
					echo "r $sacfile1" >> interpsac.m
					echo "merge $sacfile2" >> interpsac.m
				fi
				if [ $nsac2 == 3 ]; then
					echo "r $sacfile1" >> interpsac.m
					echo "merge $sacfile2" >> interpsac.m
					echo "merge $sacfile3" >> interpsac.m
				fi
				if [ $nsac2 == 4 ]; then
					echo "r $sacfile1" >> interpsac.m
					echo "merge $sacfile2" >> interpsac.m
					echo "merge $sacfile3" >> interpsac.m
					echo "merge $sacfile4" >> interpsac.m
				fi
				if [ $nsac2 == 5 ]; then
					echo "r $sacfile1" >> interpsac.m
					echo "merge $sacfile2" >> interpsac.m
					echo "merge $sacfile3" >> interpsac.m
					echo "merge $sacfile4" >> interpsac.m
					echo "merge $sacfile5" >> interpsac.m
				fi
				if [ $nsac2 == 6 ]; then
					echo "r $sacfile1" >> interpsac.m
					echo "merge $sacfile2" >> interpsac.m
					echo "merge $sacfile3" >> interpsac.m
					echo "merge $sacfile4" >> interpsac.m
					echo "merge $sacfile5" >> interpsac.m
					echo "merge $sacfile6" >> interpsac.m
				fi
				if [ $nsac2 == 7 ]; then
					echo "r $sacfile1" >> interpsac.m
					echo "merge $sacfile2" >> interpsac.m
					echo "merge $sacfile3" >> interpsac.m
					echo "merge $sacfile4" >> interpsac.m
					echo "merge $sacfile5" >> interpsac.m
					echo "merge $sacfile6" >> interpsac.m
					echo "merge $sacfile7" >> interpsac.m
				fi
				if [ $nsac2 == 8 ]; then
					echo "r $sacfile1" >> interpsac.m
					echo "merge $sacfile2" >> interpsac.m
					echo "merge $sacfile3" >> interpsac.m
					echo "merge $sacfile4" >> interpsac.m
					echo "merge $sacfile5" >> interpsac.m
					echo "merge $sacfile6" >> interpsac.m
					echo "merge $sacfile7" >> interpsac.m
					echo "merge $sacfile8" >> interpsac.m
				fi
				
				#echo "echo on" >> interpsac.m
				# taper 600 s at each end, taperwid = 600s /86400 s in a day = 0.00694
				taperwid=0.01
				echo "rmean;rtr;taper w $taperwid" >> interpsac.m
				#echo "lp bu co 20 n 4 p 2" >> interpsac.m
				echo "lp bu co 0.4 n 4 p 2" >> interpsac.m

				#echo "interpolate d 0.02 b 0.02 n 4320000" >> interpsac.m
				# LHZ is sampled at 1 sps
				#echo "interpolate d 1 b 1 n 86400" >> interpsac.m
				echo "interpolate d 0.5 b 0" >> interpsac.m
				#echo "transfer from polezero subtype $respfile to none freq 0.05 0.0625 20 25" >> interpsac.m
				# down to "hum" frequency
				echo "transfer from polezero subtype $respfile to none freq 0.001 0.00125  0.5 1.0">>interpsac.m
				echo "rmean;rtr;taper w $taperwid" >> interpsac.m
				echo "w $targetsacname" >> interpsac.m
				echo "quit" >> interpsac.m

				echo '    --Done: '$targetsacname
#run sac with the macro file.
sac <<EOF>> sac.output
m interpsac.m
quit
EOF
				
				/bin/mv $targetsacname $stndir'/'
 				#exit 0
				/bin/rm *.$station.*.SAC
				# for ff in `find . -name '*.'$station'.*.SAC' -type f -print`
# 				do
# 					/bin/rm $ff #$stndir'/'
# 				done
				#/bin/rm -f *.$station.*.SAC
				#/bin/rm *.$station.*.SAC
				/bin/mv SAC_PZs* $respdir
				/bin/rm sac.output
			fi # nsac >= 1
			/bin/rm rdseed.output
		else
		   echo $stndir/$targetsacname ' exists'
		fi
		((dd1=$dd1 + 1))
		((dd2=$dd1 + 1))
		done # loop for dd1
	((year=$year + 1))
	done #end loop for years
done #loop for netstachan
/bin/rm rdseed.err_lo* interpsac.m
