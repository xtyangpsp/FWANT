#!/bin/csh
#######################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO REMOVE TIME SEGMENTS CONTAINING
# LARGE EARTHQUAKES
#
# last modified,
# Y.Shen, 02-28-2011
# Haiying Gao, 10-11-2016
# Cong Li, 8/6/2017
########################################################################
#          SET WORK DIRECTORY
###########################################
set area = Alaska
set wkdir = /Users/xiaotaoyang/Work/FWANT/$area

set sacdir = $wkdir/daily_FTN
cd $sacdir

set delEQdir = $wkdir/daily_FTN_delEQ
if ( ! ( -e $delEQdir ) ) mkdir $delEQdir
set stnlst = stn4day.lst    # tmp lst for stations with day directories
set saclst = sacfile.lst    # tmp lst for sacfile with the day directories

# EQ list in the following format
set eqlst = $wkdir/EQinfo/NEIC_2007_201709m5.5_formatted.txt
cat /dev/null > $delEQdir/log.txt
foreach day ( `/bin/ls -d 2*.*` )
# pwd
cd $day
echo $day

cat /dev/null > $stnlst
cat /dev/null > $saclst
foreach sacfile ( `/bin/ls -f ftn.*.SAC ` ) # .7D.J09.* UW.EHZ 2005-2011
set stnm = `echo $sacfile | awk -F. '{print $4"."$5}'`
echo $stnm >> $stnlst
echo $sacfile >> $saclst
# echo 
end
set year = `echo $day | awk '{print substr($1,1,4)}'`
set jday = `echo $day | awk '{print substr($1,6,3)}'`
echo $year, $jday
set targetdir = $delEQdir/$day
if ( !( -e $targetdir ) ) mkdir $targetdir

set nsac = `ls -1 ftn*.*.SAC | tail -n 1 | wc -l`
set nsac_target = `ls -1 $targetdir/ftn*.*.SAC | tail -n 1 | wc -l`

if ( $nsac != ${nsac_target} ) then  #only start the work on new days, skip existing files.
echo 'Working on day: ' $day ' ...'


calday $jday $year > temp
set mon = `cat temp | grep Calendar | awk '{print $3}'`
set dayofmon = `cat temp | grep Calendar | awk '{print $4}'`
set yearmonday = $year"-"$mon"-"$dayofmon
set eqlstday = tempEQlst
/bin/cat $eqlst | grep $yearmonday > $eqlstday # find the date of events in event catalog

set nday = $jday
set ch011 = `echo $jday | awk '{print substr($1,1,1)}'`
if ( $ch011 == 0 ) then
set nday = `echo $jday | awk '{print substr($1,2,2)}'`
endif
set ch001 = `echo $jday | awk '{print substr($1,1,2)}'`
if ( $ch001 == 00 ) then
set nday = `echo $jday | awk '{print substr($1,3,1)}'`
endif

@ jday1 = $nday - 1
calday $jday1 $year > temp
set mon = `cat temp | grep Calendar | awk '{print $3}'`
set dayofmon = `cat temp | grep Calendar | awk '{print $4}'`
set yearmonday1 = $year"-"$mon"-"$dayofmon
/bin/cat $eqlst | grep $yearmonday1 >> $eqlstday

@ neq = 0
foreach eq ( `cat $eqlstday | grep $yearmonday | awk '{print $1}'` )
@ neq = $neq + 1
if neq == 1 then
	echo $jday > jday.dat
else
	echo $jday >> jday.dat
endif
end

foreach eq ( `cat $eqlstday | grep $yearmonday1 | awk '{print $1}'` )
@ neq = $neq + 1
if neq == 1 then
        echo $jday1 > jday.dat
else
        echo $jday1 >> jday.dat
endif
end

echo "neq = " $neq



if ( $neq == 0 ) then
/bin/cp ftn*.*.SAC  $targetdir
/bin/rm $eqlstday jday.dat temp $stnlst $saclst

else
/bin/cat $eqlstday | awk '{print $3}' > elat.dat
/bin/cat $eqlstday | awk '{print $4}' > elon.dat
/bin/cat $eqlstday | awk '{print substr($2,1,2)}' > hr.dat
/bin/cat $eqlstday | awk '{print substr($2,4,2)}' > min.dat
/bin/cat $eqlstday | awk '{print substr($2,7,2)}' > sec.dat

#%%%% Check the portion that highlighted %%%%%%%%

matlab -nodesktop -nojvm -nosplash <<EOF

addpath('../../matlab_script/');

!echo $targetdir > targetdir
targetname = textread('targetdir','%s');

!echo $day > daynumber
dayname = textread('daynumber','%s');
tmp=char(dayname);
tday=str2num(tmp(6:8)); tmp=[];

!/bin/cp $stnlst stnames
stn=textread('stnames', '%s');
nst = length(stn);

!/bin/cp $saclst saclst
sacfiles=textread('saclst','%s');

load elat.dat; load elon.dat; load hr.dat; load min.dat; load sec.dat; load jday.dat
neq=length(elat);
origin = (jday-tday)*24*3600 + hr*3600 + min*60 + sec;

for m=1:nst
stname=stn(m);
sacfile = char(sacfiles(m));
f1=fopen('/Users/xiaotaoyang/Work/FWANT/Alaska/daily_FTN_delEQ/log.txt','a+');
fprintf(f1,'%s\n',sacfile);
fclose(f1);
iflag = exist(sacfile, 'file');
if iflag == 2
	data=rsac(sacfile);
	tt=data(:,1); dt=tt(2)-tt(1);nt=length(tt);
	trace=data(:,2);
	hdr=data(1:306,3);
	slat=data(32,3);slon=data(33,3);
	w=ones(nt,1);
	for k=1:neq
	dist=geo2dist(elat(k),elon(k),slat,slon);
	tmin=origin(k)+dist/14.5; % the estimated minimum arrival
	tmax=origin(k)+dist/1.5; % the estimated maximum arrival

	itmin=round(tmin/dt);
	itmax=round(tmax/dt); 
	if itmin > nt;itmin=nt;end
	if itmax > nt;itmax=nt;end
	if itmin < 1; itmin=1;end
	if itmax < 1; itmax=1;end
	N=itmax-itmin+1; % the estimated length of removal
	tukeywindow = tukeywin(N,0.3);
    w(itmin:itmax)=tukeywindow(1:N)'-1.0; %construction of filter
	end

	trace=trace.*w;

	xc(:,1)=data(:,1);
	xc(:,2)=trace;
	xc(1:306,3)=data(1:306,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	xcname = [char(targetname) '/' sacfile];
    wsac(xcname,xc);

end % iflag == 2
end % m

EOF

/bin/rm daynumber *.dat stnames temp tempEQlst targetdir saclst
/bin/rm $stnlst $saclst

endif # if neq == 0

cd ..

else
	echo $day 'data exists'
	cd ..
# 	exit 0
endif


end #foreach day


