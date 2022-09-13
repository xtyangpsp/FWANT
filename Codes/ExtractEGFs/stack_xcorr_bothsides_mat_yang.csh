#!/bin/csh
####################################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO STACK CROSS-CORRELATIONS, CALCULATE
# EGFs AND UNCERTAINTIES
#
#
# last modified,
# Y.Shen, 02-28-2011
# Haiying Gao, 10-11-2016
# Cong Li, 8/10/2017
#       Change the structure, call matlab once
######################################################################################
matlab -nodesktop -nojvm -nosplash <<EOF
addpath('/Users/congli/FWT_Research/EPM/Virginia/matlab_script/');
clc;clear;close all;
area='Virginia';
wkdir = ['/Users/congli/FWT_Research/EPM/' area];
xcorrdir = [wkdir '/xcorr'];
stackdir = [wkdir '/stacked_xcorr'];
if (exist(stackdir)~=7)
mkdir(stackdir);
end
outputdir=stackdir;
network = textread('rdseed.pl','%s');
cd(xcorrdir);
for i=1:length(network)
a=strcat(network(i),'*');
stndir=dir(char(a));
for j=1:length(stndir)
stnpair{j}=stndir(j).name;
stnpair{j}
cd(stnpair{j});
if (exist('filelist.txt','file')==2)
!/bin/rm filelist.txt
end
!/bin/ls -f xc.*  > filelist.txt
files=textread('filelist.txt','%s');
nfiles=length(files);
nmonths=floor(nfiles/30); %daily??
nptmax=16001;%??
nmonthmax=200; %??
ymonth=zeros(nptmax,nmonthmax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Daily Stack
for k=1:nfiles
imonth=floor(k/30)+1;imonthres=k-(imonth-1)*30; %daily
clear file data xx yy hdr delta
file=char(files(k));
data=rsac(file);
xxtmp=data(:,1); % lag time
yytmp=data(:,2); % xcorrelation
npts=length(xxtmp);
nptshalf=(npts-1)/2+1;
xx = xxtmp(nptshalf-(nptmax-1)/2:nptshalf-(nptmax-1)/2+nptmax-1);
yy = yytmp(nptshalf-(nptmax-1)/2:nptshalf-(nptmax-1)/2+nptmax-1);
yymax=max(abs(yy));
yy=yy/yymax; %Normalization
hdr=data(:,3); %Header
delta=hdr(1);   %dt
if(k == 1) %stack
yyy=yy;
else
yyy=yyy+yy;
end
if (imonthres == 1)
ymonth(:,imonth)=yy;
else
% if imonth==27
%     pause;
% end
ymonth(:,imonth)=ymonth(:,imonth)+yy; % stacking daily xcorr into monthly
end
end % file (files stack)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Extract the EGFs
if npts>nptmax,
npts=nptmax;
end
nptshalf=(npts-1)/2+1;
tpos=xx;dt=xx(2)-xx(1);
sig=yyy; % the stack of all daily xcorr files
egf=zeros(npts,1);
for k=2:nptshalf-1
egf(k)=-(sig(k-1)-sig(k+1))/(2.*dt); % later signal (from the center) subtract the earlier signal
end
for k=nptshalf+1:npts-1
egf(k)=-(sig(k+1)-sig(k-1))/(2.*dt); % derivative of Xcorr, get the EGFs
end

% dress the egf
egf=egf-mean(egf); %demean
w=tukeywin(length(egf),0.02);
egf=egf.*w;  %taper (2000*0.01/2=5 point at each end)

egf=egf/imonth;  %since egf err is obtained from monthly std ??

egfmonth=zeros(npts,imonth);

for im=1:imonth
% monthly stack
for k=2:nptshalf-1
egfmonth(k,im)=-(ymonth(k-1,im)-ymonth(k+1,im))/(2.*dt);
end
for k=nptshalf+1:npts-1
egfmonth(k,im)=-(ymonth(k+1,im)-ymonth(k-1,im))/(2.*dt);
end
egfmonth(:,im)=egfmonth(:,im)-mean(egfmonth(:,im)); %demean but no taper
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard error of the mean
stderr=zeros(npts,1);
for k=1:npts
stderr(k)=std(egfmonth(k,1:imonth))/sqrt(imonth); % standard deviation
end
%save the egfs and uncertainty
egfdir = [outputdir '/' char(stnpair{j})];
dirflag = exist(egfdir,'dir');
if dirflag ~= 7
mkdir(egfdir);
end
egfname = [egfdir '/egf.' char(stnpair{j}) '.ZZ.sym.SAC'];

output=zeros(npts,3);
output(:,1)=tpos;% Lagtime
output(1:306,3)=data(1:306,3); % see cal_xcorr.csh
output(6,3)=tpos(1); %B
output(7,3)=tpos(end); %E
output(80,3)=length(tpos); % npts
output(71,3)=2599; %nzyear
output(72,3)=1; %nzjday
output(73,3)=0; %nzhour
output(74,3)=0; %nzmin
output(75,3)=0; %nzsec
output(76,3)=0; %nzmsec
output(:,2)=egf;
wsac(egfname,output);pause(0.1);

egferrname = [egfdir '/egf.' char(stnpair{j}) '.ZZ.sym.err.SAC'];
output(:,2)=stderr;
wsac(egferrname,output);pause(0.1);

for im=1:imonth
if ( im < 10 ); monchar = ['00' num2str(im)];end;
if ( im >= 10 & im < 100 ); monchar = ['0' num2str(im)];end;
if ( im >= 100 ); monchar = [ num2str(im)];end

egfmonthname = [egfdir '/egf.' char(stnpair{j}) '.ZZ.sym.m' monchar '.SAC'];

output(:,2)=egfmonth(:,im);
wsac(egfmonthname,output);pause(0.1);
end
!/bin/rm filelist.txt
cd ..
end % station pair
end %network
EOF
