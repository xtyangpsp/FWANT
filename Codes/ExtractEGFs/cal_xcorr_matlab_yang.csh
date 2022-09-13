#!/bin/csh
##############################################################################################
# MODIFIED FROM YANG SHEN AND HAIYING'S CODE TO CALCULATE CROSS-CORRELATION BETWEEN PAIRS
# OF SEISMIC RECORDS
#
#  Y.Shen, 02-28-2011
#  Haiying Gao, 10-11-2016
#  Cong Li, 08/07/2017
#       Parallel computing
#       Checking exist files
##############################################################################################
#    SET WORK DIRECTORY
##############################################################################################
matlab -nodesktop -nosplash <<EOF
clc;clear;
addpath('/Users/xiaotaoyang/Work/FWANT/matlab_script/');
area='Alaska';
wkdir = ['/Users/xiaotaoyang/Work/FWANT/' area];
sacdir = [wkdir '/daily_FTN_delEQ'];
cd (sacdir);
xcorrdir = [wkdir '/xcorr']
outputdir = xcorrdir;
if (exist(xcorrdir)~=7)
mkdir(xcorrdir)
end
if (exist('dayfolder.txt','file')==2)
!/bin/rm dayfolder.txt
end
!/bin/ls -d 2* > dayfolder.txt
day=textread('dayfolder.txt','%s');
day_num=length(day);
f1=fopen('log_xcorr.txt','a+');
sacfiles_exist=textread('log_xcorr.txt','%s');
for i=1:day_num
cd (day{i});
if (exist('filelist_exist.txt','file')==2)
!/bin/rm filelist_exist.txt
end
!/bin/ls -f *.SAC > filelist_exist.txt
sacfiles=textread('filelist_exist.txt','%s');
Lia = unique(ismember(sacfiles,sacfiles_exist));
if (length(Lia)==0 | find(Lia==0))
[ntst,index]=regexp(sacfiles,'\w*(?s)\w*','match');
stn={};
for j=1:length(ntst)
stn{j}=[char(ntst{j}(4)) '.' char(ntst{j}(5))];
end
nst = length(stn);

npts=172801; %86400; % daily sac with 1 sps
w=tukeywin(npts,0.01);
sacdata=zeros(npts,3,nst);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:nst
sacfilename=char(sacfiles(m));
sacdata(:,:,m)=rsac(sacfilename);
end
tt=sacdata(:,1,1); dt=tt(2)-tt(1);

% maximum time lag ( using maximum travel times for the area )
maxlag = round(400/dt);
CoreNum=4;
parpool('local',CoreNum);
parfor m=1:nst-1
stname1=stn(m);

data1=sacdata(:,:,m);
trace1=data1(:,2).*w; % taper

for m2 = m+1:nst
stname2 = stn(m2);

data2=sacdata(:,:,m2);
trace2=data2(:,2).*w; % taper

% For two known time functions A and B and A signal arrives earlier, tests show
% xcorr(A,B) yields the peak value at negative leg, xcorr(B,A) yields the peak at positive leg.
% thus for the two stations, trace2 should be the first var to match the convention in seismology.
%	[c,lags] = xcorr(trace1,trace2,maxlag,'unbiased');
[c,lags] = xcorr(trace2,trace1,maxlag,'unbiased'); %
xc=zeros(length(lags),3);
xc(:,1)=lags'*dt;
xc(:,2)=c;
xc(1:306,3)=data1(1:306,3);%header
xc(6,3)=lags(1)*dt; %B see read_sac.m in ~/matlab
xc(7,3)=lags(length(lags))*dt; %E
xc(80,3)=length(c); %NPTS
xc(32,3)=data2(32,3);xc(33,3)=data2(33,3);xc(34,3)=data2(34,3); %stla,stlo,stel

xc(36,3)=data1(32,3);xc(37,3)=data1(33,3);xc(38,3)=data1(34,3); %evla,evlo,evel

%xc(111,3)=data1(111,3);char;xc(112,3)=data1(112,3);xc(113,3)=data1(113,3);xc(114,3)=data1(114,3);
xc(115,3)=data2(111,3);xc(116,3)=data2(112,3);xc(117,3)=data2(113,3);xc(118,3)=data2(114,3); %kstnm, 4 characters from each stn, see rsac.m

xcdir = [outputdir '/' char(stname1) '.' char(stname2)];
xcflag = exist(xcdir,'dir');
if xcflag ~= 7
mkdir(xcdir);
end
xcname = [xcdir '/xc.' char(day{i}) '.' char(stname1) '.' char(stname2) '.ZZ.SAC'];
wsac(xcname,xc); %pause(0.01);
end % m2

end % m
for k=1:length(sacfiles)
fprintf(f1,'%s\n',sacfiles{k});
end % writting the log
!/bin/rm filelist_exist.txt
else
!/bin/rm filelist_exist.txt
end
clear sacfiles sacdata data1 data2 trace1 trace2 tt w 
cd ..
end
fclose(f1);
p = gcp;
delete(p);
EOF
