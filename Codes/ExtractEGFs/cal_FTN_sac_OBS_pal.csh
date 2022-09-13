#!/bin/sh

#  cal_FTN_sac_OBS_pal.csh
#  
#
#  Created by Cong Li on 4/6/17.
#  5/5/2017 Checking existed files
matlab -nodesktop -nosplash <<EOF
clc;clear;
addpath('/Users/xiaotaoyang/Work/FWANT/matlab_script/');
area='Virginia';
wkdir = ['/Users/xiaotaoyang/Work/FWANT/Alaska'];
sacdir = [wkdir '/sac_data'];
cd (sacdir);
if (exist('tempdel')==7)
rmdir('tempdel','s');
end
sacFTNdir = [wkdir '/sac_FTN'];
if (exist(sacFTNdir)~=7)
mkdir(sacFTNdir)
end
!/bin/ls -d 5C.HRLQ > listfolder.txt %test purpose
%!/bin/ls -d * > listfolder.txt
list=textread('listfolder.txt','%s');
list_num=length(list);
for i=1:list_num
[ntst,index]=regexp(list{i},'\w*(?s)\w*','match');
ntk=char(ntst{1});
stn=char(ntst{2});
station=char(list{i});
stndir = [sacFTNdir '/' station];
if (exist(stndir)~=7)
mkdir(stndir);
end
cd(stndir);
E_files={}; Efiles={};
!/bin/ls -f *.SAC > filelist_exist.txt
Efiles=textread('filelist_exist.txt','%s');
Enfiles = length(Efiles);
for et=1:Enfiles
E_files{et}=Efiles{et}(5:end);
end
E_files=E_files';
!/bin/rm -f filelist_exist.txt
cd(sacdir);
cd(station);
station
fl=0.0033; % default low frequency end
targetdir=stndir;
%Initialize Matlab Parallel Computing Enviornment
outputstn = targetdir;
files={}; S_files={};
!/bin/ls -f *.SAC > filelist.txt
S_files=textread('filelist.txt','%s');
files=setdiff(S_files,E_files);
nfiles = length(files);
if nfiles==0
cd ..
continue;
end
targetdirname=targetdir;
% define frequency range to be studied
df=fl/4; % sampling interval in Frequency
fh=0.4; % high frequencyend T(min)~2s
nf=round((fh-fl)/df)+1;% npt in Frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define filter and parameters
N=2; b=zeros(nf,5);a=zeros(nf,5);
%%% dt0 in second must be consistent with the sampling of the sac files
dt0=0.5; sps=1/dt0;
%
b=[]; a=[];
Feff=1/dt0/2;
for kf=1:nf
fone=(fl+(kf-1)*df)/Feff;
ftwo=(fl+kf*df)/Feff;
[b(kf,:),a(kf,:)]=butter(N,[fone ftwo]);
end
nfiles
npts=86400*sps+1; % daily at 1 sps
%ftnsac=zeros(npts,3);
%tt=[1:npts]*dt0;ftnsac(:,1)=tt';
dt=dt0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define taper
taperwid = 2/fl/npts; % both sides*max period/(time of one day*sps)
w=[];
w=tukeywin(npts,taperwid);
CoreNum=4;
parpool('local',CoreNum);
parfor n=1:nfiles
% read sacfile
ftnsac=zeros(npts,3);
tt=[1:npts]*dt0;ftnsac(:,1)=tt';
sacfile=char(files{n});
data=[];
data=rsac(sacfile);
% sac file infor
ftnsac(:,3)=data(:,3);
% taper
data(:,2)=data(:,2).*w;  % taper the ends
%disp(['running time for filter:' num2str(toc)]);% ftn
%envnorm=0;
for kf=1:nf
narrowlyfiltered=[];
envnorm=[];
narrowlyfiltered=filtfilt(b(kf,:),a(kf,:),data(:,2));
% Envolope should be calculated with the following equation instead.
env=abs(hilbert(narrowlyfiltered));
envnorm=narrowlyfiltered./env;
ftnsac(:,2)=ftnsac(:,2)+envnorm;
end %foreach frequency band
%ftnsac(:,2)=envnorm;
%disp(['running time for hilbert:' num2str(toc)]);
% normalization
ftnsac(:,2)=ftnsac(:,2)/sqrt(nf);
indx=find(abs(ftnsac(:,2))>2);
ftnsac(indx,2)=0;

% taper the FTN seismogram
ftnsac(:,2)=ftnsac(:,2).*w;
%disp(['running time for normalization:' num2str(toc)]);
% write sac
ftnsacname=[outputstn '/ftn.' sacfile];
wsac(ftnsacname,ftnsac); %pause(0.5);

% re-set to zero for each new sac file
ftnsac(:,2)=0;

end % foreach sacfile
%clear
p = gcp;
delete(p);
!/bin/rm -f filelist.txt
cd ../
end
!/bin/rm -f listfolder.txt
EOF
