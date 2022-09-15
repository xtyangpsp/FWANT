function plt_data_vs_syn(datadir,syndir,evname, fband, subplotpar,my_ylim)
% matlab script to plot EGFs by distance
% Modified from Haiying Gao's plt_egfs_by_dist.m
% By Xiaotao Yang @ UMass Amherst on 11/23/2016
% Modifications:
%   1. change rsac and wsac to readsac and writesac, respectively;
%   2. plot maximum of 6 panels each run, which allows for 6
%   period/frequency bands;
%   3. reorganized the code to put global parameters at the beginning.
%   4. minimized the number of hard-coded variables.
%
%clear all;close all;
ylimtype='~';
if nargin < 5
    error('**Not enough input arguments!');
elseif nargin == 5 || numel(my_ylim) < 2
    ylimtype='auto';
end

path('/usr/local/sac/utils',path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Set up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R0=6371.0;
%datadir = '/Volumes/RAID0/hgao/cascadia_crust/stacked_xcorr_2sides';
%datadir = '/Volumes/NewEngland/stacked_xcorr_2sides';
%datadir = '/Volumes/NewEngland/stacked_xcorr_cnv_a1.5t4.5';
%fignameappend='cnv';

%sta='CN.A11';
%sta='UW.BURN';
%sta='TA.E06A'
% lon1 = 360-76.5; lon2 = 360-66.5; lat1 = 38; lat2 = 48;

snrmin=20;
% initialize
nptmax=5001; nptshalf=(nptmax-1)/2+1;
%nptmaxsyn=2262;nptshalfsyn=(nptmaxsyn-1)/2+1;
nptshalfsyn=1131;

% define signal window
cmin=2.5;cmax=4.5;
% define noise time length
noisetime = 100; % s

N=2;
%fband=[0.0067 0.01333;0.01 0.02;0.01333 0.0286;0.02 0.04;0.0286 0.0667; 0.05 .1; 0.0667 0.1333; .1 .2];
pband=flip(1./fband,2);
%fband=[0.005 0.01; 0.04 0.1];
%%%%%%%%%%%%%% END OF - Setting up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% stack of all data (P = positive time lags; N = negative time lags)
unix(['ls -f ' datadir '/' evname '*/egf.*.all.P.SAC > filelist1']);
files1=textread('filelist1' ,'%s');
% 
% unix(['ls -f ' datadir '/*' sta '*/egf.*.all.N.SAC > filelist2']);
% files2=textread('filelist2' ,'%s');

nfiles=length(files1);

% the "source"
% unix(['ls -f *' sta '*/*ZZ.sym.SAC > filelist ']);
% files=textread('filelist','%s');
% nfiles=length(files);

[nfb, ~]=size(fband);
tfmin=1./fband(:,1); % 1 period at the lower frequency limit of each frequency band

egffb_pos=zeros(nptshalf,nfiles,nfb);
synfbdata=zeros(nptshalfsyn,nfiles,nfb);
synsignalwin=nan(nptshalfsyn,nfiles,nfb);
%egffb_neg=zeros(nptshalf,nfiles,nfb);
for i=1:nfiles
    % read EGF file
    filename1=char(files1(i));
    itemp=strfind(filename1,'egf');
    itempZZ=strfind(filename1,'ZZ.all');
    
    clear stname;
    stname=filename1(itemp+5+length(evname):itempZZ-2);
    
    dd1=readsac(filename1);
    dt=dd1.DELTA;
    tt=dd1.B:dt:dd1.E+dt;
    ntrec=dd1.NPTS;
    Feff=1/dt/2; % in increment of 2

    evla(i)=dd1.EVLA;evlo(i)=dd1.EVLO;
    stlat(i)=dd1.STLA;stlon(i)=dd1.STLO;

     %[dist0(i),~]=distance(evla(i),evlo(i),stlat(i),stlon(i));
    dist0(i)=dd1.DIST;%round(R0*dist0(i)*pi/180.0);
    azim(i)=dd1.AZ;
    tmaximum=nptshalf*dt;
    tmin=dist0(i)/cmax;tmax=dist0(i)/cmin;itmin=round(tmin/dt);itmax=round(tmax/dt);
%     if itmax > ntrec; itmax = ntrec; end
%     if itmin <1; itmin=1; end
    for k=1:nfb
        taperfraction=tfmin(k)/tmaximum*2; % taper tminimum at the ends
        w=tukeywin(ntrec,taperfraction); 
        egffb_pos(:,i,k)=dd1.DATA1.*w;
%         egffb_neg(:,i,k)=dd2.DATA1.*w;

        fone=fband(k,1)/Feff; ftwo=fband(k,2)/Feff;
        [b,a]=butter(N,[fone ftwo]);
        egffb_pos(:,i,k)=filtfilt(b,a,egffb_pos(:,i,k));
%         egffb_neg(:,i,k)=filtfilt(b,a,egffb_neg(:,i,k));

        % define signal-to-noise ratio
        tempdata=egffb_pos(:,i,k);
        itmin=round((tmin-mean(pband(k,:)))/dt);itmax=round((tmax+mean(pband(k,:)))/dt);
        if itmax > ntrec; itmax = ntrec; end
        if itmin <1; itmin=1;end
        %signalwin=itmax-itmin+1;
        signal=max(abs(tempdata(itmin:itmax)));
%         if itmin > 25
            noise=max(std(abs(tempdata(1:itmin))),std(abs(tempdata(itmax:end))));
%         else
%             noise=std(abs(tempdata(itmax:end-25)));
%         end
        %noise=max(abs(tempdata(end-25-(itmax-itmin):end-25)));  %coda at the end as an indicator of noise
        mysnr(k,i)=signal/noise;
        clear tempdata;
    end
    %TA.G61A.to.TA.O61A.fz.Uz.SAC
    synfile=strcat(syndir,'/',evname,'/',evname,'.to.',stname,'.fz.Uz.SAC');
    dd3=readsac(synfile);
    dt3=dd3.DELTA;
    tt3=dd3.B:dt3:dd3.E;
    ntrec3=dd3.NPTS;
    Feff3=1/dt3/2; % in increment of 2
    
    tmaximum3=nptshalfsyn*dt3;
    
    clear w;
    for j=1:nfb
        taperfraction3=tfmin(j)/tmaximum3*2; % taper tminimum at the ends
        w=tukeywin(ntrec3,taperfraction3); 
%         size(w)
%         nptshalfsyn
        synfbdata(:,i,j)=dd3.DATA1.*w;
%         egffb_neg(:,i,k)=dd2.DATA1.*w;

        fone=fband(j,1)/Feff3; ftwo=fband(j,2)/Feff3;
        [b,a]=butter(N,[fone ftwo]);
        synfbdata(:,i,j)=filtfilt(b,a,synfbdata(:,i,j));
        
        itmin3=round((tmin-mean(pband(j,:)))/dt3);itmax3=round((tmax+mean(pband(j,:)))/dt3);
        if itmax3 > ntrec3; itmax3 = ntrec3; end
        if itmin3 <1; itmin3=1;end
        synsignalwin(itmin3:itmax3,i,j)=tt3(itmin3:itmax3);
    end
    
%     ftemp=corrcoef(egffb_pos(:,i,k),synfbdata(:,i,k));
%     mycorrcoe(k,i)=ftemp(1,2);
end

%idx = find(mysnr(1,:)>=snrmin & stlat>=lat1 & stlat<=lat2 & (stlon+360)>=lon1 & (stlon+360)<=lon2);

%%
figure('Position',[400 400 800 500]);
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) '};
nplotmax=8;
for k=1:nfb
    clear idx
    idx = find(mysnr(k,:)>=snrmin);
    %idx=find(dist0>=400 & azim>=225 & azim<=270);
    %nfig=1;
    if(k>nplotmax), error(['Only plot ',num2str(nplotmax),' frequency band for now.']), end
%     if k==1
%         figure(nfig);
        %subplot('Position',[0.5 0.55 0.45 0.35]),
        subplot(subplotpar(1),subplotpar(2),k),
        axis on, hold on, box on
        xlabel('Cross-correlation time (s)  ','FontSize',14);
        ylabel('Inter-station distance (km)  ','FontSize',14);
        title([figlabel{k},' EGFs v.s. Synthetics at: ', num2str(round(pband(k,1))), ...
            '-', num2str(round(pband(k,2))) ' s'],'FontSize',14)
        scale=10;
        %set(gca,'XLim',[-500 500])
        set(gca,'XLim',[0 300])
        if strcmp(ylimtype,'auto')
            %set(gca,'YLim',round([220 270]))
            set(gca,'YLim',round([max(mean(pband(k,:))*cmin,min(dist0)) 700]))
        else
            set(gca,'YLim',round(my_ylim))
        end

    for jj=1:length(idx)
        i=idx(jj);
        tmin(i)=dist0(i)/cmax;tmax(i)=dist0(i)/cmin;
        if k==1
          if mean(egffb_pos(round(tmin(i):(tmin(i)+tmax(i))/2),i,k))<=mean(egffb_pos(round((tmin(i)+(tmax(i))/2):tmax(i)),i,k)), continue, end
          %if mean(egffb_neg(round(tmin(i):(tmin(i)+tmax(i))/2),i,k))<=mean(egffb_neg(round((tmin(i)+(tmax(i))/2):tmax(i)),i,k)), continue, end
        end
        egffb_pos(:,i,k)=egffb_pos(:,i,k)/max(egffb_pos(:,i,k));
        %egffb_neg(:,i,k)=egffb_neg(:,i,k)/max(egffb_neg(:,i,k));
        plot(tt,scale*egffb_pos(:,i,k)+dist0(i),'k'); 
        synfbdata(:,i,k)=synfbdata(:,i,k)/max(synfbdata(:,i,k));
        plot(synsignalwin(:,i,k),scale*synfbdata(:,i,k)+dist0(i),'r-','Linewidth',1.5);
        %plot(-tt,scale*egffb_neg(:,i,k)+dist0(i),'k'); 
    end
    
    %plot([0 0],[min(dist0)-50 max(dist0)+50],'k--','LineWidth',2)
    hold off;
    legend('Observed','Synthetic','Location','southeast');
    set(gca,'FontSize',14);
    
    drawnow;
end

% figname = strcat('EGFs_', char(sta), '_',  num2str(snrmin), '.eps');
% set(gcf,'PaperPositionMode','auto');   
% eval(['print -depsc ' figname ]);
%unix(['mv ' figname ' /Users/hgao/FWT/ANT/Proj/cascadia/figures/EGFs/']);


%cd /Users/hgao/FWT/ANT/Proj/cascadia/matlab/
end