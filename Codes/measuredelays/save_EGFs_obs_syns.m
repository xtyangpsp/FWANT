% script to measure phase delays between the observed and synthetic waveforms
% modified from measure_phase_delay.csh
% Further modified by Xiaotao Yang @ UMASS
% Modifications:
% 1. changed to power relation ship time window length as a funtion of
% period;
% 2. SNR-based xcorrelation shift. if SNRs for the positive and negative
% data differ too much, use the larger SNR result. Otherwise, use weighted
% by SNR values to get the final shift (maximum correlation coefficient).
% 3. Removed cycle-skipping related selection criterion. This can be
% applied later when inspecting the overall measurements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Set up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wkdir = '/Users/xiaotaoyang/Work/FWANT/Alaska/ite_0.05deg_06';
%wkdir = '/Volumes/MacHD320/Work/UMass_PD/NewEnglandANT/FWANT/eUS/ite_0.015deg_02';
egfdir = '/Volumes/ProcessedXY1/Alaska/data/stacked_xcorr_cnv_a3t9';
syndir = [wkdir '/syn.seismograms'];
outdir = [wkdir '/measure/EGFs_compare'];
plotdir = strcat(wkdir,'/measure/plots_test');
system(['mkdir ' outdir]);
%%%%% define parameters %%%%
max_dV = 0.2; % maximum average velocity perturbation along the path (7%)
max_dT = 15; % maximum delay in second

fig_flag = 1;  % 1 == plot figure; else no figure
savefig=0;
saveresult=1; %1, save measurements to output file.
% define time windows
tminimum=5; % shortest period in second, should be the taper window length
tmaximum=1000; % length of synthetic green's function
tmaximumegf=1000; % length of egfs

%%% interpolate and decimate both syn and egf to the same sampling rate
dtuni=0.1;
cmin=2.5;cmax=4.5; % km/s group velocity
synerr=1.3; % about half of the grid spacing divided by group velocity (~0.05*100 km / 4 km/s)

snr_cutoff = 4.0; % snr limit
xcoeff_cutoff = 0.75; % cross correlation coefficient limit

%%% define the filter and taper for the egf and syn
taperfraction=tminimum/tmaximumegf*2; % taper tminimum at the ends
syntaperfraction=tminimum/tmaximum*2; % taper tminimum at the ends
N=2;
%pband=[75 150;50 100;35 75;25 50;15 30;10 20;7.5 15;5 10];
%fband=flip(1./pband,2);
% fband=[0.0067 0.01333;0.01 0.02;0.01333 0.0286;0.02 0.04;0.0286 0.0667; 0.05 .1; 0.0667 0.1333; .1 .2];
%Alaska fband
fband=[0.0067 .04];
%fband=[0.0067 0.01333];
%fband=[0.05 .1; 0.0667 0.1333; .1 .2];
tfmin=1./fband(:,1); % 1 period at the lower frequency limit of each frequency band

%addpath('~/FWT/ANT/Proj/cascadia/matlab/');
%warning off all
%stalist = {'TA', 'US', 'XY'};        
stalist = {'TA'}; 
%%%%%%%%%%%%%% END OF - Setting up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% poolobj=parpool('local');

for cc = 1:length(stalist)
    verylargenumber = 1.e9; 
        %if snr(k) > verylargenumber, it indicates err(k) = 0, so the estimate is invalid
    filelist=[];
    filelist=dir([egfdir '/' char(stalist{cc}) '*']);
    filenum=length(filelist);
    for ii=1:filenum

        stnpairname=filelist(ii).name;
        outfileobs_pos=[char(stnpairname) '_obs_pos'];
        outfileobs_neg=[char(stnpairname) '_obs_neg'];
        outfilesyn_pos=[char(stnpairname) '_syn_pos'];
        outfilesyn_neg=[char(stnpairname) '_syn_neg'];
        disp([num2str(ii),' --> ',stnpairname]);
%         if 1  % skip existing results. turn this off for reprocessing of the whole dataset.
        if exist([outdir '/' char(stnpairname) '.dat'],'file') ==2
            disp(['     ',outdir '/' char(stnpairname) '.dat   exists. Skipped!']);
        else
%         end
            tmp=strfind(stnpairname,'.');

            srcname=stnpairname(1:tmp(2)-1);
            rcvname=stnpairname(tmp(2)+1:end);

            %%%%%% INPUT EGFs %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% stack of all data (P = positive time lags; N = negative time lags)
            files=[]; sacdata=[]; egf_pos=[];egf_neg=[];tt=[];
            unix(['ls -f ' egfdir '/' stnpairname '/egf.*.all.P.SAC > filelist3_' stnpairname '_' num2str(ii)]);
            files=textread(['filelist3_' stnpairname '_' num2str(ii)] ,'%s');
            sacdata=readsac(char(files));
            egf_pos=sacdata.DATA1;
            %hdr=data(:,3);
            w=tukeywin(length(egf_pos),taperfraction);
            egf_pos=egf_pos.*w;

            files=[]; sacdata=[];
            unix(['ls -f ' egfdir '/' stnpairname '/egf.*.all.N.SAC > filelist3_' stnpairname '_' num2str(ii)]);
            files=textread(['filelist3_' stnpairname '_' num2str(ii)],'%s');
            sacdata=readsac(char(files));
            egf_neg=sacdata.DATA1.*w;
            ntrec=sacdata.NPTS; dt=sacdata.DELTA; tt=sacdata.B:dt:sacdata.E+dt;
            %tt=sacdata(:,1);ntrec=length(tt);dt=tt(2)-tt(1);

            unix(['rm filelist3_' stnpairname '_' num2str(ii)]);
            unix(['rm filelist2_' stnpairname '_' num2str(ii)]);

            %%% source-receiver information
            evla=sacdata.EVLA;evlo=sacdata.EVLO;stla=sacdata.STLA;stlo=sacdata.STLO;
            dist=geo2dist(evla,evlo,stla,stlo);
            dist_ellipse=geo2dist_ellipse(evla,evlo,stla,stlo);
            %lldistkm(evla,evlo,stla,stlo)
            ec=(dist-dist_ellipse)/cmax;  %correction for ellipticity for surface waves with cmax
            %ec=0;

            % define time windows
            tmin0=dist/cmax ;
            %if tmin0 < tminimum; tmin0=tminimum;end;
            if tmin0 < tminimum 
            %   disp('WARNING: Arrival time smaller than the minimum taper window time');
              syn_flag = 0;
            else
              syn_flag = 1;
            end
            tmin=tmin0-tfmin(1); % arrival time minus the longest period
            if tmin < tminimum; tmin=tminimum-4; end;
            if tmin > tmaximum; tmin=tmaximum-4; end;
            tmin0minustmin=tmin0-tmin;
            % use only arrivals within the simulation window-400 s to account for finite period, taper and source time function
            if tmin0 > tmaximum-100
            %   disp('WARNING: Arrival time greater than the maximum simulation time');
              syn_flag = 0;
            else
              syn_flag = 1;
            end
            %tmax=tmin0+2.0*tfmin(1);
            tmax=tmin0+2.0*tfmin(1);
            if tmax > tmaximum; tmax=tmaximum;end;
            tlen=tmax-tmin;
            %%%%%% END of EGF INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%% INPUT SYN displacement %%%%%%%%%%%%%%%%%%
            %if strmatch(sta, char(srcname))
            synfile = [wkdir '/syn.seismograms/' char(srcname) '/' char(srcname) '.to.' char(rcvname) '.fz.Uz.SAC'];
            %else
            %synfile = [wkdir '/syn.seismograms/' char(rcvname) '/' char(rcvname) '.to.' char(srcname) '.fz.Uz.SAC'];
            %endi

            % check if the synthetic data exists!
            synfid = fopen(synfile,'r');
            if synfid<=0, syn_flag = 0;
            else syn_flag = 1;fclose(synfid);
            end

            %%%%%%%%%
            if syn_flag == 1

                sacdata=[];ttsyn=[];syn=[];
                sacdata=readsac(synfile);
                %ttsyn=sacdata(:,1);dtsyn=ttsyn(2)-ttsyn(1);ntsyn=length(ttsyn);
                syn=sacdata.DATA1;
                ntsyn=sacdata.NPTS; dtsyn=sacdata.DELTA; ttsyn=sacdata.B:dtsyn:sacdata.E+dtsyn;

                if ~isempty(find(isnan(syn), 1)), continue, end

                %%%%%% END of SYN INPUT %%%%%%%%%%%%%%%%%%%%%%%%

                %%% filter EGFs
                [nfb, nc]=size(fband);
                egffb_pos=[];egfmonthfb_pos=[];
                egffb_neg=[];egfmonthfb_neg=[];

                itmin=round(tmin/dt);
                itmax=round(tmax/dt);
                Feff=1/dt/2; % in increment of 2

                for k=1:nfb
                  fone=fband(k,1)/Feff; ftwo=fband(k,2)/Feff;
                  [b,a]=butter(N,[fone ftwo]);
                  egffb_pos(:,k)=filtfilt(b,a,egf_pos);
                  egffb_neg(:,k)=filtfilt(b,a,egf_neg);
                end

                %%% filter syn 
                % demean and taper
                syn=syn-mean(syn);

                w=[];
                w=tukeywin(ntsyn,syntaperfraction); 
                syn=syn.*w;

                Feff=1/dtsyn/2;
                synfb=[];
                for k=1:nfb
                  fone=fband(k,1)/Feff; ftwo=fband(k,2)/Feff;
                  [b,a]=butter(N,[fone ftwo]);
                  synfb(:,k)=filtfilt(b,a,syn);
                end

                %%%%%% Cut and resample
                %%% cut egf to the length of the synthetics
                itsyn=round(ttsyn(end)/dt)+1;  % length of synthetics divided by dt of EGF
                egf_pos(itsyn+1:end)=[];
                egf_neg(itsyn+1:end)=[];
                egffb_pos(itsyn+1:end,:)=[];
                egffb_neg(itsyn+1:end,:)=[];
                tt(itsyn+1:end)=[];ntrec=length(tt);

                %%% interpolate and decimate both syn and egf to the same sampling rate
                % Note: select a dtnui that is smaller than data uncertainty but not too small to increase computation cost
                %dtuni=0.05; 
                %dtuni=0.2; % for a dtsyn of 3.5 s, dtsyn/0.02 is not a integer.  When interpolated with 
                % a roundoff number 18, dtuni for the synthetics is in fact 0.1944 s not 0.2 s, effectively stretching 
                % the synthetics.

                display(['interpolation rate of EGFs: ' num2str(dt/dtuni)])
                display(['interpolation rate of synthetics: ' num2str(dtsyn/dtuni)])
                rsyn=round(dtsyn/dtuni); regf=round(dt/dtuni);  %interpolation rate
                if abs(dtsyn/dtuni - rsyn) >= 0.001
                disp('WARNING: interpolation rate dtsyn/dtuni is not an integer and would result in roundoff errors');
                end
                if abs(dt/dtuni - regf) >= 0.001
                disp('WARNING: interpolation rate dt/dtuni is not an integer and would result in roundoff error');
                end

                tmpsyn=[];
                tmpegf_pos=[];
                tmpegf_neg=[];
                for k=1:nfb
                  tmpsyn(:,k)=interp(synfb(:,k),rsyn);
                  tmpegf_pos(:,k)=interp(egffb_pos(:,k),regf);
                  tmpegf_neg(:,k)=interp(egffb_neg(:,k),regf);
                end
                synfb=[];efgfb_pos=[];egffb_neg=[];
                synfb=tmpsyn;
                egffb_pos=tmpegf_pos; 
                egffb_neg=tmpegf_neg; 

                tmp_pos=[];
                tmp_neg=[];

                ttuni=0:dtuni:ttsyn(end);ntuni=length(ttuni);
                synfb(ntuni+1:end,:)=[];
                egffb_pos(ntuni+1:end,:)=[];
                egffb_neg(ntuni+1:end,:)=[];

                %%%%%% Define arrival of interest
                itmin=round(tmin/dtuni);
                itmax=round(tmax/dtuni);
                egfsig_pos=egffb_pos(itmin:itmax,:);
                egfsig_neg=egffb_neg(itmin:itmax,:);

                synsig=synfb(itmin:itmax,:);
                ntsig=length(egfsig_pos(:,1));

                t1=[]; t2=[];  tw_eff1=[]; tw_eff2=[];
                it1=[]; it2=[]; snr_pos=[]; snr_neg=[]; snr=[];
                for k=1:nfb
                  w=[]; w(1:ntsig,1)=0; %initialize taper of same length
                  % adjust the time window for each freq
                  %if tmin0minustmin > tfmin(k)/2; 
                  if tmin0minustmin > tfmin(k);
                  %t1(k)=tmin0minustmin - tfmin(k)/2;
                   t1(k)=tmin0minustmin - tfmin(k);
                  else
                    t1(k)=1;
                  end
                  if t1(k) < dtuni; t1(k)=dtuni;end
                  %t2(k)=tmin0minustmin+(1.0+1.5*k/nfb)*tfmin(k); % empirical 
                  %t2(k)=tmin0minustmin+(2.0+3.0*k/nfb)*tfmin(k); % empirical
                  t2(k)=tmin0minustmin+50*tfmin(k).^(-2/3)*tfmin(k)-25; % empirical
                  if ( t2(k) > tmax - tmin )
                    t2(k) = tmax - tmin;
                  end
                  if t1(k) > t2(k); t1(k)=t2(k)-1; end
                  tw_eff1(k)=tmin+t1(k);tw_eff2(k)=tmin+t2(k);
                  it1(k)=round(t1(k)/dtuni);it2(k)=round(t2(k)/dtuni);t1t2len=it2(k)-it1(k)+1;

                  w(it1(k):it2(k),1)=tukeywin(t1t2len,taperfraction);
                  egfsig_pos(:,k)=egfsig_pos(:,k).*w; % taper the selected time window
                  egfsig_neg(:,k)=egfsig_neg(:,k).*w;
                  synsig(:,k)=synsig(:,k).*w;

                end

                maxdelay = min(max_dV*tmin0,max_dT);
                %maxdelay=max_dT;
                if maxdelay < 1; maxdelay=1; end; % avoid a situation when tmin0 ~ 0
                maxlag=round(maxdelay/dtuni);
                ttc=(-maxlag:maxlag)*dtuni;

                xc_pos=[]; xc_neg=[]; xcm_pos=[]; xcm_neg=[];
                c_pos=[]; c_neg=[]; cm_pos=[]; cm_neg=[];
                for k=1:nfb
                    xc_pos(:,k)=xcorr(egfsig_pos(:,k),synsig(:,k),maxlag);
                    xc_neg(:,k)=xcorr(egfsig_neg(:,k),synsig(:,k),maxlag);
                    end
                end

                phase=[]; phaseerr=[];
                for k=1:nfb
                  [~, ixcmax_pos]=max(xc_pos(:,k));
                  [dum, ixcmax_neg]=max(xc_neg(:,k));
                  %select by snr
%                   if snr_neg(k) > verylargenumber || snr_pos(k) > verylargenumber
                      phase(k)=(ttc(ixcmax_pos)+ttc(ixcmax_neg))/2; 
                      % phase delay is the average of the positive and negative lags
%                   else
%                       if snr_neg(k)/snr_pos(k)<5/6
%                           phase(k)=ttc(ixcmax_pos);
%                       elseif snr_pos(k)/snr_neg(k)<5/6
%                           phase(k)=ttc(ixcmax_neg);
%                       else
%                           phase(k)=ttc(ixcmax_pos)*snr_pos(k)/(snr_pos(k)+snr_neg(k))+...
%                               ttc(ixcmax_neg)*snr_neg(k)/(snr_pos(k)+snr_neg(k)); %weight by snr
%                           %phase(k)=(ttc(ixcmax_pos)+ttc(ixcmax_neg))/2; % phase delay is the average of the positive and negative lags
%                       end
%                   end
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  phase(k)=phase(k)+ec; % with ellipticity correction!
                  %phase(k)
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  phaseerr(k)=0;  %error of the mean + sampling error for obs & syn
                  %%% add the errors due to synthetic waveforms
                  %%% inspection of reciprical synthetic waveforms show ~0.3 s shift, caused by the fact
                  %%% that the source & receiver are not exactly on grids and interpolation is needed.
                  %%% assuming data errors and synthetic errors are independent
                  %synerr=0.3;
                %   synerr=0.35; % half of the grid spacing divided by group velocity (~22/2 km / 4 km/s)
                  phaseerr(k)=sqrt(phaseerr(k)*phaseerr(k)+synerr*synerr);  % combine observation and synthetic errors
                end

                % cross-correlation coefficient
                rxc=[]; rc_pos=[]; rc_neg=[]; synshift=[];
                for k=1:nfb
                  itshift=round(phase(k)/dtuni);
                  if itshift == 0
                    synshift=synsig;
                  elseif itshift > 0
                    synshift(1:itshift,k)=0;
                    synshift(itshift+1:ntsig,k)=synsig(1:ntsig-itshift,k);
                  elseif itshift < 0
                    synshift(ntsig+itshift:ntsig,k)=0;
                    synshift(1:ntsig+itshift,k)=synsig(-itshift+1:ntsig,k);
                  end
                  rc_pos=corrcoef(egfsig_pos(:,k),synshift(:,k));
                  rc_neg=corrcoef(egfsig_neg(:,k),synshift(:,k));
                  rxc(k)=(rc_pos(1,2)+rc_neg(1,2))/2;
                end

                %correlation of the delay curves
                %xc_temp=corrcoef(xc_pos(:,k),xc_neg(:,k));
                %xc_delay(k)=xc_temp(1,2);
                %%%%% plot figure 
                if fig_flag==1 && syn_flag==1 

                    hid=figure('Position',[100 400 1200 1300]);
                    for k=1:nfb

                    % if rxc(k) > xcoeff_cutoff & snr(k) > snr_cutoff & snr(k) < verylargenumber & abs(phase(k)) < 0.9*maxdelay & abs(phase(k)) < 0.5*tfmin(k) & tmin0 > tfmin(k)

                    subplot(nfb,2,2*k)
                    %plot(ttc,xc_pos(:,k)/max(xc_pos(:,k)),'k-');hold on
                    plot(ttc+ec,xc_pos(:,k)/max(abs(xc_pos(:,k))),'k-');hold on % for ellipticity correction
                    %plot(ttc,xc_neg(:,k)/max(xc_neg(:,k)),'k--');hold on
                    plot(ttc+ec,xc_neg(:,k)/max(abs(xc_neg(:,k))),'k--');hold on % 
                    plot(phase(k),1,'rv'); hold on
                    text(-maxdelay+0.5,0.3,['lag: ' num2str(phase(k),3) '+/-' num2str(phaseerr(k),2)],'FontSize',10);hold on
                    %text(-maxdelay+0.5,0.3,['tmin: ' num2str(tmin,3)],'FontSize',5);hold on
                    text(-maxdelay+0.5,-0.3,['xcoeff: ' num2str(rxc(k),2)],'FontSize',10); hold on
%                     text(-maxdelay+0.5,-0.9,['snr: ' num2str(snr(k),3)],'FontSize',10);hold on

    %                text(maxdelay-1,0.3,['delayxcor: ' num2str(xc_delay(k),3)],'FontSize',10);hold on
                    % ite_03
                    % if tmin > tfmin(k) & rxc(k) > 0.5 & snr(k) > 4
                    % ite_04 after ite_02, the inspected absolute delays are less than 5 s 
                    % limiting phase < 0.95*maxdelay avoid measurements at the boundries of the time window,which
                    % are likely unreliable
    %                 if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
    %                         abs(phase(k)) <= 0.9*maxdelay && abs(phase(k)) <= 0.5*tfmin(k) && tmin0 >= tfmin(k)
%                     if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
%                             abs(phase(k)) <= 0.9*maxdelay && tmin0 >= tfmin(k)
%                     % use a lower snr for island stations II.COCO and II.DGAR
%                     %if tmin0 > tfmin(k) & rxc(k) > 0.85 & snr(k) > 4 & abs(phase(k)) < 0.90*maxdelay  
%                         text(+0.5,-0.9,'thumbs UP','Color',[0 0 1]);hold on
%                     else
%                         text(+0.5,-0.9,'thumbs DOWN','Color',[1 0 0]); hold on
%                     end
                    if k==nfb
                    xlabel('Lag time, s','FontSize',12);
                    end
                    axis([ttc(1) ttc(end) -1.3 1.3]);

                    % end
                    end

                    for k=1:nfb
                    % if rxc(k) > xcoeff_cutoff & snr(k) > snr_cutoff & snr(k) < verylargenumber & abs(phase(k)) < 0.9*maxdelay & abs(phase(k)) < 0.5*tfmin(k) & tmin0 > tfmin(k)

                    subplot(nfb,2,2*(k-1)+1)
                    %for im=1:imonth
                    %plot(ttuni,egfmonthfb_pos(:,im,k)/max(egfmonthfb_pos(itmin:itmax,im,k)),'Color',[0.7 0.7 0.7]);  hold on
                    %plot(ttuni,egfmonthfb_neg(:,im,k)/max(egfmonthfb_neg(itmin:itmax,im,k)),'Color',[0.7 0.7 0.7]);  hold on
                    %end
                    % % plot(ttuni,egffb_pos(:,k)/max(egffb_pos(itmin:itmax,k)),'k-'); hold on
                    % % plot(ttuni,egffb_neg(:,k)/max(egffb_neg(itmin:itmax,k)),'k--'); hold on
                    plot(ttuni,egffb_pos(:,k)/max(egffb_pos(:,k)),'k-'); hold on
                    plot(ttuni,egffb_neg(:,k)/max(egffb_neg(:,k)),'k--'); hold on
                    %plot(ttuni,synfb(:,k)/max(synfb(itmin:itmax,k)),'b'); hold on
                    %plot(ttuni-ec,synfb(:,k)/max(synfb(itmin:itmax,k)),'b'); hold on  % with ellipticity correction

                    plot(ttuni(itmin:itmax),synsig(:,k)/max(synsig(:,k)),'r');hold on
                    %plot(ttuni(itmin:itmax),synshift(:,k),'r');hold on
                    plot([tmin tmin],[-1.5 1.5],'g-','LineWidth',2); hold on
                    plot([tmin+t1(k) tmin+t1(k)],[-1.3 1.3],'b--'); hold on
                    plot([tmax tmax],[-1.5 1.5],'g-','LineWidth',2);  hold on
                    plot([tmin+t2(k) tmin+t2(k)],[-1.3 1.3],'b--');  hold on

                    text(ttuni(end)-150,0.7,[num2str(round(1/fband(k,2))) '-' num2str(round(1/fband(k,1))) 's'],'FontSize',14)

                    if k == nfb
                        xlabel('Time, s');
                    end
                    % axis([0 ttuni(end) -1.3 1.3]);
                    axis([tmin tmax -1.3 1.3]);

                    end

                    suptitle(strcat(char(srcname), ' to  ',char(rcvname),', ',num2str(dist),' km'));
                    if(savefig)
                        fignam=strcat(plotdir,'/',char(srcname), '_to_',char(rcvname),'_snr_',num2str(snr_cutoff),...
                            '_corcoe_',num2str(xcoeff_cutoff),'.eps');
                        saveas(gca,fignam,'epsc');
                        pause;
                    else
                        pause;
                    end
                    close all;
                    %print('-depsc',figname);
                    %print('-dpng',figname);

                    % end
                end % if fig_flag
                %%%%% END of plot figure
                if saveresult
                    %%%%% save phase measurements
                     outfilename = [outdir '/' char(stnpairname) '.dat'];
                     fid=fopen(outfilename,'r');
                     if fid>0; 
                         fclose(fid);
                         unix(['rm ' outfilename]);
                     end
                     for k=1:nfb
                         xcid=[char(srcname) '/bp' num2str(fband(k,1)) '_' num2str(fband(k,2)) '/' char(srcname) '_' char(rcvname) '_BHZ.P2.CORR.T1T2.SAC'];
                         fbid=['f' num2str(k)];
                         %%%%%% Normal stations %%%%%%%%%%%%%%%%%%%%%
                    %                  if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
                    %                        abs(phase(k)) <= 0.9*maxdelay && abs(phase(k)) <= 0.5*tfmin(k) && tmin0 >= tfmin(k)
                     end
                end
            end % of syn_flag ==1
        end
     end
%     %delete(poolobj);




