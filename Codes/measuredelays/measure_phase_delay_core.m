%%% insert parameters for measuring phase delays.
set_measure_parameters;

%%% read station information
fid = fopen(stainfo);
recv = textscan(fid,'%s%s%*f%*f%*f'); %ntwk stnm lon lat elevation
fclose(fid);
nrecv=length(recv{1});
ntwk=recv{1}; stnm=recv{2}; 
sourcelist=cell(nrecv,1);
for nstsrc=1:nrecv %
    sourcelist{nstsrc}=[char(ntwk(nstsrc)) '.' char(stnm(nstsrc))];
end
nsource=length(sourcelist);
%%%
if ~exist(outdir,'dir')
    system(['mkdir ' outdir]);
end
if ~exist(plotdir,'dir')
    system(['mkdir ' plotdir]);
end
synfileext='.fz.Uz.SAC';
egf_comp='ZZ';

fband=flip(flip(1./pband),2);

tfmin=1./fband(:,1); % 1 period at the lower frequency limit of each frequency band

synerr=0.05*110/2/mean([cmin,cmax]); % about half of the grid spacing divided by group velocity (~0.05*110 km /2/ 4.5 km/s)
verylargenumber = 1.e9; %to check numerically unstable values. 

%%% define the filter and taper for the egf and syn
taperfraction=tminimum/tmaximumegf; % taper tminimum at the ends
syntaperfraction=tminimum/tmaximum; % taper tminimum at the ends
N=2;

%%%%%%%%%%%%%% END OF - Setting up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
parfor ii = 1:nsource
    source = sourcelist{ii};
    disp([num2str(ii),' --> ',source]);
    %get synthetic file list for the source.
    synfilelist=dir([syndir,'/',source,'/*',synfileext]);

    %loop through all pairs.
    npairs=length(synfilelist);
    if saveresult
        outfilename = [outdir '/' char(source) '.dat'];
        if exist(outfilename,'file')
            unix(['rm ' outfilename]);
        end
        fidout = fopen(outfilename,'a+');
    end
    for n =1:npairs
        synfile=synfilelist(n).name;
        stemp=strsplit(synfile(1:end - length(synfileext)),'.to.');
        src=stemp{1};
        rcv=stemp{2};
        
        if strcmp(src,rcv); continue;end
        stnpair=[src,'_',rcv];
        fntemp=dir([egfdir,'/',source,'/',stnpair,'*',egf_comp,'_N_stack.h5']);
        if isempty(fntemp); continue; end
        %log
        disp(['working on pair: ' stnpair '---> ' num2str(n) '/' num2str(npairs)])
      
        egffile_stack_N=fntemp.name;
        fntemp=dir([egfdir,'/',source,'/',stnpair,'*',egf_comp,'_P_stack.h5']);
        egffile_stack_P=fntemp.name;

        fntemp=dir([egfdir,'/',source,'/',stnpair,'*',egf_comp,'_N.h5']);
        egffile_N=fntemp.name;
        fntemp=dir([egfdir,'/',source,'/',stnpair,'*',egf_comp,'_P.h5']);
        egffile_P=fntemp.name;
        %read in data.
        synsacdata=readsac([syndir,'/',source,'/',synfile]);
        cdataall=extract_corrdata_asdf([egfdir,'/',source,'/',egffile_stack_N]);
        cdata_stack_N=cdataall.value{1}.(egf_comp);
        cdataall=extract_corrdata_asdf([egfdir,'/',source,'/',egffile_stack_P]);
        cdata_stack_P=cdataall.value{1}.(egf_comp);

        cdataall=extract_corrdata_asdf([egfdir,'/',source,'/',egffile_N]);
        cdata_substack_N=cdataall.value{1}.(egf_comp);
        cdataall=extract_corrdata_asdf([egfdir,'/',source,'/',egffile_P]);
        cdata_substack_P=cdataall.value{1}.(egf_comp);

        nsubstack=length(cdata_substack_N.time);
        ntrec=size(cdata_stack_N.data,1); dt=cdata_stack_N.dt; tt=cdata_stack_N.tvec;

        % apply cutoff to number of substacks.
        if nsubstack < min_substack; continue;end

        %apply taper to data.
        w=tukeywin(length(cdata_stack_P.data),taperfraction);
        egf_pos_temp=cdata_stack_P.data.*w;
        egf_neg_temp=cdata_stack_N.data.*w;

        %%% substacks     
        egf_substack_pos_temp=cdata_substack_P.data.*w;
        egf_substack_neg_temp=cdata_substack_N.data.*w;

        %%%%%%%%%%%%% IMPORTANT IF EGFs WERE CONVOLVED WITHOUT SHIFT
        %%%%%%%%%%%%% (ZERO-PHASE). Shift EGFs here with backward for 0.5*waveletshift.
        shiftsamples=int16(waveletshift/dt);
        egf_pos=zeros(length(egf_pos_temp),1);
        egf_neg=zeros(length(egf_neg_temp),1);
        egf_pos(shiftsamples+1:end)=egf_pos_temp(1:end-shiftsamples);
        egf_neg(shiftsamples+1:end)=egf_neg_temp(1:end-shiftsamples);

        egf_substack_pos = zeros(size(egf_substack_pos_temp));
        egf_substack_neg = zeros(size(egf_substack_neg_temp));
        egf_substack_pos(shiftsamples+1:end,:)=egf_substack_pos_temp(1:end-shiftsamples,:);
        egf_substack_neg(shiftsamples+1:end,:)=egf_substack_neg_temp(1:end-shiftsamples,:);
        
        %%% source-receiver information
        evla=cdata_substack_N.lat(1);evlo=cdata_substack_N.lon(1);
        stla=cdata_substack_N.lat(2);stlo=cdata_substack_N.lon(2);
        dist=cdata_substack_N.dist;
        dist_ellipse=geo2dist_ellipse(evla,evlo,stla,stlo);
        ec=(dist-dist_ellipse)/cmax;  %correction for ellipticity for surface waves with cmax
        
        %%%
        syn=synsacdata.DATA1;
        ntsyn=synsacdata.NPTS; dtsyn=synsacdata.DELTA; 
        ttsyn=synsacdata.B:dtsyn:synsacdata.E+dtsyn;

        if ~isempty(find(isnan(syn), 1)), continue, end

        %%% filter EGFs
        [nfb, nc]=size(fband);
        egffb_pos=nan(length(egf_pos),nfb);egffb_substack_pos=nan(length(egf_pos),nsubstack,nfb);
        egffb_neg=nan(length(egf_pos),nfb);egffb_substack_neg=nan(length(egf_pos),nsubstack,nfb);

        Feff=1/dt/2; % in increment of 2

        for k=1:nfb
            fone=fband(k,1)/Feff; ftwo=fband(k,2)/Feff;
            [b,a]=butter(N,[fone ftwo]);
            egffb_pos(:,k)=filtfilt(b,a,egf_pos);
            egffb_neg(:,k)=filtfilt(b,a,egf_neg);
            for im=1:nsubstack
                egffb_substack_pos(:,im,k)=filtfilt(b,a,egf_substack_pos(:,im));
                egffb_substack_neg(:,im,k)=filtfilt(b,a,egf_substack_neg(:,im));
            end
        end

        %%% filter syn 
        % demean and taper
        syn=syn-mean(syn);

        w=tukeywin(ntsyn,syntaperfraction); 
        syn=syn.*w;

        Feff=1/dtsyn/2;
        synfb=nan(length(syn),nfb);
        for k=1:nfb
          fone=fband(k,1)/Feff; ftwo=fband(k,2)/Feff;
          [b,a]=butter(N,[fone ftwo]);
          synfb(:,k)=filtfilt(b,a,syn);
        end

        %%%%%% Cut and resample
        %%% cut egf to the length of the synthetics
        itsyn=round(ttsyn(end)/dt)+1;  % length of synthetics divided by dt of EGF
        egf_substack_pos(itsyn+1:end,:)=[];egf_pos(itsyn+1:end)=[];
        egf_substack_neg(itsyn+1:end,:)=[];egf_neg(itsyn+1:end)=[];
        egffb_pos(itsyn+1:end,:)=[];egffb_substack_pos(itsyn+1:end,:,:)=[];
        egffb_neg(itsyn+1:end,:)=[];egffb_substack_neg(itsyn+1:end,:,:)=[];
        tt(itsyn+1:end)=[];
        
        %%% interpolate and decimate both syn and egf to the same sampling rate
        % Note: select a dtnui that is smaller than data uncertainty but not too small to increase computation cost
        %dtuni=0.05; 
        %dtuni=0.2; % for a dtsyn of 3.5 s, dtsyn/0.02 is not a integer.  When interpolated with 
        % a roundoff number 18, dtuni for the synthetics is in fact 0.1944 s not 0.2 s, effectively stretching 
        % the synthetics.

        display(['interpolation rate of EGFs: ' num2str(dt/dt_resample)])
        display(['interpolation rate of synthetics: ' num2str(dtsyn/dt_resample)])
        rsyn=round(dtsyn/dt_resample); regf=round(dt/dt_resample);  %interpolation rate
        if abs(dtsyn/dt_resample - rsyn) >= 0.001
        disp('WARNING: interpolation rate dtsyn/dtuni is not an integer and would result in roundoff errors');
        end
        if abs(dt/dt_resample - regf) >= 0.001
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
        for k=1:nfb
          for im=1:nsubstack
            tmp_pos(:,im,k)=interp(egffb_substack_pos(:,im,k),regf);
            tmp_neg(:,im,k)=interp(egffb_substack_neg(:,im,k),regf);
          end
        end
        egffb_substack_pos=[];egffb_substack_pos=tmp_pos; 
        egffb_substack_neg=[];egffb_substack_neg=tmp_neg; 

        ttuni=0:dt_resample:ttsyn(end);ntuni=length(ttuni);
        synfb(ntuni+1:end,:)=[];
        egffb_pos(ntuni+1:end,:)=[];egffb_substack_pos(ntuni+1:end,:,:)=[];
        egffb_neg(ntuni+1:end,:)=[];egffb_substack_neg(ntuni+1:end,:,:)=[];

        %%%%%% Define arrival of interest
        % define time windows
        tmin0=dist/cmax;

        tmin=tmin0-tfmin(1); % arrival time minus the longest period
        if tmin < tminimum; tmin=tminimum-4; end
        if tmin > tmaximum; tmin=tmaximum-4; end
        tmin0minustmin=tmin0-tmin;
        % use only arrivals within the simulation window-400 s to account for finite period, taper and source time function
        if tmin0 < pband(1,1) || tmin0 > tmaximum-tfmin(end)
            continue
        end

        tmax=dist/cmin + tfmin(1);%tmin0+2.0*tfmin(1);
        if tmax > tmaximum; tmax=tmaximum;end

        itmin=round(tmin/dt_resample);
        itmax=round(tmax/dt_resample);
        egfsig_pos=egffb_pos(itmin:itmax,:);
        egfsig_neg=egffb_neg(itmin:itmax,:);

        egfsubstacksig_pos=egffb_substack_pos(itmin:itmax,:,:);
        egfsubstacksig_neg=egffb_substack_neg(itmin:itmax,:,:);
        synsig=synfb(itmin:itmax,:);
        ntsig=length(egfsig_pos(:,1));

        t1=nan(nfb,1); t2=nan(nfb,1);  tw_eff1=nan(nfb,1); tw_eff2=nan(nfb,1);
        it1=nan(nfb,1); it2=nan(nfb,1); snr_pos=nan(nfb,1); snr_neg=nan(nfb,1); snr=nan(nfb,1);
        for k=1:nfb
          w=[]; w(1:ntsig,1)=0; %initialize taper of same length
          % adjust the time window for each freq
          %if tmin0minustmin > tfmin(k)/2; 
          if tmin0minustmin > tfmin(k)
            t1(k)=tmin0minustmin - tfmin(k);
          else
            t1(k)=1;
          end
          if t1(k) < dt_resample; t1(k)=dt_resample;end

%           t2(k)=tmin0minustmin+50*tfmin(k).^(-2/3)*tfmin(k) + 2*tfmin(k); % empirical
          t2(k)=tmin0minustmin+2*tfmin(k)+k*40;
          if ( t2(k) > tmax - tmin )
            t2(k) = tmax - tmin;
          end
          if t1(k) > t2(k); t1(k)=t2(k)-1; end
          tw_eff1(k)=tmin+t1(k);tw_eff2(k)=tmin+t2(k);
          it1(k)=round(t1(k)/dt_resample);it2(k)=round(t2(k)/dt_resample);t1t2len=it2(k)-it1(k)+1;

          w(it1(k):it2(k),1)=tukeywin(t1t2len,taperfraction);
          egfsig_pos(:,k)=egfsig_pos(:,k).*w; % taper the selected time window
          egfsig_neg(:,k)=egfsig_neg(:,k).*w;
          synsig(:,k)=synsig(:,k).*w;
          for im=1:nsubstack
            egfsubstacksig_pos(:,im,k)=egfsubstacksig_pos(:,im,k).*w;
            egfsubstacksig_neg(:,im,k)=egfsubstacksig_neg(:,im,k).*w;
          end

          err_pos=zeros(ntsig,1);err_neg=zeros(ntsig,1);
          for i=1:ntsig
            err_pos(i)=std(egfsubstacksig_pos(i,:,k))/sqrt(nsubstack);
            err_neg(i)=std(egfsubstacksig_neg(i,:,k))/sqrt(nsubstack);
          end
          snr_pos(k)=max(abs(egfsig_pos(:,k)))/mean(abs(egffb_pos(:,k)));
          snr_neg(k)=max(abs(egfsig_neg(:,k)))/mean(abs(egffb_neg(:,k)));
%           snr_pos(k)=max(abs(egfsig_pos(:,k)))/max(err_pos);
%           snr_neg(k)=max(abs(egfsig_neg(:,k)))/max(err_neg);

          snr(k)=max([snr_pos(k),snr_neg(k)]);%(snr_pos(k)+snr_neg(k))/2; % average snr from positive and negative lags
        end
        %mean snr for all frequency bands.
        snr_pos_mean=mean(snr_pos);
        snr_neg_mean=mean(snr_neg);
        use_side='m';
%         if snr_neg_mean < verylargenumber && snr_pos_mean < verylargenumber
%             if snr_neg_mean<snr_pos_mean
%                 use_side='p';
%             elseif snr_pos_mean<snr_neg_mean 
%                 use_side='n';
%             end
%             if snr_neg_mean/snr_pos_mean<2/3
%                 use_side='p';
%             elseif snr_pos_mean/snr_neg_mean < 2/3
%                 use_side='n';
%             else
%                 use_side='m';
%             end
%         end
        maxdelay = min(max_dV*tmin0,max_dT);
%         maxdelay=max_dT;
        if maxdelay < 1; maxdelay=1; end % avoid a situation when tmin0 ~ 0
        maxlag=round(maxdelay/dt_resample);
        ttc=(-maxlag:maxlag)*dt_resample;

        xc_pos=[]; xc_neg=[]; xcm_pos=[]; xcm_neg=[];
        c_pos=[]; c_neg=[]; cm_pos=[]; cm_neg=[];
        phase=nan(nfb,1); phaseerr=nan(nfb,1);
        
        rxc=nan(nfb,1); rc_pos=[]; rc_neg=[]; 
        synshift=synsig;
        for k=1:nfb
            xc_pos(:,k)=xcorr(egfsig_pos(:,k),synsig(:,k),maxlag);
            xc_neg(:,k)=xcorr(egfsig_neg(:,k),synsig(:,k),maxlag);
            for im=1:nsubstack
                xcm_pos(:,im,k)=xcorr(egfsubstacksig_pos(:,im,k),synsig(:,k),maxlag);
                xcm_neg(:,im,k)=xcorr(egfsubstacksig_neg(:,im,k),synsig(:,k),maxlag);
            end
        
            [~, ixcmax_pos]=max(xc_pos(:,k));
            [~, ixcmax_neg]=max(xc_neg(:,k));
            %use all frequency mean snr
            phase_temp=[];
            if strcmp(use_side,'m')
                phase_temp=mean([ttc(ixcmax_pos),ttc(ixcmax_neg)]);
            elseif strcmp(use_side,'p')
                phase_temp=ttc(ixcmax_pos);
            elseif strcmp(use_side,'n')
                phase_temp=ttc(ixcmax_neg);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            phase(k)=phase_temp+ec; % with ellipticity correction!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            phasem=nan(nsubstack,1);
            for im=1:nsubstack
                [~, ixcmax_pos]=max(xcm_pos(:,im,k));
                [~, ixcmax_neg]=max(xcm_neg(:,im,k));
%                 phasem(im)=(ttc(ixcmax_pos)+ttc(ixcmax_neg))/2;
                if strcmp(use_side,'m')
                    phasem(im)=mean([ttc(ixcmax_pos),ttc(ixcmax_neg)]);
                elseif strcmp(use_side,'p')
                    phasem(im)=ttc(ixcmax_pos);
                elseif strcmp(use_side,'n')
                    phasem(im)=ttc(ixcmax_neg);
                end
            end
            phaseerr_temp=std(phasem)/sqrt(nsubstack)+dt_resample;  
            %error of the mean + sampling error for obs & syn
            %%% add the errors due to synthetic waveforms
            %%% inspection of reciprical synthetic waveforms show ~0.3 s shift, caused by the fact
            %%% that the source & receiver are not exactly on grids and interpolation is needed.
            %%% assuming data errors and synthetic errors are independent
            phaseerr(k)=sqrt(phaseerr_temp*phaseerr_temp+synerr*synerr);  % combine observation and synthetic errors
        
            % cross-correlation coefficient
            itshift=round(phase(k)/dt_resample);
            if itshift > 0
                synshift(1:itshift,k)=0;
                synshift(itshift+1:ntsig,k)=synsig(1:ntsig-itshift,k);
            elseif itshift < 0
                synshift(ntsig+itshift:ntsig,k)=0;
                synshift(1:ntsig+itshift,k)=synsig(-itshift+1:ntsig,k);
            end
            rc_pos=corrcoef(egfsig_pos(:,k),synshift(:,k));
            rc_neg=corrcoef(egfsig_neg(:,k),synshift(:,k));
            rxc(k)=max([rc_pos(1,2),rc_neg(1,2)]); %(rc_pos(1,2)+rc_neg(1,2))/2;
        end

        %correlation of the delay curves
        %%%%% plot figure 
        if fig_flag==1  
            if savefig
                hid=figure('Position',[100 400 1200 1300],'visible','off');
            else
                hid=figure('Position',[100 400 1200 1300]);
            end
            for k=1:nfb
                subplot(nfb,2,2*(k-1)+1)
                
                ph1=plot(ttuni,egffb_pos(:,k)/max(abs(egffb_pos(:,k))),'k-'); hold on
                ph2=plot(ttuni,egffb_neg(:,k)/max(abs(egffb_neg(:,k))),'b--'); hold on
                
                ph3=plot(ttuni(itmin:itmax),synshift(:,k)/max(abs(synshift(:,k))),'r');hold on
                plot([tmin tmin],[-1.5 1.5],'k-','LineWidth',1.5); hold on
                plot([tmin+t1(k) tmin+t1(k)],[-1.3 1.3],'m--'); hold on
                plot([tmax tmax],[-1.5 1.5],'k-','LineWidth',1.5);  hold on
                plot([tmin+t2(k) tmin+t2(k)],[-1.3 1.3],'m--');  hold on

                if k==1
                    legend([ph1,ph2,ph3],'egf pos','egf neg','syn','Orientation','horizontal')
                end
                if k == nfb
                    xlabel('Time (s)');
                end
                ylabel([num2str(pband(nfb-k+1,1)),' - ',num2str(pband(nfb-k+1,2)),' s'])
                axis([tmin tmax -1.3 1.3]);

                %%%%%%%%%%%%%%%%%% measurements
                subplot(nfb,2,2*k)
                plot(ttc+ec,xc_pos(:,k)/max(abs(xc_pos(:,k))),'k-');hold on % for ellipticity correction
                plot(ttc+ec,xc_neg(:,k)/max(abs(xc_neg(:,k))),'b--');hold on % 
                plot(phase(k),1,'rv'); hold on
                text(-maxdelay+0.5,0.3,['lag: ' num2str(phase(k),3) '+/-' num2str(phaseerr(k),2)],'FontSize',10);hold on
                text(-maxdelay+0.5,-0.3,['xcoeff: ' num2str(rxc(k),2)],'FontSize',10); hold on
                text(-maxdelay+0.5,-0.9,['snr: ' num2str(snr(k),3)],'FontSize',10);hold on

                % limiting phase < 0.95*maxdelay avoid measurements at 
                % the boundries of the time window,which are likely unreliable
                if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
                        abs(phase(k)) <= 0.95*maxdelay && tmin0 >= tfmin(k) 
                    text(+0.5,-0.9,'thumbs UP','Color',[0 0 1]);hold on
                else
                    text(+0.5,-0.9,'thumbs DOWN','Color',[1 0 0]); hold on
                end
                if k==nfb
                    xlabel('Phase difference (s)','FontSize',12);
                end
                axis([ttc(1) ttc(end) -1.3 1.3]);

            end

            sgtitle([strrep(stnpair,'_','-'),': ',num2str(dist,'%.1f'),' km: side = ',use_side]);
            if(savefig)
                fignam=strcat(plotdir,'/',char(src), '_to_',char(rcv),'_snr_',num2str(snr_cutoff),...
                    '_corcoe_',num2str(xcoeff_cutoff),'.png');
                saveas(gca,fignam,'png');
            else
                pause;
            end
            close all;
        end % if fig_flag
        %%%%% END of plot figure
        if saveresult
            %%%%% save phase measurements
             for k=1:nfb
                 xcid=[char(src) '/bp' num2str(fband(k,1)) '_' num2str(fband(k,2)) '/' stnpair '_BHZ.P2.CORR.T1T2.SAC'];
                 fbid=['f' num2str(k)];
                 %%%%%% Normal stations %%%%%%%%%%%%%%%%%%%%%
                 if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
                        abs(phase(k)) <= 0.95*maxdelay && tmin0 >= tfmin(k)
                     
                     fprintf(fidout,'%s %6.3f %2.0f %s %7.2f %7.2f %s %5.2f %5.2f %5.1f\n',...
                         xcid,phase(k),1,'RL',tw_eff1(k),tw_eff2(k),fbid,phaseerr(k),rxc(k),snr(k)); 
                 end
             end
        end
    end
    if saveresult;fclose(fidout);end
end

