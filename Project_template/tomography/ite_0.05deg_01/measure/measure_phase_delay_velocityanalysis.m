% script to measure phase delays between the observed and synthetic waveforms
% modified from measure_phase_delay.csh
% Further modified by Xiaotao Yang
% Modifications: 
% 1. changed to use moveout-awared time window, with specified search range. - 2/20/2026
% 2. removed option to shift the EGFs, leave this to data preprocessing. -
% 2/20/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Set up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise/src'))
% PROJHOME = '/depot/xtyang/data/projects/xtyang/FloridanAquifer/RiverRise';
PROJHOME = '/Users/xtyang/Work/Research/Projects/FloridanAquifer/RiverRise';
ite = '/ite_0.00001deg_01';
wkdir = [PROJHOME ite];
egfdir = [PROJHOME '/data/data_riverrise/PAIRS_TWOSIDES_gaussian_a0.005t0.015'];
syndir = [wkdir '/syn.seismograms'];
outdir = [wkdir '/measure/measure_stnpair'];
plotdir = strcat(wkdir,'/measure/plots_test');
stainfo = [PROJHOME '/STinfo/station_riverrise_withdata_formatok.txt'];
%%%%%%%%%
%%%%%%%%%%%%%%%%%END OF BEHAVIOR CONTROL PARAMETERS
%period and frequency band
fband=[4 8; 6 10; 7 12; 9 15; 12 20; 15 25];
pband = flip(flip(1./fband),2);
[nfb, nc]=size(fband);
%%%%% define parameters %%%%
max_dV = 0.25;
max_dT = 0.3; % maximum delay in second
snr_cutoff = 4; % snr limit
xcoeff_cutoff = 0.6; % cross correlation coefficient limit
min_wavelength = 0.75; %minimum near-field wavelength to satisfy.

min_substack = 3; % minimum number of substacks.
%%% interpolate and decimate both syn and egf to the same sampling rate
dt_resample=0.002;
dt_egf = 0.01;
% define time windows
t_taper=10*max(dt_egf,dt_resample);% use 10 samples here. 
%0.2*min(min(pband)); % shortest period in second, should be the taper window length
tmaximum=4.0; % length of synthetic green's function
tmaximumegf=5.0; % length of egfs

%%%%%%%%%
v_search_grid = 0.2:0.01:2.0; % Velocity range to search (km/s) in velocity analysis

model_grid_spacing=0.00001; %simulation model grid spacing in degrees.
verylargenumber = 1.e9; %to check numerically unstable values. 

%%% define the filter and taper for the egf and syn
taperfraction=t_taper/tmaximumegf; % taper t_taper at the ends
syntaperfraction=t_taper/tmaximum; % taper t_taper at the ends
N=2;
%%%%%%%%%%%%%%%%%BEHAVIOR CONTROL PARAMETERS
fig_flag = 1;  % 1 == plot figure; else no figure
plot_vanalysis=fig_flag; %plot figs for velocity analysis with moveout plots.
savefig=0;
saveresult=1; %1, save measurements to output file.
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

tfmin=1./fband(:,1); % 1 period at the lower frequency limit of each frequency band

%
close all;
%%%%%%%%%%%%%% END OF - Setting up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for ii = 1%:nsource
    source = sourcelist{ii};
    disp([num2str(ii),' --> ',source]);
    %get synthetic file list for the source.
    synfilelist=dir([syndir,'/',source,'/*',synfileext]);

    %loop through all pairs.
    npairs=length(synfilelist);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- NEW: MOVEOUT-BASED VELOCITY ANALYSIS ---
    disp('    Performing velocity analysis ... ');
    
    % 1. Collect all valid EGFs for this source to build a moveout gather
    all_egfs = cell(npairs, 1);
    all_dists = zeros(npairs, 1);
    valid_pair_idx = [];

    for n = 1:npairs
        np=n;
        synfile = synfilelist(np).name;
        stemp = strsplit(synfile(1:end - length(synfileext)), '.to.');
        src = stemp{1}; rcv = stemp{2};
        if strcmp(src, rcv); continue; end
        
        stnpair = [src, '_', rcv];
        % Check for observed data files
        fntemp_N = dir([egfdir, '/', source, '/', stnpair, '*', egf_comp, '_N_stack.h5']);
        fntemp_P = dir([egfdir, '/', source, '/', stnpair, '*', egf_comp, '_P_stack.h5']);
        
        if ~isempty(fntemp_N) && ~isempty(fntemp_P)
            % Load and stack the two sides (Positive and Negative lags)
            d_N = extract_corrdata_asdf([egfdir, '/', source, '/', fntemp_N.name]);
            d_P = extract_corrdata_asdf([egfdir, '/', source, '/', fntemp_P.name]);
            
            % Symmetric stack to enhance the surface wave signal
            % Note: Negative lag must be flipped to align with positive lag
            combined_data = (d_P.value{1}.(egf_comp).data + flipud(d_N.value{1}.(egf_comp).data)) / 2;
            
            all_egfs{np} = combined_data;
            all_dists(np) = d_P.value{1}.(egf_comp).dist;
            valid_pair_idx = [valid_pair_idx; np];
        end
    end
    % subset to only store the valid pairs
    all_egfs=all_egfs(valid_pair_idx);
    all_dists=all_dists(valid_pair_idx);
    npairs_valid=length(valid_pair_idx);
    npts_egf=length(all_egfs{1});
    taxis_egf=(0:npts_egf-1)*dt_egf;

    %
    disp(['    Num of EGF pairs: ',num2str(npairs_valid)])
    if npairs_valid < 1; disp('    -> skipped.'); continue; end
    w_taper_egf = tukeywin(npts_egf,taperfraction);
    %
    if plot_vanalysis
        % plot moveout
        figure; hold on;
        
        for ip = 1:npairs_valid
            goodidx=ip;
            d_temp = all_egfs{goodidx}.*w_taper_egf;
            plot(taxis_egf,all_dists(goodidx)+0.01*max(all_dists)*d_temp./max(abs(d_temp)),'k')
        end
        hold off;
        box on;
        title(['moveout: ',source])
        xlabel('time (s)')
        ylabel('distance (km)')
        set(gca,'TickDir','out')
        drawnow;

    end
        
    % 2. Grid search for best velocity per frequency band
    best_v_source = zeros(nfb, 1);
    if plot_vanalysis
        figure('Position',[400 400 1000 650]);
    end
    for kf = 1:nfb
        k=kf;
        
        energy_results = zeros(length(v_search_grid), 1);
        [b, a] = butter(N, [fband(k,1) fband(k,2)] / (1/dt_egf/2));
        
        for iv = 1:length(v_search_grid)
            v_test = v_search_grid(iv);
            stack_energy = 0;
            
            for idx = 1:length(valid_pair_idx)
                np = idx;
                temp_data = filtfilt(b, a, all_egfs{np});
                temp_data = temp_data.*w_taper_egf;
                
                % Determine window based on moveout
                t_arrival = all_dists(np) / v_test;
                it_arrival = round(t_arrival / dt_egf);
                win_samples = round(tfmin(k) / dt_egf);
                
                % Calculate energy (envelope) in the predicted window
                if it_arrival > win_samples && (it_arrival + win_samples) < length(temp_data)
                    env = abs(hilbert(temp_data(it_arrival-win_samples : it_arrival+win_samples)));
                    stack_energy = stack_energy + max(env);
                end
            end
            energy_results(iv) = stack_energy;
        end
        [~, max_idx] = max(energy_results);
        best_v_source(k) = v_search_grid(max_idx);

        %
        if plot_vanalysis
            subplot(2,3,k); hold on;
            
            for idx = 1:npairs_valid
                np = idx;
                temp_data = filtfilt(b, a, all_egfs{np});
                temp_data = temp_data.*w_taper_egf;
                plot(taxis_egf,all_dists(np)+0.01*max(all_dists)*temp_data./max(abs(temp_data)),'k')
            end
            %plot the best velocity line
            plot(taxis_egf,best_v_source(k)*taxis_egf,'r')
            ylim([0 max(all_dists)])
            box on;
            hold off
            title([source,' at ',num2str(fband(k,1)),'-',num2str(fband(k,2)),' Hz: ',num2str(best_v_source(k)),' km/s'])
            xlabel('time (s)')
            ylabel('distance (km)')
            set(gca,'TickDir','out')
%             drawnow;
        end
    end
    disp(['  Source-specific moveout velocities (km/s): ', num2str(best_v_source')]);
    % --- END OF VELOCITY ANALYSIS ---

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if saveresult
        outfilename = [outdir '/' char(source) '.dat'];
        if exist(outfilename,'file')
            unix(['rm ' outfilename]);
        end
        fidout = fopen(outfilename,'a+');
    end
    for n =1:npairs
        np=n;
        synfile=synfilelist(np).name;
        stemp=strsplit(synfile(1:end - length(synfileext)),'.to.');
        src=stemp{1};
        rcv=stemp{2};
        
        if strcmp(src,rcv); continue;end
        stnpair=[src,'_',rcv];
        fntemp=dir([egfdir,'/',source,'/',stnpair,'*',egf_comp,'_N_stack.h5']);
        if isempty(fntemp); continue; end
        %log
        disp(['working on pair: ' stnpair '---> ' num2str(np) '/' num2str(npairs)])
      
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

        %%% source-receiver information
        evla=cdata_substack_N.lat(1);evlo=cdata_substack_N.lon(1);
        stla=cdata_substack_N.lat(2);stlo=cdata_substack_N.lon(2);
        dist=cdata_substack_N.dist;
        dist_ellipse=geo2dist_ellipse(evla,evlo,stla,stlo);
        
        %%%
        syn=synsacdata.DATA1;
        ntsyn=synsacdata.NPTS; dtsyn=synsacdata.DELTA; 
        ttsyn=synsacdata.B:dtsyn:synsacdata.E+dtsyn;

        if ~isempty(find(isnan(syn), 1)), continue, end

        %%% filter EGFs
        
        egffb_pos=nan(length(egf_pos),nfb);egffb_substack_pos=nan(length(egf_pos),nsubstack,nfb);
        egffb_neg=nan(length(egf_pos),nfb);egffb_substack_neg=nan(length(egf_pos),nsubstack,nfb);

        Feff=1/dt/2; % in increment of 2

        for kf=1:nfb
            k=kf;
            fone=fband(k,1)/Feff; ftwo=fband(k,2)/Feff;
            [b,a]=butter(N,[fone ftwo]);
            egffb_pos(:,k)=filtfilt(b,a,egf_pos);
            egffb_neg(:,k)=filtfilt(b,a,egf_neg);
            for im0=1:nsubstack
                im=im0;
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
        for kf=1:nfb
            k=kf;
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
        for kf=1:nfb
            k=kf;
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
        for kf=1:nfb
            k=kf;
              for im0=1:nsubstack
                  im=im0;
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
        tmin0=dist/max(best_v_source);

        tmin=tmin0-2.5*tfmin(1); % arrival time minus the longest period
        if tmin < t_taper; tmin=0.9*t_taper; end
        if tmin > tmaximum; tmin=tmaximum; end
        tmin0minustmin=tmin0-tmin;
        % use only arrivals within the simulation window-400 s to account for finite period, taper and source time function
        if tmin0 < pband(1,1) || tmin0 > tmaximum-tfmin(end)
            continue
        end

        tmax=dist/min(best_v_source) + 2.5*max(tfmin);%tmin0+2.0*tfmin(1);
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
        for kf=1:nfb
            k=kf;
            % Use the specific velocity found for this frequency band and source
            t_center = dist / best_v_source(k);
            
            % Set dynamic window relative to the moveout arrival
            t1(k) = t_center - 2 * tfmin(k) - tmin; % Start window 1.5 periods early
            t2(k) = t_center + 2.5 * tfmin(k) -tmin; % End window 2.5 periods late
            
            % Constraints for valid indexing
            if t1(k) < dt_resample; t1(k) = dt_resample; end
            if t2(k) > (tmax-tmin); t2(k) = tmax-tmin; end
            
            tw_eff1(k) = tmin + t1(k); % These map to your SAC header outputs
            tw_eff2(k) = tmin + t2(k);
            
            it1(k) = round(t1(k) / dt_resample);
            it2(k) = round(t2(k) / dt_resample);
            t1t2len=it2(k)-it1(k)+1;
            
            %
            w=zeros(ntsig,1); %initialize taper of same length
            w(it1(k):it2(k),1)=tukeywin(t1t2len,taperfraction);
            w_anti = max(w) - w; %this is used to filter out the signal window and only use the outside data to compute noise.
            temp_noise_pos = egfsig_pos(:,k).*w_anti;
            temp_noise_neg = egfsig_neg(:,k).*w_anti;
            egfsig_pos(:,k)=egfsig_pos(:,k).*w; % taper the selected time window
            egfsig_neg(:,k)=egfsig_neg(:,k).*w;
            
            synsig(:,k)=synsig(:,k).*w;
            for imm=1:nsubstack
                im=imm;
                egfsubstacksig_pos(:,im,k)=egfsubstacksig_pos(:,im,k).*w;
                egfsubstacksig_neg(:,im,k)=egfsubstacksig_neg(:,im,k).*w;
            end
            
            err_pos=zeros(ntsig,1);err_neg=zeros(ntsig,1);
            for i2=1:ntsig
                i=i2;
                err_pos(i)=std(egfsubstacksig_pos(i,:,k))/sqrt(nsubstack);
                err_neg(i)=std(egfsubstacksig_neg(i,:,k))/sqrt(nsubstack);
            end
            snr_pos(k)=max(abs(egfsig_pos(:,k)))/mean(abs(temp_noise_pos));%mean(abs(egffb_pos(:,k)));
            snr_neg(k)=max(abs(egfsig_neg(:,k)))/mean(abs(temp_noise_neg));%mean(abs(egffb_neg(:,k)));
            %           snr_pos(k)=max(abs(egfsig_pos(:,k)))/max(err_pos);
            %           snr_neg(k)=max(abs(egfsig_neg(:,k)))/max(err_neg);
            
            snr(k)=max([snr_pos(k),snr_neg(k)]);%(snr_pos(k)+snr_neg(k))/2; % average snr from positive and negative lags
        end
        %mean snr for all frequency bands.
        snr_pos_mean=mean(snr_pos);
        snr_neg_mean=mean(snr_neg);
        use_side='m';
        if snr_neg_mean < verylargenumber && snr_pos_mean < verylargenumber
%             if snr_neg_mean<snr_pos_mean
%                 use_side='p';
%             elseif snr_pos_mean<snr_neg_mean 
%                 use_side='n';
%             end
            if snr_neg_mean/snr_pos_mean<2/3
                use_side='p';
            elseif snr_pos_mean/snr_neg_mean < 2/3
                use_side='n';
            else
                use_side='m';
            end
        end
        maxdelay = double(min(max_dV*tmin0,max_dT));
        if maxdelay < 5*dt_resample; maxdelay=5*dt_resample; end % avoid a situation when tmin0 ~ 0
        maxlag=round(maxdelay/dt_resample);
        ttc=(-maxlag:maxlag)*dt_resample;
        xc_pos=[]; xc_neg=[]; xcm_pos=[]; xcm_neg=[];
        c_pos=[]; c_neg=[]; cm_pos=[]; cm_neg=[];
        phase=nan(nfb,1); phaseerr=nan(nfb,1);
        
        rxc=nan(nfb,1); rc_pos=[]; rc_neg=[]; 
        synshift=synsig;
        ec=zeros(nfb,1);
        for kf=1:nfb
            k=kf;
            xc_pos(:,k)=xcorr(egfsig_pos(:,k),synsig(:,k),maxlag);
            xc_neg(:,k)=xcorr(egfsig_neg(:,k),synsig(:,k),maxlag);
            for im0=1:nsubstack
                im=im0;
                xcm_pos(:,im,k)=xcorr(egfsubstacksig_pos(:,im,k),synsig(:,k),maxlag);
                xcm_neg(:,im,k)=xcorr(egfsubstacksig_neg(:,im,k),synsig(:,k),maxlag);
            end
        
            [~, ixcmax_pos]=max(xc_pos(:,k));
            [~, ixcmax_neg]=max(xc_neg(:,k));
%             disp([ttc(ixcmax_pos),ttc(ixcmax_neg)])
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
            ec(k)=(dist-dist_ellipse)/max(best_v_source);  %correction for ellipticity for surface waves with cmax
            phase(k)=phase_temp+ec(k); % with ellipticity correction!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            phasem=nan(nsubstack,1);
            for im0=1:nsubstack
                im=im0;
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
            synerr=model_grid_spacing*110/2/max(best_v_source); % about half of the grid spacing divided by group velocity (~0.05*110 km /2/ 4.5 km/s)
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
                hid=figure('Position',[100 0 1100 800]);
            end
            for kf=1:nfb
                k=kf;
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
                ylabel([num2str(fband(k,1)),' - ',num2str(fband(k,2)),' Hz'])
                % ylabel([num2str(pband(nfb-k+1,1)),' - ',num2str(pband(nfb-k+1,2)),' s'])
                axis([tmin tmax -1.3 1.3]);

                %%%%%%%%%%%%%%%%%% measurements
                subplot(nfb,2,2*k)
                plot(ttc+ec(k),xc_pos(:,k)/max(abs(xc_pos(:,k))),'k-');hold on % for ellipticity correction
                plot(ttc+ec(k),xc_neg(:,k)/max(abs(xc_neg(:,k))),'b--');hold on % 
                plot(phase(k),1,'rv'); hold on
                text(-maxdelay+0.005,0.3,['lag: ' num2str(phase(k),3) '+/-' num2str(phaseerr(k),2)],'FontSize',10);hold on
                text(-maxdelay+0.005,-0.3,['xcoeff: ' num2str(rxc(k),2)],'FontSize',10); hold on
                text(-maxdelay+0.005,-0.9,['snr: ' num2str(snr(k),3)],'FontSize',10);hold on

                % limiting phase < 0.95*maxdelay avoid measurements at 
                % the boundries of the time window,which are likely unreliable
                if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
                        abs(phase(k)) <= min_wavelength*maxdelay && tmin0 >= tfmin(k) 
                    text(+0.011,-0.1,'thumbs UP','Color',[0 0 1]);hold on
                else
                    text(+0.011,-0.1,'thumbs DOWN','Color',[1 0 0]); hold on
                end
                if k==nfb
                    xlabel('Phase difference (s)','FontSize',12);
                end
                axis([ttc(1) ttc(end) -1.3 1.3]);

            end

            sgtitle([strrep(stnpair,'_','-'),': ',num2str(1000*dist,'%.1f'),' m: side = ',use_side]);
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
             for kf=1:nfb
                 k=kf;
                 xcid=[char(src) '/bp' num2str(fband(k,1)) '_' num2str(fband(k,2)) '/' stnpair '_BHZ.P2.CORR.T1T2.SAC'];
                 fbid=['f' num2str(k)];
                 %%%%%% Normal stations %%%%%%%%%%%%%%%%%%%%%
                 if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
                        abs(phase(k)) <= min_wavelength*maxdelay && tmin0 >= tfmin(k)
                     
                     fprintf(fidout,'%s %6.3f %2.0f %s %7.2f %7.2f %s %5.2f %5.2f %5.1f\n',...
                         xcid,phase(k),1,'RL',tw_eff1(k),tw_eff2(k),fbid,phaseerr(k),rxc(k),snr(k)); 
                 end
             end
        end
    end
    if saveresult;fclose(fidout);end
end
