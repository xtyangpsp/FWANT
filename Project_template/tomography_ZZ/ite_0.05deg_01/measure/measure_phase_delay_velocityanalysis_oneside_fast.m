% script to measure phase delays between the observed and synthetic waveforms
% Modifications by Xiaotao Yang:
% 1. changed to use moveout-awared time window, with specified search range. - 2/20/2026
% 2. removed option to shift the EGFs, leave this to data preprocessing. - 2/20/2026
% 3. Added SNR-weighted side combining logic - 2/21/2026.
%    - Uses a bootstrap approach: Grid search finds the best velocity.
%    - A refined frequency-dependent window and adaptive noise window (3-cycle min)
%      calculate SNR for each side.
%    - Sides are combined using SNR-squared weighting to maximize signal quality.
%    - The same weights are applied to substacks for consistent error estimation.
% 4. Made the max cross-correlation lag (cycle-skip guard) period-specific,
%    instead of a single lag bound shared by all frequency bands. - 10/02/2026
% 5. Speed/memory version (results unchanged): the EGF stack files are listed
%    once per source instead of two dir() calls per pair (twice); every EGF
%    is filtered/tapered/normalised once per band in the velocity analysis
%    instead of once per test velocity; big per-pair arrays are released;
%    per-source timing is printed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Set up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise/src'))

%%%%%%%%%%%%%%% Parameter section 1: users should modify for each project
%%%%%%%%%%%%%%% and iteration
PROJHOME = '/depot/xtyang/data/projects/xtyang/RioGrande/FWANT_ZZ';
%PROJHOME = '/Volumes/xtyang/data/projects/xtyang/RioGrande/FWANT_TT';
ite = '/ite_0.025deg_01';
wkdir = [PROJHOME ite];
egfdir = [PROJHOME '/data/data_riogrande/PAIRS_TWOSIDES_gaussian_a1t3'];
syndir = [wkdir '/syn.seismograms'];
outdir = [wkdir '/measure/measure_stnpair'];
plotdir = strcat(wkdir,'/measure/plots_test');
stainfo = [PROJHOME '/STinfo/station_riogrande_withdata_formatok.txt'];

%wave types to process.

synfileext='.fz.Uz.SAC'; %check the displacement file extension in the syn.seismograms folder.
egf_comp='ZZ'; %TT (Love waves) or ZZ (Rayleigh waves)

pband=[7.5 15; 10 20; 15 30; 20 40; 30 55; 40 70];
fband = flip(1./pband,2);

% max_dV and max_dT are combined in checking the phase delays to make sure
% they are within the allowed perturbation range and to avoid
% cycle-skipping
max_dV = 0.2;
max_dT = 10;

% Additional cycle-skip guard: cap the max cross-correlation lag (and the
% QC delay-acceptance bound) at this fraction of each band's own mean
% period, in addition to max_dV/max_dT above. This keeps the lag search
% range from reaching into a neighboring-cycle sidelobe of the
% correlation function, which is the mechanism behind cycle skipping and
% is otherwise most likely to bite the longer-period bands (their mean
% period is comparable to, or smaller than, a fixed max_dT).
maxlag_period_fraction = 0.5;

% QC to save phase delays.
snr_cutoff = 5;
xcoeff_cutoff = 0.7;
min_wavelength = 1.; %
maxdelay_scaling_buffer=0.99; %measured delay needs to be within this scaling for maxdelay.
min_substack = 3;

% dt in seconds.
dt_resample=0.1;
dt_egf = 0.2;

tmaximum=800.0;   %max length of the synthetic EGFs
tmaximumegf=1000.0; % length of the observed EGFs

d_taper_period_scale = 0.5; % taper length in time scales with the dominant period.

v_search_grid = 2.:0.1:5.5; %velocity range for velocity analysis to decide the best signal window
cmax_default=linspace(3,5,size(pband,1)); %if velocity analysis failed, fall back to these default velocities to search for the phase delays.
search_window_dist_scaling = 0.17; %increasing window width with distance to account for the spreadout of wavetrain.

syn_shaping_tshift = 3; %time shift in shaping wavelet, same as in synthetics.
% Scaling factors for each frequency band (one for each row in fband)
% Increase these for higher frequencies to capture scattering/coda
win_scales = [5, 5, 2.5, 2, 1.5, 1];
model_grid_spacing=0.025;

%plot control
fig_flag = 1;  %turn off for parfor.
save_fig=0;

%save the measured phase delay or not.
save_delays=1;

%speed / diagnostics
print_timing = 1;         % print timing and client memory per source

%%%%%%%%%%% Parameter section 2: derived parameters or those that users
%%%%%%%%%%% don't normally need to change.
d_taper_time=zeros(size(pband,1),1);
d_taper_samples_egf=zeros(size(pband,1),1);
taperfraction_egf=zeros(size(pband,1),1);
taperfraction_syn=zeros(size(pband,1),1);
for i=1:size(pband,1)
    d_taper_time(i) = d_taper_period_scale*mean(pband(i,:));
    d_taper_samples_egf(i) = round(d_taper_time(i)/dt_egf);

    taperfraction_egf(i)=d_taper_time(i)/tmaximumegf;
    taperfraction_syn(i)=d_taper_time(i)/tmaximum;
end

verylargenumber = 1.e9;

N=2; %number of poles in butterworth filter design

if strcmp(egf_comp,'TT')
	delay_label = 'LV';
	chan_label = 'BHT';
elseif strcmp(egf_comp,'ZZ')
	delay_label = 'RL';
	chan_label = 'BHZ';
else
	error('wrong egf_comp, has to be TT or ZZ')
end

tfmin=1./fband(:,1);
[nfb, nc]=size(fband);

%read in station information.
fid = fopen(stainfo);
recv = textscan(fid,'%s%s%*f%*f%*f');
fclose(fid);
nrecv=length(recv{1});
ntwk=recv{1}; stnm=recv{2};
sourcelist=cell(nrecv,1);
for nstsrc=1:nrecv
    sourcelist{nstsrc}=[char(ntwk(nstsrc)) '.' char(stnm(nstsrc))];
end
nsource=length(sourcelist);

if ~exist(outdir,'dir'); system(['mkdir ' outdir]); end
if ~exist(plotdir,'dir'); system(['mkdir ' plotdir]); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:nsource %change the range here to process a subset of sources
    t_source = tic;
    source = sourcelist{ii};
    disp([num2str(ii),' --> ',source]);
    synfilelist=dir([syndir,'/',source,'/*',synfileext]);
    npairs=length(synfilelist);

    disp('    Performing velocity analysis ... ');
    all_egfs = cell(npairs, 1);
    all_dists = zeros(npairs, 1);
    valid_pair_idx = [];

    % list the EGF stack files of this source ONCE; file name of each pair
    % (empty if the pair has no N or P stack), used here and in the parfor
    dN = dir([egfdir, '/', source, '/*', egf_comp, '_N_stack.h5']);
    dP = dir([egfdir, '/', source, '/*', egf_comp, '_P_stack.h5']);
    keyN = cell(numel(dN),1); keyP = cell(numel(dP),1);
    for q = 1:numel(dN)
        t = strsplit(dN(q).name, '_');
        if numel(t) >= 2; keyN{q} = [t{1}, '_', t{2}]; end
    end
    for q = 1:numel(dP)
        t = strsplit(dP(q).name, '_');
        if numel(t) >= 2; keyP{q} = [t{1}, '_', t{2}]; end
    end
    egfname_N = cell(npairs, 1); egfname_P = cell(npairs, 1);
    for n = 1:npairs
        stemp = strsplit(synfilelist(n).name(1:end - length(synfileext)), '.to.');
        if strcmp(stemp{1}, stemp{2}); continue; end
        iN = find(strcmp(keyN, [stemp{1}, '_', stemp{2}]), 1);
        iP = find(strcmp(keyP, [stemp{1}, '_', stemp{2}]), 1);
        if ~isempty(iN) && ~isempty(iP)
            egfname_N{n} = dN(iN).name;
            egfname_P{n} = dP(iP).name;
        end
    end

    for n = 1:npairs
        synfile = synfilelist(n).name;
        stemp = strsplit(synfile(1:end - length(synfileext)), '.to.');
        src = stemp{1}; rcv = stemp{2};
        if strcmp(src, rcv); continue; end
        stnpair = [src, '_', rcv];

        if ~isempty(egfname_N{n}) && ~isempty(egfname_P{n})
            d_N = extract_corrdata_asdf([egfdir, '/', source, '/', egfname_N{n}]);
            d_P = extract_corrdata_asdf([egfdir, '/', source, '/', egfname_P{n}]);
            if isempty(d_N.value{1}.(egf_comp).data) || isempty(d_P.value{1}.(egf_comp).data)
                disp(['data is empty, skipped! Pair: ',stnpair])
                continue
            end
            % Symmetric stack to enhance the surface wave signal for moveout search
            combined_data = (d_P.value{1}.(egf_comp).data + d_N.value{1}.(egf_comp).data) / 2;
            all_egfs{n} = combined_data;
            all_dists(n) = d_P.value{1}.(egf_comp).dist;
            valid_pair_idx = [valid_pair_idx; n];
        end
    end

    all_egfs_subset = all_egfs(valid_pair_idx);
    all_dists_subset = all_dists(valid_pair_idx);
    npairs_valid = length(valid_pair_idx);
    if npairs_valid < 1; continue; end
    npts_egf = length(all_egfs_subset{1});
    if npts_egf - 1 > round(tmaximumegf/dt_egf)
        disp('npts of egf is inconsistent with the specified length.')
    end
    taxis_egf = (0:npts_egf-1)*dt_egf;


    % A velocity (moveout) analysis needs more than one pair. With a single
    % valid pair, skip it and use the default group velocities for all bands
    % (same behavior as the original script).
    if npairs_valid == 1
        best_v_source = cmax_default(:);
        vel_bands = [];
        disp('Only one pair available. Use default group velocities for all frequencies.')
    else
        best_v_source = zeros(nfb, 1);
        vel_bands = 1:nfb;
        if fig_flag
            fig_1=figure('Position',[400 400 1000 650]);
            fig_2=figure('Position',[400 400 1000 650]);
        end
    end
    for k = vel_bands
        w_taper_egf = tukeywin(npts_egf,taperfraction_egf(k));
        energy_results = zeros(length(v_search_grid), 1);
        [b, a] = butter(N, [fband(k,1) fband(k,2)] / (1/dt_egf/2));
        % filter, taper and normalise every EGF ONCE for this band; the
        % result does not depend on the test velocity
        filt_egf = cell(npairs_valid, 1);
        for idx = 1:npairs_valid
            if size(all_egfs_subset{idx},1) > npts_egf; continue; end
            temp_data = filtfilt(b, a, all_egfs_subset{idx}) .* w_taper_egf;
            filt_egf{idx} = temp_data/max(abs(temp_data));
        end
        for iv = 1:length(v_search_grid)
            v_test = v_search_grid(iv);
            stack_energy = 0;
            for idx = 1:npairs_valid
                if size(all_egfs_subset{idx},1) > npts_egf; continue; end
                temp_data = filt_egf{idx};
                t_arrival = syn_shaping_tshift + all_dists_subset(idx) / v_test;
                it_arrival = round(t_arrival / dt_egf);
                win_samples = round((all_dists_subset(idx)*search_window_dist_scaling + 0.5*win_scales(k)*tfmin(k)) / dt_egf);
                if it_arrival > win_samples && (it_arrival + win_samples) < length(temp_data)
                    env = abs(hilbert(temp_data(it_arrival-win_samples : it_arrival+win_samples)));
                    stack_energy = stack_energy + sum(env);
                end
            end
            energy_results(iv) = stack_energy;
        end
        [~, max_idx] = max(energy_results);
        if max(energy_results) == 0
            disp(['1st round velocity analysis returned 0 energy. Fall back to default velocity for source ',source, ' at frequency index ',num2str(k)])
            best_v_alldist = cmax_default(k);
        else
            best_v_alldist = v_search_grid(max_idx);
        end

        %second search to redo it after excluding the near-source stations,
        %defined by min_wavelength and the 1st round velocity.
        %in this search, we also QC to only consider SNR>snr_cutoff
        %signal is computed based on the velocity anlaysis result from
        %round one

        for iv = 1:length(v_search_grid)
            v_test = v_search_grid(iv);
            stack_energy = 0;
            for idx = 1:npairs_valid
                if size(all_egfs_subset{idx},1) > npts_egf; continue; end
                min_dist = min_wavelength * tfmin(k) * best_v_alldist;
                if all_dists_subset(idx) < min_dist; continue; end %skip near source stations.

                temp_data = filt_egf{idx};

                t_arrival = syn_shaping_tshift + all_dists_subset(idx) / v_test;
                it_arrival = round(t_arrival / dt_egf);
                win_samples = round((all_dists_subset(idx)*search_window_dist_scaling + 0.5*win_scales(k)*tfmin(k)) / dt_egf);
                if it_arrival > win_samples && (it_arrival + win_samples) < length(temp_data)
                    %compute SNR
                    temp_signal = temp_data(it_arrival-win_samples+1 : it_arrival+win_samples);
                    temp_noise = [temp_data(1:max(2,it_arrival - win_samples));...
                        temp_data(min(it_arrival + win_samples + 1,npts_egf-win_samples):npts_egf)];

                    temp_snr = rms(abs(temp_signal))/rms(abs(temp_noise));
                    if temp_snr < snr_cutoff; continue;end
                    env = abs(hilbert(temp_signal));
                    stack_energy = stack_energy + sum(env);
                end
            end
            energy_results(iv) = stack_energy;
        end
        [~, max_idx] = max(energy_results);
        if max(energy_results) == 0
            disp(['2nd round velocity analysis returned 0 energy. Fall back to default velocity for source ',source, ' at frequency index ',num2str(k)])
            best_v_source(k) = cmax_default(k);
        else
            best_v_source(k) = v_search_grid(max_idx);
        end
        %plot moveout for each frequency
        if fig_flag
            %energy curves.
            figure(fig_1);
            subplot(2,ceil(nfb/2),k); hold on;
            plot(v_search_grid,energy_results)
            hold on
            plot([best_v_source(k) best_v_source(k)],[min(energy_results) max(energy_results)],'r')
            hold off
            xlabel('test velocities (km/s)')
            ylabel('stacked amplitudes')
            box on;
            title([source,' at ',num2str(pband(k,1)),'-',num2str(pband(k,2)),' s: ',num2str(best_v_source(k)),' km/s'])
            if save_fig
                saveas(gcf, fullfile(plotdir, [source, '_moveout_energy.png']))
            end
            %moveout wiggles
            figure(fig_2)
            subplot(2,ceil(nfb/2),k); hold on;
            wiggle_scale=0.01*range(all_dists_subset);
            for idx = 1:npairs_valid
                temp_data = filt_egf{idx};
                if isempty(temp_data); continue; end
                plot(taxis_egf,all_dists_subset(idx)+wiggle_scale*temp_data,'k')
            end
            %plot the best velocity line
            plot(taxis_egf + syn_shaping_tshift,best_v_source(k)*taxis_egf,'r')
            ylim([min(all_dists_subset)-wiggle_scale max(all_dists_subset)+wiggle_scale])
            box on;
            hold off
            title([source,' at ',num2str(pband(k,1)),'-',num2str(pband(k,2)),' s: ',num2str(best_v_source(k)),' km/s'])
            xlim([0 min([2.5*max(all_dists_subset)/best_v_source(k),max(taxis_egf)])])
            xlabel('time (s)')
            ylabel('distance (km)')
            set(gca,'TickDir','out')

            if save_fig
                saveas(gcf, fullfile(plotdir, [source, '_moveout_velocity.png']))
            end
        end
    end
    disp(['Best velocity for source ',source])
    best_v_source
    % result_strings collects output to avoid FID conflicts
    % error('test')
    all_egfs = []; all_egfs_subset = []; filt_egf = [];   % free client memory before the parfor
    t_vel_s = toc(t_source);
    t_par = tic;
    result_strings = cell(npairs, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% PARFOR for parallel here%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor n = 1:npairs

        best_v_pair=zeros(nfb,1); %store refined velocity for each pair for each frequency.

        synfile = synfilelist(n).name;
        stemp = strsplit(synfile(1:end - length(synfileext)), '.to.');
        src = stemp{1}; rcv = stemp{2};
        if strcmp(src, rcv); continue; end
        stnpair = [src, '_', rcv];

        if isempty(egfname_N{n}) || isempty(egfname_P{n}); continue; end
        f_N_stack = struct('name', egfname_N{n});
        f_P_stack = struct('name', egfname_P{n});

        disp(['working on pair: ' stnpair ' (' num2str(n) '/' num2str(npairs) ')'])

        d_N = extract_corrdata_asdf([egfdir,'/',source,'/',f_N_stack.name]);
        d_P = extract_corrdata_asdf([egfdir,'/',source,'/',f_P_stack.name]);
        d_N_sub = extract_corrdata_asdf([egfdir,'/',source,'/',strrep(f_N_stack.name,'_stack','')]);
        d_P_sub = extract_corrdata_asdf([egfdir,'/',source,'/',strrep(f_P_stack.name,'_stack','')]);

        npts_egf = length(d_P.value{1}.(egf_comp).data);
        taxis_egf = (0:npts_egf-1)*dt_egf;
        dist = d_P.value{1}.(egf_comp).dist;
        nsubstack = length(d_P_sub.value{1}.(egf_comp).time);
        if nsubstack < min_substack; continue; end

        synsacdata = readsac([syndir,'/',source,'/',synfile]);
        syn_raw = synsacdata.DATA1 - mean(synsacdata.DATA1);

        ttuni = 0:dt_resample:tmaximum;
        ntuni = length(ttuni);
        egffb_weighted = zeros(ntuni, nfb);
        egffb_sub_weighted = zeros(ntuni, nsubstack, nfb);
        synfb = zeros(ntuni, nfb);

        snr_pos = zeros(nfb, 1); snr_neg = zeros(nfb, 1);
        t1_band = zeros(nfb, 1); t2_band = zeros(nfb, 1);
        snr=nan(nfb,1);
        for k = 1:nfb
            % Pre-filter buffers
            w_taper_egf = tukeywin(npts_egf, taperfraction_egf(k));
            egf_pos = d_P.value{1}.(egf_comp).data .* w_taper_egf;
            egf_neg = d_N.value{1}.(egf_comp).data .* w_taper_egf;
            egf_sub_pos = d_P_sub.value{1}.(egf_comp).data .* w_taper_egf;
            egf_sub_neg = d_N_sub.value{1}.(egf_comp).data .* w_taper_egf;

            [b, a] = butter(N, [fband(k,1) fband(k,2)] / (1/dt_egf/2));
            f_pos = filtfilt(b, a, egf_pos);
            f_neg = filtfilt(b, a, egf_neg);

            %normalize the data
            f_pos = f_pos/max(abs(f_pos));
            f_neg = f_neg/max(abs(f_neg));
            % --- NEW REFINEMENT LOGIC START ---
            % 1. Moveout-aware "First Guess"
            t_moveout = syn_shaping_tshift + dist / best_v_source(k);

            % 2. Define a search range around the moveout (+/- 0.5*scaling cycles)
            % to find the actual observed arrival peak.
            search_half_width = dist*search_window_dist_scaling + win_scales(k) * tfmin(k); %narrower scaling for refinement search of the EGF arrivals.
            it_search = max(d_taper_samples_egf(k),round((t_moveout - search_half_width)/dt_egf)) : ...
                        min(npts_egf,round((t_moveout + search_half_width)/dt_egf));
%             it_noise = npts_egf - length(it_search): npts_egf;
%             it_noise = it_noise - d_taper_samples;
            it_noise = [d_taper_samples_egf(k):max(d_taper_samples_egf(k) + 1,it_search(1)-1),...
                it_search(end) + 1: npts_egf - d_taper_samples_egf(k)];
            % 3. Find the peak of the envelope on the combined observed data
            % Using the weighted combined signal helps stability.
            % (Note: Pre-calculating combined_f for search here)
            snr_p_tmp = max(abs(f_pos(it_search))) / rms(abs(f_pos(it_noise))); % Rough SNR for peak search
            snr_n_tmp = max(abs(f_neg(it_search))) / rms(abs(f_neg(it_noise)));
            w_p_tmp = snr_p_tmp^2 / (snr_p_tmp^2 + snr_n_tmp^2 + eps);
            w_n_tmp = snr_n_tmp^2 / (snr_p_tmp^2 + snr_n_tmp^2 + eps);
            combined_search = w_p_tmp * f_pos + w_n_tmp * f_neg;

            [~, local_idx] = max(combined_search(it_search));
            if isempty(local_idx)
                t_refined_center = t_moveout; % Fallback
            else
                t_refined_center = taxis_egf(it_search(local_idx));
                v_candidate = dist/(t_refined_center - syn_shaping_tshift);
                abs_tol = min(max_dT, 0.5*mean(pband(k,:))); % never drift more than ~half a period, or max_dT, from the prediction
                if v_candidate <= 0 ...
                        || abs(1/v_candidate - 1/best_v_source(k)) > max_dV/best_v_source(k) ...
                        || abs(t_refined_center - t_moveout) > abs_tol
                    t_refined_center = t_moveout; % reject implausible pick, fall back to moveout prediction
                end
            end

            %retrieve best pair velocity from the refined center.
            best_v_pair(k) = dist/(t_refined_center - syn_shaping_tshift);

            % 4. Center the measurement window on the REFINED arrival
            % Focus on the main packet to avoid scattering/coda decorrelation.
            t1_band(k) = max(dt_resample, t_refined_center - search_half_width);
            t2_band(k) = min(tmaximumegf, t_refined_center + search_half_width);
            % --- NEW REFINEMENT LOGIC END ---
            % Refined window for SNR and measurement
            it_sig = max(d_taper_samples_egf(k),round(t1_band(k)/dt_egf)):min(npts_egf-d_taper_samples_egf(k),round(t2_band(k)/dt_egf));

            % Adaptive Noise Window
            it_noise = [d_taper_samples_egf(k):max(d_taper_samples_egf(k) + 1,it_sig(1)-1),...
                it_sig(end) + 1: npts_egf - d_taper_samples_egf(k)];

            snr_p = max(abs(f_pos(it_sig))) / rms(abs(f_pos(it_noise)));
            snr_n = max(abs(f_neg(it_sig))) / rms(abs(f_neg(it_noise)));
            snr_pos(k) = snr_p; snr_neg(k) = snr_n;

            % Weighted Combination (SNR Squared)
            w_p = snr_p^2 / (snr_p^2 + snr_n^2);
            w_n = snr_n^2 / (snr_p^2 + snr_n^2);
            combined_f = w_p * f_pos + w_n * f_neg;

            % Apply to stacks and resample
            egffb_weighted(:,k) = interp1(taxis_egf, combined_f, ttuni, 'linear', 0);

            [b_s, a_s] = butter(N, [fband(k,1) fband(k,2)] / (1/synsacdata.DELTA/2));
            syn = syn_raw .* tukeywin(length(syn_raw), taperfraction_syn(k));
            f_syn = filtfilt(b_s, a_s, syn);
            synfb(:,k) = interp1(synsacdata.B + (0:synsacdata.NPTS-1)*synsacdata.DELTA, f_syn, ttuni, 'linear', 0);

            % Apply weights to substacks
            for im = 1:nsubstack
                f_sub_p = filtfilt(b, a, egf_sub_pos(:,im));
                f_sub_n = filtfilt(b, a, egf_sub_neg(:,im));

                snr_p_sub = max(abs(f_sub_p(it_sig))) / rms(abs(f_sub_p(it_noise)));
                snr_n_sub = max(abs(f_sub_n(it_sig))) / rms(abs(f_sub_n(it_noise)));

                % Weighted Combination (SNR Squared)
                w_p_sub = snr_p_sub^2 / (snr_p_sub^2 + snr_n_sub^2);
                w_n_sub = snr_n_sub^2 / (snr_p_sub^2 + snr_n_sub^2);

                combined_sub = w_p_sub * f_sub_p + w_n_sub * f_sub_n;
                egffb_sub_weighted(:,im,k) = interp1(taxis_egf, combined_sub, ttuni, 'linear', 0);
            end

            % get SNR for combined data.
            snr(k) = rms(abs(combined_f(it_sig))) / rms(abs(combined_f(it_noise)));
        end

        % Prepare segments for measurement
        % Define general shared window
        tmin_arrival_all = dist / max(best_v_pair);
        tmax_arrival_all = dist / min(best_v_pair);
        tmin = max(min(d_taper_time), tmin_arrival_all - max(win_scales.*(tfmin')));
        tmax = min(tmaximum, tmax_arrival_all + max(win_scales.*(tfmin')));

        itmin = round(tmin/dt_resample); itmax = round(tmax/dt_resample);
        egfsig = egffb_weighted(itmin:itmax, :);
        synsig = synfb(itmin:itmax, :);
        egfsubstacksig = egffb_sub_weighted(itmin:itmax, :, :);
        egffb_sub_weighted = []; egf_sub_pos = []; egf_sub_neg = []; d_N_sub = []; d_P_sub = [];   % release memory
        ntsig = size(egfsig, 1);

        % --- PERIOD-SPECIFIC LAG CAP ---
        % maxdelay_band(k)/maxlag_npts_k now replace the single, band-
        % independent maxdelay_allband/maxlag_npts that used to be
        % computed once per pair and shared by every frequency band. A
        % lag range that is fine for a short-period band can be wide
        % enough, relative to a longer-period band's own period, to let
        % the correlation search lock onto a neighboring-cycle sidelobe
        % (cycle skipping). Capping each band's lag at a fraction of its
        % own mean period removes that ambiguity, independent of the
        % max_dV/max_dT bounds (which are still applied and generally
        % dominate at short distances).
        phase = nan(nfb,1); phaseerr = nan(nfb,1);
        rxc = nan(nfb,1);
        synshift=synsig;
        ec=zeros(nfb,1);
        maxdelay_band = nan(nfb,1);
        xc = cell(nfb,1);
        ttc = cell(nfb,1);
        current_pair_string="";
        for k = 1:nfb
            maxdelay_band(k) = min([max_dV * tmax_arrival_all, max_dT, maxlag_period_fraction*mean(pband(k,:))]);
            maxlag_npts_k = round(maxdelay_band(k)/dt_resample);
            ttc{k} = (-maxlag_npts_k:maxlag_npts_k)*dt_resample;

            t1_rel = t1_band(k) - tmin;
            t2_rel = t2_band(k) - tmin;
            it1 = max(1, round(t1_rel/dt_resample));
            it2 = min(ntsig, round(t2_rel/dt_resample));

            w_win = zeros(ntsig, 1);
            win_len_sec = max(t2_band(k) - t1_band(k), eps);
            taperfraction_win = min(1, d_taper_time(k) / win_len_sec);
            w_win(it1:it2) = tukeywin(it2-it1+1, taperfraction_win);

            d_win = egfsig(:,k) .* w_win;
            s_win = synsig(:,k) .* w_win;

            [xc{k}, lags] = xcorr(d_win, s_win, maxlag_npts_k);
            [~, i_max] = max(xc{k});

            % Correction for ellipticity
            ec(k) = (dist - geo2dist_ellipse(d_N.value{1}.(egf_comp).lat(1),d_N.value{1}.(egf_comp).lon(1),...
                d_N.value{1}.(egf_comp).lat(2),d_N.value{1}.(egf_comp).lon(2))) / best_v_pair(k);
            phase(k) = lags(i_max) * dt_resample + ec(k);

            % Substacks for error
            phasem = zeros(nsubstack, 1);
            for im = 1:nsubstack
                sub_win = egfsubstacksig(:,im,k) .* w_win;
                [xc_sub, lags_sub] = xcorr(sub_win, s_win, maxlag_npts_k);
                [~, i_sub] = max(xc_sub);
                phasem(im) = lags_sub(i_sub) * dt_resample + ec(k);
            end

            % about half of the grid spacing divided by the correponding best group velocity
            synerr = model_grid_spacing*110/2/best_v_source(k);
            %error of the mean + sampling error for obs & syn
            %%% add the errors due to synthetic waveforms
            %%% inspection of reciprical synthetic waveforms show ~0.3 s shift, caused by the fact
            %%% that the source & receiver are not exactly on grids and interpolation is needed.
            %%% assuming data errors and synthetic errors are independent
            phaseerr(k) = sqrt((std(phasem)/sqrt(nsubstack) + dt_resample)^2 + synerr^2);

            % cross-correlation coefficient
            itshift=round(phase(k)/dt_resample);
            if itshift > 0
                synshift(1:itshift,k)=0;
                synshift(itshift+1:ntsig,k)=synsig(1:ntsig-itshift,k);
            elseif itshift < 0
                synshift(ntsig+itshift:ntsig,k)=0;
                synshift(1:ntsig+itshift,k)=synsig(-itshift+1:ntsig,k);
            end
            rxc_temp=corrcoef(egfsig(:,k),synshift(:,k));
            rxc(k)=rxc_temp(1,2);

            %form phase delay string:
            if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
                       abs(phase(k)) <= maxdelay_scaling_buffer*maxdelay_band(k) && ...
                       dist >= min_wavelength*best_v_source(k)*tfmin(k)
                xcid=[char(src) '/bp' num2str(fband(k,1)) '_' num2str(fband(k,2)) '/' stnpair '_' chan_label '.P2.CORR.T1T2.SAC'];
                outstring=sprintf('%s %8.4f %2.0f %s %8.4f %8.4f %s %8.4f %8.4f %5.1f\n',...
                    xcid, phase(k), 1, delay_label, tmin+t1_band(k), tmin+t2_band(k), ['f' num2str(k)], phaseerr(k), rxc(k), snr(k));
                current_pair_string = current_pair_string + outstring;
            end
        end

        result_strings{n} = current_pair_string;

        if fig_flag
            figure('Position',[100 0 1100 800]);
            for k = 1:nfb
                subplot(nfb, 2, 2*k-1); hold on;
                plot(ttuni(itmin:itmax), synshift(:,k)/max(abs(synshift(:,k))), 'r', 'LineWidth', 1.2);
                plot(ttuni, egffb_weighted(:,k)/max(abs(egffb_weighted(:,k))), 'k');
                line([t1_band(k) t1_band(k)], [-1 1], 'Color', 'm', 'LineStyle', '--');
                line([t2_band(k) t2_band(k)], [-1 1], 'Color', 'm', 'LineStyle', '--');
                xlim([tmin tmax]); ylabel([num2str(pband(k,1)) '-' num2str(pband(k,2)) 's']);
                if k==1; legend('Syn','Combined EGF'); end
                box on;
                axis on;

                subplot(nfb, 2, 2*k);
                plot(ttc{k}+ec(k),xc{k}/max(abs(xc{k})),'k-');hold on % for ellipticity correction
                plot(phase(k),1,'rv'); hold on
                text(-maxdelay_band(k)+0.005,0.3,['lag: ' num2str(phase(k),3) '+/-' num2str(phaseerr(k),2)],'FontSize',10);hold on
                text(-maxdelay_band(k)+0.005,-0.3,['xcoeff: ' num2str(rxc(k),2)],'FontSize',10); hold on
                text(-maxdelay_band(k)+0.005,-0.9,['snr: ' num2str(snr(k),3)],'FontSize',10);hold on

                if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && snr(k) < verylargenumber && ...
                       abs(phase(k)) <= maxdelay_scaling_buffer*maxdelay_band(k) && ...
                       dist >= min_wavelength*best_v_source(k)*tfmin(k)
                    text(+0.011,-0.1,'thumbs UP','Color',[0 0 1]);hold on
                else
                    text(+0.011,-0.1,'thumbs DOWN','Color',[1 0 0]); hold on
                end
                if k==nfb
                    xlabel('Phase difference (s)','FontSize',12);
                end
                axis([ttc{k}(1) ttc{k}(end) -1.3 1.3]);
            end
            sgtitle([strrep(stnpair,'_','-') ' | One-sided | ',num2str(dist),' km']);
            if save_fig
                saveas(gcf, fullfile(plotdir, [stnpair '.png']))
            elseif isempty(getCurrentTask()) %if not in parfor worker.
                pause
            end
            close all;
        end
    end
    t_par_s = toc(t_par);
    if save_delays
        outfilename = [outdir '/' char(source) '.dat'];
        fidout = fopen(outfilename, 'w');
        for n = 1:npairs
            if ~isempty(result_strings{n})
                fprintf(fidout, '%s', result_strings{n});
            end
        end
        fclose(fidout);
    end
    if fig_flag && npairs_valid > 1
        if ishandle(fig_1); close(fig_1); end
        if ishandle(fig_2); close(fig_2); end
    end
    if print_timing
        rss_mb = NaN;
        try
            [st_ps, out_ps] = system(sprintf('ps -o rss= -p %d', feature('getpid')));
            if st_ps == 0; rss_mb = str2double(out_ps)/1024; end
        catch
        end
        fprintf('TIMING source %s (%d/%d): pairs %d, valid %d | velocity analysis %.1f s | parfor %.1f s | total %.1f s | client RSS %.0f MB\n', ...
            source, ii, nsource, npairs, npairs_valid, t_vel_s, t_par_s, toc(t_source), rss_mb);
    end
end
