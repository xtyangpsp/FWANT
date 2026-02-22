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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Set up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise/src'))

%%%%%%%%%%%%%%% Parameter section 1: users should modify for each project
%%%%%%%%%%%%%%% and iteration 
PROJHOME = '/Users/xtyang/Work/Research/Projects/FloridanAquifer/RiverRise';
ite = '/ite_0.00001deg_01';
wkdir = [PROJHOME ite];
egfdir = [PROJHOME '/data/data_riverrise/PAIRS_TWOSIDES_gaussian_a0.005t0.015'];
syndir = [wkdir '/syn.seismograms'];
outdir = [wkdir '/measure/measure_stnpair'];
plotdir = strcat(wkdir,'/measure/plots_test');
stainfo = [PROJHOME '/STinfo/station_riverrise_withdata_formatok.txt'];

fband=[4 8; 6 10; 7 12; 9 15; 12 20; 15 25];
pband = flip(1./fband,2);

% max_dV and max_dT are combined in checking the phase delays to make sure
% they are within the allowed perturbation range and to avoid
% cycle-skipping
max_dV = 0.25;
max_dT = 0.3; 

% QC to save phase delays.
snr_cutoff = 4;
xcoeff_cutoff = 0.6; 
min_wavelength = 0.75; %
min_substack = 3; 

% dt in seconds.
dt_resample=0.002;
dt_egf = 0.01;
d_taper_samples=10; %number of observed EGF samples to taper
d_taper_time=10*max(dt_egf,dt_resample);
tmaximum=4.0;   %max length of the synthetic EGFs
tmaximumegf=5.0; % length of the observed EGFs

v_search_grid = 0.2:0.01:2.0; %velocity range for velocity analysis to decide the best signal window
model_grid_spacing=0.00001; 

%plot control
fig_flag = 1;  
save_fig=0;

%save the measured phase delay or not.
save_delays=1; 

%%%%%%%%%%% Parameter section 2: derived parameters or those that users
%%%%%%%%%%% don't normally need to change.
verylargenumber = 1.e9; 
taperfraction=d_taper_time/tmaximumegf; 
syntaperfraction=d_taper_time/tmaximum; 
N=2; %number of poles in butterworth filter design

synfileext='.fz.Uz.SAC';
egf_comp='ZZ';
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

parfor ii = 1:nsource
    source = sourcelist{ii};
    disp([num2str(ii),' --> ',source]);
    synfilelist=dir([syndir,'/',source,'/*',synfileext]);
    npairs=length(synfilelist);

    disp('    Performing velocity analysis ... ');
    all_egfs = cell(npairs, 1);
    all_dists = zeros(npairs, 1);
    valid_pair_idx = [];

    for n = 1:npairs
        synfile = synfilelist(n).name;
        stemp = strsplit(synfile(1:end - length(synfileext)), '.to.');
        src = stemp{1}; rcv = stemp{2};
        if strcmp(src, rcv); continue; end
        stnpair = [src, '_', rcv];
        fntemp_N = dir([egfdir, '/', source, '/', stnpair, '*', egf_comp, '_N_stack.h5']);
        fntemp_P = dir([egfdir, '/', source, '/', stnpair, '*', egf_comp, '_P_stack.h5']);
        
        if ~isempty(fntemp_N) && ~isempty(fntemp_P)
            d_N = extract_corrdata_asdf([egfdir, '/', source, '/', fntemp_N.name]);
            d_P = extract_corrdata_asdf([egfdir, '/', source, '/', fntemp_P.name]);
            % Symmetric stack to enhance the surface wave signal for moveout search
            combined_data = (d_P.value{1}.(egf_comp).data + flipud(d_N.value{1}.(egf_comp).data)) / 2;
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
    taxis_egf = (0:npts_egf-1)*dt_egf;
    w_taper_egf = tukeywin(npts_egf,taperfraction);
    
    %
    if fig_flag
        figure('Position',[400 400 1000 650]);
    end
    best_v_source = zeros(nfb, 1);
    for k = 1:nfb
        energy_results = zeros(length(v_search_grid), 1);
        [b, a] = butter(N, [fband(k,1) fband(k,2)] / (1/dt_egf/2));
        for iv = 1:length(v_search_grid)
            v_test = v_search_grid(iv);
            stack_energy = 0;
            for idx = 1:npairs_valid
                temp_data = filtfilt(b, a, all_egfs_subset{idx}) .* w_taper_egf;
                t_arrival = all_dists_subset(idx) / v_test;
                it_arrival = round(t_arrival / dt_egf);
                win_samples = round(tfmin(k) / dt_egf);
                if it_arrival > win_samples && (it_arrival + win_samples) < length(temp_data)
                    env = abs(hilbert(temp_data(it_arrival-win_samples : it_arrival+win_samples)));
                    stack_energy = stack_energy + max(env);
                end
            end
            energy_results(iv) = stack_energy;
        end
        [~, max_idx] = max(energy_results);
        best_v_source(k) = v_search_grid(max_idx);

        %plot moveout for each frequency
        if fig_flag
            subplot(2,ceil(nfb/2),k); hold on;
            wiggle_scale=0.01*range(all_dists_subset);
            for idx = 1:npairs_valid
                np = idx;
                temp_data = filtfilt(b, a, all_egfs_subset{np});
                temp_data = temp_data.*w_taper_egf;
                plot(taxis_egf,all_dists_subset(np)+wiggle_scale*temp_data./max(abs(temp_data)),'k')
            end
            %plot the best velocity line
            plot(taxis_egf,best_v_source(k)*taxis_egf,'r')
            ylim([min(all_dists_subset)-wiggle_scale max(all_dists_subset)+wiggle_scale])
            box on;
            hold off
            title([source,' at ',num2str(fband(k,1)),'-',num2str(fband(k,2)),' Hz: ',num2str(best_v_source(k)),' km/s'])
            xlabel('time (s)')
            ylabel('distance (km)')
            set(gca,'TickDir','out')

            if save_fig
                saveas(gcf, fullfile(plotdir, [source, '_moveout_velocity.png']))
            end
        end
    end

    if save_delays
        outfilename = [outdir '/' char(source) '.dat'];
        if exist(outfilename,'file'); delete(outfilename); end
        fidout = fopen(outfilename,'a+');
    end

    for n = 1:npairs
        synfile = synfilelist(n).name;
        stemp = strsplit(synfile(1:end - length(synfileext)), '.to.');
        src = stemp{1}; rcv = stemp{2};
        if strcmp(src, rcv); continue; end
        stnpair = [src, '_', rcv];
        
        f_N_stack = dir([egfdir,'/',source,'/',stnpair,'*',egf_comp,'_N_stack.h5']);
        f_P_stack = dir([egfdir,'/',source,'/',stnpair,'*',egf_comp,'_P_stack.h5']);
        if isempty(f_N_stack) || isempty(f_P_stack); continue; end
        
        disp(['working on pair: ' stnpair ' (' num2str(n) '/' num2str(npairs) ')'])
        
        d_N = extract_corrdata_asdf([egfdir,'/',source,'/',f_N_stack.name]);
        d_P = extract_corrdata_asdf([egfdir,'/',source,'/',f_P_stack.name]);
        d_N_sub = extract_corrdata_asdf([egfdir,'/',source,'/',strrep(f_N_stack.name,'_stack','')]);
        d_P_sub = extract_corrdata_asdf([egfdir,'/',source,'/',strrep(f_P_stack.name,'_stack','')]);
        
        dist = d_P.value{1}.(egf_comp).dist;
        nsubstack = length(d_P_sub.value{1}.(egf_comp).time);
        if nsubstack < min_substack; continue; end
        
        synsacdata = readsac([syndir,'/',source,'/',synfile]);
        syn = synsacdata.DATA1 - mean(synsacdata.DATA1);
        syn = syn .* tukeywin(length(syn), syntaperfraction);

        % Pre-filter buffers
        egf_pos = d_P.value{1}.(egf_comp).data .* w_taper_egf;
        egf_neg = flipud(d_N.value{1}.(egf_comp).data) .* w_taper_egf;
        egf_sub_pos = d_P_sub.value{1}.(egf_comp).data .* w_taper_egf;
        egf_sub_neg = flipud(d_N_sub.value{1}.(egf_comp).data) .* w_taper_egf;

        % Define general shared window
        tmin0 = dist / max(best_v_source);
        tmin = max(0.9 * d_taper_time, tmin0 - 3 * tfmin(1));
        tmax = min(tmaximum, dist / min(best_v_source) + 3.5 * max(tfmin));
        
        ttuni = 0:dt_resample:tmaximum;
        ntuni = length(ttuni);
        egffb_weighted = zeros(ntuni, nfb);
        egffb_sub_weighted = zeros(ntuni, nsubstack, nfb);
        synfb = zeros(ntuni, nfb);
        
        snr_pos = zeros(nfb, 1); snr_neg = zeros(nfb, 1);
        t1_band = zeros(nfb, 1); t2_band = zeros(nfb, 1);

        for k = 1:nfb
            [b, a] = butter(N, [fband(k,1) fband(k,2)] / (1/dt_egf/2));
            f_pos = filtfilt(b, a, egf_pos);
            f_neg = filtfilt(b, a, egf_neg);
            
            % Refined window for SNR and measurement
            t_center = dist / best_v_source(k);
            t1_band(k) = t_center - 3.0 * tfmin(k);
            t2_band(k) = t_center + 3.0 * tfmin(k);
            it_sig = round(t1_band(k)/dt_egf):round(t2_band(k)/dt_egf);
            it_sig = it_sig(it_sig > 0 & it_sig <= length(f_pos));
            
            % Adaptive Noise Window
            min_n = round(3 * tfmin(k) / dt_egf);
            if (length(f_pos) - max(it_sig)) >= min_n
                n_idx = (length(f_pos)-min_n+1):length(f_pos);
            elseif min(it_sig) > (min_n + round(0.1/dt_egf))
                n_idx = round(0.1/dt_egf):(round(0.1/dt_egf)+min_n);
            else
                n_idx = round(0.9*length(f_pos)):length(f_pos);
            end
            
            snr_p = max(abs(f_pos(it_sig))) / rms(f_pos(n_idx));
            snr_n = max(abs(f_neg(it_sig))) / rms(f_neg(n_idx));
            snr_pos(k) = snr_p; snr_neg(k) = snr_n;
            
            % Weighted Combination (SNR Squared)
            w_p = snr_p^2 / (snr_p^2 + snr_n^2);
            w_n = snr_n^2 / (snr_p^2 + snr_n^2);
            combined_f = w_p * f_pos + w_n * f_neg;
            
            % get SNR for combined data.
            snr(k) = max(abs(combined_f(it_sig))) / rms(combined_f(n_idx));

            % Apply to stacks and resample
            egffb_weighted(:,k) = interp1(taxis_egf, combined_f, ttuni, 'linear', 0);
            
            [b_s, a_s] = butter(N, [fband(k,1) fband(k,2)] / (1/synsacdata.DELTA/2));
            f_syn = filtfilt(b_s, a_s, syn);
            synfb(:,k) = interp1(synsacdata.B + (0:synsacdata.NPTS-1)*synsacdata.DELTA, f_syn, ttuni, 'linear', 0);
            
            % Apply same weights to substacks
            for im = 1:nsubstack
                f_sub_p = filtfilt(b, a, egf_sub_pos(:,im));
                f_sub_n = filtfilt(b, a, egf_sub_neg(:,im));
                combined_sub = w_p * f_sub_p + w_n * f_sub_n;
                egffb_sub_weighted(:,im,k) = interp1(taxis_egf, combined_sub, ttuni, 'linear', 0);
            end
        end

        % Prepare segments for measurement
        itmin = round(tmin/dt_resample); itmax = round(tmax/dt_resample);
        egfsig = egffb_weighted(itmin:itmax, :);
        synsig = synfb(itmin:itmax, :);
        egfsubstacksig = egffb_sub_weighted(itmin:itmax, :, :);
        ntsig = size(egfsig, 1);
        
        maxdelay = min(max_dV * tmin0, max_dT);
        maxlag = round(maxdelay/dt_resample);
        ttc = (-maxlag:maxlag)*dt_resample;
        phase = nan(nfb,1); phaseerr = nan(nfb,1); rxc = nan(nfb,1); snr = nan(nfb,1);
        ec=zeros(nfb,1);
        xc=[];
        for k = 1:nfb
            t1_rel = t1_band(k) - tmin;
            t2_rel = t2_band(k) - tmin;
            it1 = max(1, round(t1_rel/dt_resample));
            it2 = min(ntsig, round(t2_rel/dt_resample));
            
            w_win = zeros(ntsig, 1);
            w_win(it1:it2) = tukeywin(it2-it1+1, taperfraction);
            
            d_win = egfsig(:,k) .* w_win;
            s_win = synsig(:,k) .* w_win;
            
            [xc(:,k), lags] = xcorr(d_win, s_win, maxlag);
            [xc_val, i_max] = max(xc(:,k));
            
            % Correction for ellipticity
            ec(k) = (dist - geo2dist_ellipse(d_N.value{1}.(egf_comp).lat(1),d_N.value{1}.(egf_comp).lon(1),...
                d_N.value{1}.(egf_comp).lat(2),d_N.value{1}.(egf_comp).lon(2))) / best_v_source(k);
            phase(k) = lags(i_max) * dt_resample + ec(k);
            
            % Substacks for error
            phasem = zeros(nsubstack, 1);
            for im = 1:nsubstack
                sub_win = egfsubstacksig(:,im,k) .* w_win;
                [xc_sub, lags_sub] = xcorr(sub_win, s_win, maxlag);
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
            rxc(k) = xc_val / (norm(d_win)*norm(s_win));
        end

        if fig_flag
            figure('Position',[100 0 1100 800]);
            for k = 1:nfb
                subplot(nfb, 2, 2*k-1); hold on;
                plot(ttuni, synfb(:,k)/max(abs(synfb(:,k))), 'r', 'LineWidth', 1.2);
                plot(ttuni, egffb_weighted(:,k)/max(abs(egffb_weighted(:,k))), 'k');
                line([t1_band(k) t1_band(k)], [-1 1], 'Color', 'm', 'LineStyle', '--');
                line([t2_band(k) t2_band(k)], [-1 1], 'Color', 'm', 'LineStyle', '--');
                xlim([tmin tmax]); ylabel([num2str(fband(k,1)) '-' num2str(fband(k,2)) 'Hz']);
                if k==1; legend('Syn','SNR-Weighted EGF'); end
                box on;
                axis on;
                
                subplot(nfb, 2, 2*k);
                plot(ttc+ec(k),xc(:,k)/max(abs(xc(:,k))),'k-');hold on % for ellipticity correction
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
            sgtitle([strrep(stnpair,'_','-') ' | One-sided | ',num2str(1000*dist),' m']);
            if save_fig
                saveas(gcf, fullfile(plotdir, [stnpair '.png']))
            else
                pause
            end
            close all;
        end

        if save_delays
            for k = 1:nfb
                if rxc(k) >= xcoeff_cutoff && snr(k) >= snr_cutoff && abs(phase(k)) <= min_wavelength*maxdelay
                    xcid=[char(src) '/bp' num2str(fband(k,1)) '_' num2str(fband(k,2)) '/' stnpair '_BHZ.P2.CORR.T1T2.SAC'];
                    fprintf(fidout,'%s %6.3f %2.0f %s %7.2f %7.2f %s %5.2f %5.2f %5.1f\n',...
                        xcid, phase(k), 1, 'RL', tmin+t1_band(k), tmin+t2_band(k), ['f' num2str(k)], phaseerr(k), rxc(k), snr(k)); 
                end
            end
        end
    end
    if save_delays; fclose(fidout); end
end