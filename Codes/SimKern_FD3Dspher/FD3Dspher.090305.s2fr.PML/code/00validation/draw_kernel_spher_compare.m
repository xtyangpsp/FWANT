% draw_kernel_spher_compare
%
% Plot kernels from two runs side by side (old code vs new code), one row per
% kernel variable, optionally with a third column showing the difference.
% Run this script from the folder that contains sim.kernel and
% sim.kernel_validate (and sim.input, sim.station).

close all

% =====================================================================
%  USER-SPECIFIED PARAMETERS
% =====================================================================

% ---- kernel runs to compare (left -> right) ----
% Each run is a LIST of folders plus a list of weights; the plotted kernel is
% sum_j ker_w{k}(j) * kernel(ker_dir{k}{j}).
%  - ZZ check : one folder per run, weight 1 (default below).
%  - TT check : the old code cannot mix the fx and fy source runs, but it is linear in the
%    source run: run the old code twice (event side = fx run -> sim.kernel_fx, event side =
%    fy run -> sim.kernel_fy; same conf line, same mixed fT synthetics, source weights NOT
%    given) and combine with the weights [sin(alpha) cos(alpha)] of the pair; the new code
%    (conf line WITH the two weights) is run once in sim.kernel_validate.
ker_dir   = {{'sim.kernel'}, {'sim.kernel_validate'}};   % parent folder of <source>/<receiver>/...
ker_w     = {[1],            [1]};
ker_label = {'old code',     'new code'};                % panel labels
% TT example:
% ker_dir = {{'sim.kernel_fx','sim.kernel_fy'}, {'sim.kernel_validate'}};
% ker_w   = {[sind(40) cosd(40)],               [1]};   % sin(alpha), cos(alpha) of the pair
% ker_label = {'old code (fx,fy combined)',     'new code'};

% ---- pair, component, frequency band ----
source    = 'TA.P26A';
receiver  = 'TA.P23A';
comp      = 'BHZ';      % kernel component folder (BHZ for ZZ; BHT for TT if you use that label)
freq_id   = 1;
kern_tag  = 'T1T2.P2';  % kernel subfolder name
pband     = [7.5 15; 10 20; 15 30; 20 40; 30 55; 40 70];   % period bands (s)

% ---- variables to plot and color limits (one row each) ----
var_list = {}; scl_caxis_list = {};
var_list{end+1}='phase_Vp'; scl_caxis_list{end+1}=[-1.5e-14,1.5e-14];
var_list{end+1}='phase_Vs'; scl_caxis_list{end+1}=[-1.5e-13,1.5e-13];
%var_list{end+1}='Kbpz'; scl_caxis_list{end+1}=[-1e-11,1e-11];
%var_list{end+1}='Kbqz'; scl_caxis_list{end+1}=[-2e-11,2e-11];

% ---- difference column (new - old), color limits = diff_frac * kernel limits ----
flag_diff = 1;          % 1: add a third column with (run2 - run1)
diff_frac = 0.1;

% ---- model / snapshot ----
RUN_ROOT = 'sim.station/skel/fx/';
dir_coord = './sim.input';
snap_id = 1;                         % snapshot id (see SeisFD3D.conf)
subs=[1 1 42]; subc=[-1 -1 1]; subt=[1 1 1];   % [x,y,z] start / count / stride (z=41: horizontal slice)

% ---- plot options ----
flag_jetwr  = 1;        % use 'jetwr' colormap
flag_km     = 1;        % axes in km
flag_print  = 1;        % save png
MFILEROOT = '../mfiles';
view_angle = [0 0];     % view(az,el)

% =====================================================================
%  DERIVED PARAMETERS (no need to edit below)
% =====================================================================

path([MFILEROOT '/fun-spool'],path);

fnm_conf   = [RUN_ROOT '/' 'SeisFD3D.conf'];
dir_metric = dir_coord;

nrun = length(ker_dir);
nvar = length(var_list);
ncol = nrun + (flag_diff==1 && nrun==2);

fband = flip(1./pband,2);
freq_key = ['f' num2str(freq_id)];
period_tag = [num2str(1/fband(freq_id,2)),' - ',num2str(1/fband(freq_id,1)),' s'];

% kernel folders of each run
pnm_out = cell(1,nrun);
for k=1:nrun
    if ischar(ker_dir{k}), ker_dir{k} = {ker_dir{k}}; end
    if length(ker_w{k})~=length(ker_dir{k}), error('ker_w{%d} must match ker_dir{%d}',k,k); end
    for j=1:length(ker_dir{k})
        pnm_out{k}{j} = [ker_dir{k}{j},'/',source,'/',receiver,'/',comp,'/',num2str(freq_id),'/',kern_tag,'/'];
        if ~exist(pnm_out{k}{j},'dir')
            error('kernel folder not found: %s', pnm_out{k}{j});
        end
    end
end

% -------------------- coordinates --------------------
[snapinfo] = locate_snap(fnm_conf,snap_id,'start',subs,'count',subc,'stride',subt);
[X,Y,Z] = gather_coord(snapinfo,'coorddir',dir_metric);
R=Z; CLAT=X; LON=Y;

[y,x,z] = sph2cart(LON,pi/2-CLAT,R);
str_unit = 'm';
if flag_km
    x=x/1e3; y=y/1e3; z=z/1e3; str_unit='km';
end
xp = squeeze(permute(y,[2 1 3]));
yp = squeeze(permute(x,[2 1 3]));
zp = squeeze(permute(z,[2 1 3]));

% depth of horizontal surface in km
dep = (6371000-R(1,1))/1.e3;

% -------------------- plot kernels --------------------
figure('position',[100 100 450*ncol 380*nvar]);
for n=1:nvar
    % load the kernel of each run
    V = cell(1,nrun);
    for k=1:nrun
        V{k} = 0;
        for j=1:length(pnm_out{k})
            [Vk,~] = gather_dist(snapinfo,snap_id,var_list{n},'outdir',pnm_out{k}{j});
            V{k} = V{k} + ker_w{k}(j)*double(Vk);
        end
    end

    if strcmp(var_list{n},'phase_Vs')
        phase = 'Vs';
    elseif strcmp(var_list{n},'phase_Vp')
        phase = 'Vp';
    else
        phase = strrep(var_list{n},'_','\_');
    end

    % print how close the two runs are
    if nrun==2
        d = V{2}-V{1};
        fprintf('%s: max|old|=%.3e  max|new|=%.3e  max|new-old|=%.3e  (rel. to max|old|: %.2e)\n', ...
            var_list{n}, max(abs(V{1}(:))), max(abs(V{2}(:))), max(abs(d(:))), ...
            max(abs(d(:)))/max(abs(V{1}(:))));
    end

    for c=1:ncol
        subplot(nvar,ncol,(n-1)*ncol+c)
        if c<=nrun
            Vplot = V{c};
            clim = scl_caxis_list{n};
            ttl = [ker_label{c} ': ' phase ', ' num2str(round(dep)) ' km, ' period_tag];
        else
            Vplot = V{2}-V{1};
            clim = diff_frac*scl_caxis_list{n};
            ttl = ['new - old (x' num2str(1/diff_frac) '): ' phase];
        end

        surf(xp,yp,zp,squeeze(permute(Vplot,[2 1 3])));
        shading flat;
        if flag_jetwr==1
            colormap(gca,'jetwr');
        end
        caxis(clim);
        xlabel(['X (' str_unit ')']);
        view(view_angle(1),view_angle(2));
        set(gca,'box','off');
        camlight(0,0,'local');
        lighting phong
        material dull
        grid off;

        cb = colorbar;
        cb.Location = 'eastoutside';
        cpos = cb.Position;
        cpos(2) = 1.1*cpos(2);
        cpos(3:4) = 0.5*cpos(3:4);
        cb.Position = cpos;

        title(ttl,'FontSize',12);
    end
    drawnow
end

sgtitle([source ' -> ' receiver ', ' comp ', ' freq_key]);

if flag_print==1
    print(gcf,'-dpng',['kernel_compare_',source,'_',receiver,'_',freq_key,'.png']);
end
