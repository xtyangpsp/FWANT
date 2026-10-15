% rotate_fx_fy_TT_RR: combine the fx and fy FD3Dspher runs for one virtual
% source and build the transverse-transverse (TT, Love) and radial-radial
% (RR, Rayleigh) responses at every receiver on a line.
%
% RUN THIS FROM THE PARENT FOLDER that contains the fx/ and fy/ run folders:
%
%     parent/
%       fx/  (SeisFD3D.conf, SeisSource.conf, input/, output/)
%       fy/  (same layout)
%       rotate_fx_fy_TT_RR.m   <-- this script
%
% Frame of the code: x = colatitude direction (points SOUTH), y = EAST,
% z = radial (up).  Run "fx" has a unit force toward +x, run "fy" toward +y.
%
% Notation: G_AB = component A at the receiver due to a force in direction B
% at the source (A,B in T,R,Z).  With alpha = source->receiver azimuth and
% beta = receiver->source (back) azimuth, both clockwise from north, the unit
% vectors in the (x,y) = (south,east) frame are
%
%   at source:    R_s = (-cos a,  sin a)   propagation direction (toward receiver)
%                 T_s = ( sin a,  cos a)   clockwise from R_s
%   at receiver:  R_r = ( cos b, -sin b)   propagation direction (away from source)
%                 T_r = (-sin b, -cos b)   clockwise from R_r
%
% A force along the unit vector f gives  u(f) = f(1)*u_fx + f(2)*u_fy  (linearity),
% so  G_TT = T_r . u(T_s),   G_RR = R_r . u(R_s),  and the "leakage" terms
% G_RT = R_r . u(T_s),  G_TR = T_r . u(R_s),  G_ZT = Z of u(T_s).
%
% Both transverse directions use the same handedness relative to the
% propagation direction, so the sign convention matches rotating E/N
% cross-correlations with obspy's rotate_ne_rt at both ends.
%
% Requires the FD3Dspher MATLAB helpers used by draw_seismo.m
% (locate_seismo, retrieve_seismo, retrieve_station, get_fnm_station, nc_varget)
% on the path.  Needs MATLAB R2016b or newer (implicit expansion, local functions
% in scripts) and the Signal Processing Toolbox (butter, filtfilt).

clear; close all

%% ===================== user parameters =====================
run_x    = 'fx';            % folder of the run with force toward +x (south)
run_y    = 'fy';            % folder of the run with force toward +y (east)
line_id  = 2;               % receiver line id (ryid in draw_seismo.m)
sta_ids  = 2:80;            % receiver indices on that line (stid in draw_seismo.m)

mfiles_root = '';           % folder holding fun-spool/ and saclab/; '' = use set_mfiles_path

Fc            = [0.05 0.1]; % band-pass corners in Hz; [] = no filter
N             = 4;          % Butterworth order (butter(N,[lo hi]) -> 2N-pole band-pass)
use_filtfilt  = false;      % false = causal filter like draw_seismo.m, true = zero phase
out_type      = 'disp';     % 'disp' (cumtrapz of Vx/Vy/Vz as in draw_seismo.m) or 'vel'
norm_mode     = 'trace';    % wiggle scaling: 'trace' = each trace to 1, 'global' = one common max
wig_scale     = 0.5;
coord_units   = 'auto';     % units of the 'grid' variable for x,y: 'auto', 'rad' or 'deg'
R_earth_km    = 6371;
save_mat      = true;       % writes TT_RR_<run_x>_<run_y>.mat

%% ===================== paths =====================
if ~isempty(mfiles_root)
    addpath(genpath(mfiles_root));
else
    try, set_mfiles_path; catch, end
end
if exist('locate_seismo','file') ~= 2
    error('locate_seismo not found: set mfiles_root to the folder containing fun-spool/.');
end

runs = {run_x, run_y};
for r = 1:2
    if ~exist(fullfile(runs{r},'SeisFD3D.conf'),'file')
        error('Cannot find %s. Run this script from the parent of %s and %s.', ...
              fullfile(runs{r},'SeisFD3D.conf'), run_x, run_y);
    end
end

%% ===================== source: location, force direction, amplitude =====================
[csx, lsx, f0x, fvx] = read_anchor_force(fullfile(run_x,'SeisSource.conf'));
[csy, lsy, f0y, fvy] = read_anchor_force(fullfile(run_y,'SeisSource.conf'));

if abs(csx-csy) > 1e-6 || abs(lsx-lsy) > 1e-6
    error('fx and fy runs have different source locations (%.4f,%.4f) vs (%.4f,%.4f).', csx,lsx,csy,lsy);
end
if any(fvx(2:3) ~= 0) || fvx(1) == 0
    error('%s/SeisSource.conf should have force direction (nonzero, 0, 0); found [%g %g %g].', run_x, fvx);
end
if any(fvy([1 3]) ~= 0) || fvy(2) == 0
    error('%s/SeisSource.conf should have force direction (0, nonzero, 0); found [%g %g %g].', run_y, fvy);
end
scale_x = f0x*fvx(1);       % responses are divided by this so both runs are "per unit force"
scale_y = f0y*fvy(2);
if scale_x ~= scale_y
    fprintf('Note: force amplitudes differ (%g vs %g); each run is normalised to unit force.\n', scale_x, scale_y);
end
colat_s = csx;  lon_s = lsx;  lat_s = 90 - colat_s;
fprintf('Source: colat %.4f  lon %.4f  (lat %.4f)\n', colat_s, lon_s, lat_s);

%% ===================== load seismograms from both runs =====================
nrecv = numel(sta_ids);
D = cell(1,2);
for r = 1:2
    root     = runs{r};
    fnm_conf = fullfile(root,'SeisFD3D.conf');
    dir_sta  = fullfile(root,'input');
    dir_out  = fullfile(root,'output');

    for n = 1:nrecv
        si = locate_seismo(fnm_conf, line_id, sta_ids(n), dir_sta, @get_fnm_station);
        g  = retrieve_station(si,'grid','stationdir',dir_sta);
        Sx = retrieve_seismo(si,'Vx','outdir',dir_out);
        Sy = retrieve_seismo(si,'Vy','outdir',dir_out);
        [Sz,T] = retrieve_seismo(si,'Vz','outdir',dir_out);
        if n == 1
            nt = numel(T);
            D{r}.Vx = zeros(nrecv,nt);  D{r}.Vy = zeros(nrecv,nt);  D{r}.Vz = zeros(nrecv,nt);
            D{r}.grid = zeros(nrecv,3); D{r}.T = T(:)';
        end
        D{r}.Vx(n,:) = Sx(:)';  D{r}.Vy(n,:) = Sy(:)';  D{r}.Vz(n,:) = Sz(:)';
        D{r}.grid(n,:) = g(:)';
    end
    fprintf('Loaded %d receivers from %s\n', nrecv, root);
end

if numel(D{1}.T) ~= numel(D{2}.T) || max(abs(D{1}.T - D{2}.T)) > 1e-6
    error('fx and fy runs have different time axes.');
end
if max(max(abs(D{1}.grid - D{2}.grid))) > 1e-9
    error('Receiver positions differ between the fx and fy runs.');
end
T = D{1}.T;

% per-unit-force responses
Fx.Vx = D{1}.Vx/scale_x;  Fx.Vy = D{1}.Vy/scale_x;  Fx.Vz = D{1}.Vz/scale_x;
Fy.Vx = D{2}.Vx/scale_y;  Fy.Vy = D{2}.Vy/scale_y;  Fy.Vz = D{2}.Vz/scale_y;

%% ===================== geometry: azimuths and distances =====================
colat_r = D{1}.grid(:,1);  lon_r = D{1}.grid(:,2);
switch coord_units
    case 'auto', is_rad = max(colat_r) <= pi;     % colatitude in radians is <= pi
    case 'rad',  is_rad = true;
    otherwise,   is_rad = false;
end
if is_rad, colat_r = colat_r*180/pi;  lon_r = lon_r*180/pi; end
lat_r = 90 - colat_r;
fprintf('Receivers: colat %.3f..%.3f, lon %.3f..%.3f (units interpreted as %s)\n', ...
        min(colat_r), max(colat_r), min(lon_r), max(lon_r), ternary(is_rad,'radians','degrees'));

alpha = sph_azimuth(lat_s, lon_s, lat_r, lon_r);   % source -> receiver, deg clockwise from north
beta  = sph_azimuth(lat_r, lon_r, lat_s, lon_s);   % receiver -> source (back azimuth)
gcd   = acosd(min(1, max(-1, sind(lat_s)*sind(lat_r) + cosd(lat_s)*cosd(lat_r).*cosd(lon_r-lon_s))));
dist_km = gcd*pi/180*R_earth_km;
alpha(dist_km < 1) = NaN;  beta(dist_km < 1) = NaN;     % azimuth undefined at the source
fprintf('Epicentral distance %.1f .. %.1f km\n', min(dist_km), max(dist_km));

Ts = [ sind(alpha)  cosd(alpha)];     % transverse force direction at source, (x,y)
Rs = [-cosd(alpha)  sind(alpha)];     % radial force direction at source
Tr = [-sind(beta)  -cosd(beta)];      % transverse direction at receiver
Rr = [ cosd(beta)  -sind(beta)];      % radial direction at receiver

%% ===================== combine the two runs and rotate =====================
% response to a force along T_s (transverse) and along R_s (radial)
uxT = Ts(:,1).*Fx.Vx + Ts(:,2).*Fy.Vx;
uyT = Ts(:,1).*Fx.Vy + Ts(:,2).*Fy.Vy;
uzT = Ts(:,1).*Fx.Vz + Ts(:,2).*Fy.Vz;
uxR = Rs(:,1).*Fx.Vx + Rs(:,2).*Fy.Vx;
uyR = Rs(:,1).*Fx.Vy + Rs(:,2).*Fy.Vy;
uzR = Rs(:,1).*Fx.Vz + Rs(:,2).*Fy.Vz;

TT = Tr(:,1).*uxT + Tr(:,2).*uyT;     % Love
RR = Rr(:,1).*uxR + Rr(:,2).*uyR;     % Rayleigh, radial
RT = Rr(:,1).*uxT + Rr(:,2).*uyT;     % leakage: R from transverse force
TR = Tr(:,1).*uxR + Tr(:,2).*uyR;     % leakage: T from radial force
ZT = uzT;                             % leakage: Z from transverse force
ZR = uzR;                             % Rayleigh, vertical

% integrate to displacement and filter, the same way draw_seismo.m does
pp = @(W) postproc(W, T, out_type, Fc, N, use_filtfilt);
TT = pp(TT);  RR = pp(RR);  RT = pp(RT);  TR = pp(TR);  ZT = pp(ZT);  ZR = pp(ZR);

%% ===================== diagnostics =====================
pk = @(W) max(abs(W),[],2);
leak_RT = pk(RT)./pk(TT);             % should be << 1 if frames/signs are right
leak_ZT = pk(ZT)./pk(TT);
leak_TR = pk(TR)./pk(RR);

% azimuth estimated from the vertical components alone: Z(fx) ~ -cos(a)*A, Z(fy) ~ sin(a)*A
Zx = pp(Fx.Vz);  Zy = pp(Fy.Vz);
a_est = zeros(nrecv,1);
for k = 1:nrecv
    M = [Zx(k,:)*Zx(k,:)', Zx(k,:)*Zy(k,:)'; Zx(k,:)*Zy(k,:)', Zy(k,:)*Zy(k,:)'];
    [V,Dg] = eig(M);  [~,im] = max(diag(Dg));  v = V(:,im);
    a_est(k) = mod(atan2d(v(2), -v(1)), 180);
end
a_err = mod(a_est - mod(alpha,180) + 90, 180) - 90;      % deg, mod-180 difference

fprintf('\n--- diagnostics (medians over receivers) ---\n');
fprintf('peak|RT|/peak|TT| = %.3f   peak|ZT|/peak|TT| = %.3f   peak|TR|/peak|RR| = %.3f\n', ...
        median(leak_RT,'omitnan'), median(leak_ZT,'omitnan'), median(leak_TR,'omitnan'));
fprintf('azimuth from Z vs geodesic: median |error| = %.1f deg, max = %.1f deg\n', ...
        median(abs(a_err),'omitnan'), max(abs(a_err)));
fprintf('Leakage ratios near 1 or an azimuth error near 0/90 deg with sign flips means the x,y frame or a sign is off.\n');

%% ===================== figure 1: TT and RR =====================
idx = 1:nrecv;  ytk = 1:max(1,round(nrecv/10)):nrecv;
ytl = arrayfun(@num2str, sta_ids(ytk), 'UniformOutput', false);
figure('Position',[100 100 1300 750]);
ax1 = subplot(1,2,1);  wiggle(T, TT, wig_scale, normvec(TT,norm_mode));
title('TT (Love)');  xlabel('Time (s)');  ylabel('Receiver');  set(ax1,'YTick',ytk,'YTickLabel',ytl);  grid on
ax2 = subplot(1,2,2);  wiggle(T, RR, wig_scale, normvec(RR,norm_mode));
title('RR (Rayleigh, radial)');  xlabel('Time (s)');  set(ax2,'YTick',ytk,'YTickLabel',ytl);  grid on
linkaxes([ax1 ax2],'xy');

%% ===================== figure 2: geometry and checks =====================
figure('Position',[100 100 1300 420]);
subplot(1,3,1);
plot(idx,alpha,'-o',idx,beta,'-s');  grid on;  legend('\alpha (src\rightarrowrec)','\beta (rec\rightarrowsrc)','Location','best');
xlabel('Receiver (position in list)');  ylabel('Azimuth (deg from N)');  title('Path azimuths');
subplot(1,3,2);
semilogy(idx,leak_RT,'-o',idx,leak_ZT,'-s',idx,leak_TR,'-^');  grid on;
legend('RT/TT','ZT/TT','TR/RR','Location','best');
xlabel('Receiver (position in list)');  ylabel('Peak amplitude ratio');  title('Leakage after rotation');
subplot(1,3,3);
plot(idx,mod(alpha,180),'k-',idx,a_est,'r.','MarkerSize',12);  grid on;
legend('geodesic \alpha mod 180','from Z(fx), Z(fy)','Location','best');
xlabel('Receiver (position in list)');  ylabel('Azimuth mod 180 (deg)');
title(sprintf('Azimuth check (median |err| %.1f deg)', median(abs(a_err),'omitnan')));

%% ===================== figure 3: leakage waveforms, scaled to the main phase =====================
figure('Position',[100 100 1300 750]);
ax3 = subplot(1,3,1);  wiggle(T, RT, wig_scale, pk(TT));  title('RT (R from transverse force), scaled by max|TT|');
xlabel('Time (s)');  ylabel('Receiver');  set(ax3,'YTick',ytk,'YTickLabel',ytl);  grid on
ax4 = subplot(1,3,2);  wiggle(T, ZT, wig_scale, pk(TT));  title('ZT, scaled by max|TT|');
xlabel('Time (s)');  set(ax4,'YTick',ytk,'YTickLabel',ytl);  grid on
ax5 = subplot(1,3,3);  wiggle(T, TR, wig_scale, pk(RR));  title('TR (T from radial force), scaled by max|RR|');
xlabel('Time (s)');  set(ax5,'YTick',ytk,'YTickLabel',ytl);  grid on
linkaxes([ax3 ax4 ax5],'xy');

%% ===================== save =====================
if save_mat
    fnm_out = sprintf('TT_RR_%s_%s.mat', run_x, run_y);
    save(fnm_out,'T','TT','RR','RT','TR','ZT','ZR','alpha','beta','dist_km', ...
         'colat_r','lon_r','colat_s','lon_s','sta_ids','line_id','out_type','Fc','N');
    fprintf('Saved %s\n', fnm_out);
end

%% ===================== local functions =====================
function [colat,lon,f0,fvec] = read_anchor_force(fnm)
% First force line after <anchor_force>: x y z start f0 fx fy fz  (x = colatitude, y = longitude, degrees)
fid = fopen(fnm,'r');
if fid < 0, error('Cannot open %s', fnm); end
found = false;  colat = NaN;  lon = NaN;  f0 = NaN;  fvec = [NaN NaN NaN];
while true
    s = fgetl(fid);
    if ~ischar(s), break; end
    s = strtrim(s);
    if ~found
        found = strncmp(s,'<anchor_force>',14);
        continue
    end
    ic = find(s == '#', 1);  if ~isempty(ic), s = strtrim(s(1:ic-1)); end
    if isempty(s), continue; end
    v = sscanf(s,'%f');
    if numel(v) < 8, error('Unexpected force line in %s: %s', fnm, s); end
    colat = v(1);  lon = v(2);  f0 = v(5);  fvec = v(6:8)';
    break
end
fclose(fid);
if isnan(colat), error('No force line found after <anchor_force> in %s', fnm); end
end

function az = sph_azimuth(lat1, lon1, lat2, lon2)
% Azimuth (deg clockwise from north) at point 1 toward point 2 on a sphere.
dl = lon2 - lon1;
y  = sind(dl).*cosd(lat2);
x  = cosd(lat1).*sind(lat2) - sind(lat1).*cosd(lat2).*cosd(dl);
az = mod(atan2d(y,x), 360);
end

function W = postproc(W, T, out_type, Fc, N, use_filtfilt)
% Rows are receivers. Optional integration to displacement, then optional band-pass.
dt = T(2) - T(1);
if strcmp(out_type,'disp'), W = cumtrapz(W,2)*dt; end
if ~isempty(Fc)
    [b,a] = butter(N, Fc/(0.5/dt));
    if use_filtfilt, W = filtfilt(b,a,W')'; else, W = filter(b,a,W,[],2); end
end
end

function nv = normvec(W, mode)
m = max(abs(W),[],2);
if strcmp(mode,'global'), m = max(m)*ones(size(m)); end
nv = m;
end

function wiggle(T, W, scale, nv)
nv(nv == 0 | isnan(nv)) = 1;
hold on
for k = 1:size(W,1)
    plot(T, scale*W(k,:)/nv(k) + k, 'b');
end
hold off
ylim([0 size(W,1)+1]);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
