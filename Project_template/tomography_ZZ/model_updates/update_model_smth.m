% matlab script to update the velocity model using inversion results
% last modified from, n.cascadia/update_model_smth.m, Y.Shen, 03-13-2011
% script to update the velocity model
% Restructured by Xiaotao Yang for better user experiences.
%
% Revision notes by Xiaotao Yang on 20261124:
%  - all smoothing, clipping, and file-handling settings are in the parameter block below
%  - model smoothing can be 'box' (as before) or 'gaussian', and an optional EXTRA vertical smoothing can be
%    applied only for selected iterations (for example once at ite1) so it cannot be forgotten later
%  - smoothing ignores the cells that are kept at the reference values ("fixed" cells) and handles the grid edges
%  - the safety checks use the final (smoothed) model; the stability trap no longer calls undefined variables
%  - settings and diagnostics are written to a log file next to the updated media
%  - existing updated media can be backed up automatically before they are overwritten
%  - the coordinate/block loops are vectorized (same ordering as the loops, so identical results)
%  - the maxdv clip is now an option (apply_maxdv_clip), OFF by default
%  - dv_scale scales the perturbations from the inversion before they are applied (1 = no scaling)
%  - the sigma of the Gaussian distance weighting is renamed weight_sigma, so it is not confused with the
%    smoothing parameters (smooth_sigma, extra_vsmooth_sigma)
% With the smoothing settings box [7 7 1], no extra smoothing, no fixed cells, apply_maxdv_clip = 1, and
% dv_scale = 1, the result is identical to the previous version of this script, EXCEPT for the min_vs_cut
% fix (see the parameter block).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up global parameters in this section.

close all; clear all
addpath(genpath('/depot/xtyang/data/codes/mexcdf'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise'))
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))
MFILEROOT='../mfiles';
path([MFILEROOT '/fun-spool'],path);

% ---- current model iteration and which inversion result to use ----
ite_nm = 'ite_0.025deg_01';
inv_homedir='inv.structure.kerVpVs';
inv_workdir='FreqAll.Z.model';
subsample=0; %0: single data set inversion. 1: subsample bootstrapping.
if subsample
	inv_modeldir='subsamples/model_statistics';
	modelfile_namebase='model_mean'; %'model_mean' for subsample inversion results. 'try' for normal results.
else
	inv_modeldir='result.1th';
	modelfile_namebase='try'; %'model_mean' for subsample inversion results. 'try' for normal results.
end

% regularization of the INVERSION (part of the result file name). Not the model smoothing below.
sel_smooth=8; % inversion smoothing
sel_damp=10;  % inversion damping (confirm: ZZ and TT runs may use different values)

% ---- limits and scaling of the model update ----
% minimum Vp and Vs allowed. model values below these will be reverted back to the input model.
min_vp_cut = 2000; % in meters/second. Applied to the REFERENCE Vp ("water layer"): those cells are never updated.
min_vs_cut = 1000; % in m/s. Applied to the UPDATED Vs: cells whose updated Vs falls below it are reverted.
% NOTE: earlier versions of this script compared the updated Vp (not Vs) with min_vs_cut in this rule, which
% almost never triggered. Now Vs is compared, as the comment intended. Cells with updated Vs < min_vs_cut
% (for example slow sediments) are reverted to the reference model, i.e. not updated, so results differ from
% iterations run with the old rule. The number of cells reverted by each rule is printed and logged.
% Lower min_vs_cut if it blocks basins you want to update.

% Scale of the perturbations (dVp/Vp and dVs/Vs) read from the inversion result, applied before the update.
% 1 = no scaling. For example 1.2 over-relaxes the step, 0.8 under-relaxes it. The density perturbation
% follows the scaled dVp. The scaling is applied before the maxdv clip, so maxdv (if used) bounds the
% perturbation that is actually applied.
dv_scale = 1;

% Optional clip of the applied perturbation to +/- maxdv. OFF by default. Earlier iterations were run with the
% clip ON (maxdv = 0.1); set apply_maxdv_clip = 1 to reproduce them. With the clip off, the number of cells
% above maxdv is still printed and logged, and the stability trap below still stops the script if Vp, Vs, or
% density become zero, negative, or NaN.
apply_maxdv_clip = 0;
maxdv = 0.1; % maximum allowed perturbation when apply_maxdv_clip = 1. 0.3=30% + or -.

drhdvp = 0.297; % density perturbation = drhdvp * Vp perturbation (drho/rho = 0.297 dVp/Vp)

% Optional Gaussian distance weighting of the perturbations (boundary stability). Unrelated to the smoothing
% parameters below. Cells farther from the nearest inversion point than weight_sigma are down-weighted.
apply_gaussian_weight = 0; % if 1, apply the distance-based gaussian weight
weight_sigma = 10000;      % width of the weighting in meters (adjust based on grid spacing)

% ---- smoothing of the updated model (applied to mu, lambda, and rho, as in the original script) ----
apply_smoothing = 1;
smooth_type  = 'box';   % 'box' (original) or 'gaussian'
smooth_box   = [7 7 3]; % box: size in grid cells along [dim1 dim2 dim3] = [lateral lateral vertical]; each odd
smooth_sigma = [2 2 1]; % gaussian: standard deviation in grid cells along [dim1 dim2 dim3]; 0 = no smoothing along it.
                        % [2 2 0] has about the same variance as the box [7 7 1] (sigma = N/sqrt(12))
smooth_preserve_fixed = 1; % 1: cells kept at reference values (Vp or Vs below the cuts) are excluded from the
                           % smoothing and restored exactly. 0: smooth them with their neighbours (original behavior)

% Optional EXTRA vertical smoothing, applied after the smoothing above, ONLY when ite_nm is in the list below.
% For example, a one-time larger vertical window at the ite1 update, to soften the interfaces inherited from
% the starting model. Set to {} to disable. The window is in grid cells; choose it from the vertical
% resolution (checkerboard/spike tests). The printout below shows its physical size.
extra_vsmooth_ites  = {'ite_0.025deg_01'};
extra_vsmooth_type  = 'gaussian'; % 'gaussian' or 'box'
extra_vsmooth_sigma = 3;          % gaussian: vertical sigma in cells (SET THIS from the resolution tests)
extra_vsmooth_box   = 9;          % box: vertical size in cells (odd)

% ---- file handling ----
backup_existing_updated = 1; % 1: copy an existing updated_input folder to <folder>_backup_<time> before writing
                             % (needs disk space for one extra copy of the media files)
write_model = 1;             % 0: run all steps and checks but do not write media files

% read previous model -----------------------------------
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./input_' ite_nm];
dir_metric=['./input_' ite_nm];
dir_media=['./input_' ite_nm];
dir_media_updated=['./updated_input_' ite_nm];

%read-in key configuration parameters
id_stress = 1; %the snap id for the stress tensor.
confinfo=read_fdconf(fnm_conf,id_stress);

%blockfile dimension
mx=confinfo.snap_subc(1);
my=confinfo.snap_subc(2);
mz=confinfo.snap_subc(3);

% set up dimensions for output
MPII=confinfo.dims(1);
MPIJ=confinfo.dims(2);
MPIK=confinfo.dims(3); %

assert(all(mod(smooth_box,2) == 1), 'smooth_box sizes must be odd');
assert(mod(extra_vsmooth_box,2) == 1, 'extra_vsmooth_box must be odd');
assert(any(strcmp(smooth_type,{'box','gaussian'})), 'smooth_type must be box or gaussian');
assert(any(strcmp(extra_vsmooth_type,{'box','gaussian'})), 'extra_vsmooth_type must be box or gaussian');
assert(isscalar(dv_scale) && isnumeric(dv_scale) && dv_scale >= 0, 'dv_scale must be a non-negative scalar');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Read and reshape simulation model. This takes a few minutes ...')

id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];
id{end+1} = 0; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,-1];subt{end+1}=[1,1,1];
               indxem{end+1}=[];
               indxkp{end+1}=[];
n=1;
[snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
[XSIM{n},YSIM{n},ZSIM{n}]=gather_coord(snapinfo{n},'coorddir',dir_coord);
[NX,NY,NZ]=size(XSIM{n});
SIZ=NX*NY*NZ;
% grid points as an N x 3 list. MATLAB's column-major ordering (i fastest, then j, then k) is the same as
% the ordering of the original triple loop.
XYZ=[XSIM{n}(:), YSIM{n}(:), ZSIM{n}(:)];

mrh=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp=((mla+2*mmu)./mrh).^0.5;
mvs=(mmu./mrh).^0.5;

% read inversion results ---------------------------------------
disp('Read inversion results and interpolate at simulation grids ...')

blockfile=strcat('block.',num2str(mx),'x',num2str(my),'x',num2str(mz),'.1x1x1.1x1x1.nc');
fnm_blk=strcat('../',ite_nm,'/',inv_homedir,'/',blockfile);

XB=nc_varget(fnm_blk,'x');
YB=nc_varget(fnm_blk,'y');
ZB=nc_varget(fnm_blk,'z');
xb=reshape(XB,[mx,my,mz]);
yb=reshape(YB,[mx,my,mz]);
zb=reshape(ZB,[mx,my,mz]);
clear XB YB ZB

% read velocity perturbation model
inv_modeldirfull=strcat('../',ite_nm,'/',inv_homedir,'/',inv_workdir,'/',inv_modeldir);

fnm_try=strcat(inv_modeldirfull,'/',modelfile_namebase,'.damp',num2str(sel_damp),'.smot',...
    num2str(sel_smooth),'.st0.ev0.lo0.dat');

M=load(fnm_try);
MV=reshape(M,mx,my,mz,2);
rvp=squeeze(MV(:,:,:,1));
rvs=squeeze(MV(:,:,:,2));
clear M MV

% add an imaginary layer with wave speeds of the top layer above the top layer,
% so the simulation grids near the surface can be updated (no NaN in interpolation)
xb(:,:,mz+1)=xb(:,:,mz);
yb(:,:,mz+1)=yb(:,:,mz);
zb(:,:,mz+1)=2*zb(:,:,mz)-zb(:,:,mz-1);
rvp(:,:,mz+1)=rvp(:,:,mz);
rvs(:,:,mz+1)=rvs(:,:,mz);

% put the inversion block locations in a three column array (same ordering as the original loops)
xyzb=double([xb(:), yb(:), zb(:)]);
vp=double(rvp(:)); vs=double(rvs(:));

% Note: TriScatteredInterp is replaced by scatteredInterpolant for stability
disp('Constructing interpolation function...');

% Use 'none' for extrapolation to identify points outside the inversion block
FP = scatteredInterpolant(xyzb, vp, 'linear', 'none');
FS = scatteredInterpolant(xyzb, vs, 'linear', 'none');

itp_vp_raw = FP(XYZ);
itp_vs_raw = FS(XYZ);

% 3. Initialize and Map perturbations
itp_vp_3d = zeros(NX, NY, NZ);
itp_vs_3d = zeros(NX, NY, NZ);

% Identify valid (non-NaN) indices
valid_idx = ~isnan(itp_vp_raw) & ~isnan(itp_vs_raw);
itp_vp_3d(valid_idx) = itp_vp_raw(valid_idx);
itp_vs_3d(valid_idx) = itp_vs_raw(valid_idx);

% 2. Optional Gaussian distance weighting of the perturbations, for boundary stability
if apply_gaussian_weight
	disp('Calculating Gaussian weighting mask...');
% Find distance to nearest inversion point.
% If knnsearch is missing, use: dist = pdist2(xyzb, XYZ, 'euclidean', 'Smallest', 1);
	[~, dist] = knnsearch(xyzb, XYZ, 'K', 1);

	weight = exp(-(dist.^2) / (2 * weight_sigma^2));
	weight(dist > 3*weight_sigma) = 0; % Hard cut-off for far-field
	weight_3d = reshape(weight, NX, NY, NZ);

% This ensures a smooth transition and no "cliff" edges
	itp_vp = itp_vp_3d .* weight_3d;
	itp_vs = itp_vs_3d .* weight_3d;
else
	itp_vp = itp_vp_3d;
	itp_vs = itp_vs_3d;
end

% 4. Scale the perturbations (dv_scale = 1: no scaling), then optionally clip them to +/- maxdv
itp_vp = dv_scale * itp_vp;
itp_vs = dv_scale * itp_vs;

nexc_vp = nnz(abs(itp_vp) > maxdv); % cells above maxdv (clipped only if apply_maxdv_clip = 1)
nexc_vs = nnz(abs(itp_vs) > maxdv);
if apply_maxdv_clip
	itp_vp = max(min(itp_vp, maxdv), -1.0*maxdv);
	itp_vs = max(min(itp_vs, maxdv), -1.0*maxdv);
	fprintf('maxdv clip ON (+/-%.2f): cells clipped: Vp %d (%.3f%%), Vs %d (%.3f%%)\n', maxdv, ...
	    nexc_vp, 100*nexc_vp/SIZ, nexc_vs, 100*nexc_vs/SIZ);
else
	fprintf('maxdv clip OFF: cells with |dv| > %.2f (not clipped): Vp %d (%.3f%%), Vs %d (%.3f%%)\n', maxdv, ...
	    nexc_vp, 100*nexc_vp/SIZ, nexc_vs, 100*nexc_vs/SIZ);
end
fprintf('dv_scale = %g. Applied perturbation: max |dv| Vp %.4f, Vs %.4f; rms Vp %.4f, Vs %.4f\n', dv_scale, ...
    max(abs(itp_vp(:))), max(abs(itp_vs(:))), sqrt(mean(itp_vp(:).^2)), sqrt(mean(itp_vs(:).^2)));

% 5. Update the model
outvp = mvp .* (1 + itp_vp);
outvs = mvs .* (1 + itp_vs);
outrh = mrh .* (1 + drhdvp * itp_vp);

% --- CRITICAL STABILITY TRAP ---
% This will stop the script BEFORE it writes the bad files
% (with the maxdv clip off, this is the protection against non-physical updates)
bad = (outvp <= 0) | (outvs <= 0) | (outrh <= 0) | isnan(outvp) | isnan(outvs) | isnan(outrh);
if any(bad(:))
    [iv,jv,kv] = ind2sub(size(outvp), find(bad,1));
    fprintf('!!! STABILITY ERROR DETECTED !!!\n');
    fprintf('Zero, negative, or NaN at i=%d, j=%d, k=%d (%d cells in total)\n', iv, jv, kv, nnz(bad));
    error('Aborting: The updated model contains values that will cause dt=0.');
end

% finally identify the water layer in the reference model and fix the velocity.
% "fixed" marks every cell that is kept at the reference values.
fixed = false(size(mvp));

indx = find(mvp <= min_vp_cut);
outvp(indx)=mvp(indx);
outvs(indx)=mvs(indx);
outrh(indx)=mrh(indx);
fixed(indx)=true;
clear indx

nrev_vp = nnz(fixed);                           % cells reverted by the water-layer rule (reference Vp <= min_vp_cut)
% updated Vs below the cut (earlier versions compared the updated Vp here, apparently a typo)
indx = find(outvs < min_vs_cut & ~fixed);
nrev_vs = numel(indx);
nrev_vs_refbelow = nnz(mvs(indx) < min_vs_cut); % of these, cells whose REFERENCE Vs is already below the cut
outvp(indx)=mvp(indx);
outvs(indx)=mvs(indx);
outrh(indx)=mrh(indx);
fixed(indx)=true;
clear indx
fprintf('cells reverted: reference Vp <= %g: %d; updated Vs < %g: %d (of which reference Vs already below the cut: %d)\n', ...
    min_vp_cut, nrev_vp, min_vs_cut, nrev_vs, nrev_vs_refbelow);
fprintf('cells kept at the reference values: %d (%.3f%%)\n', nnz(fixed), 100*nnz(fixed)/SIZ);

% convert vp,vs to mu,la
outmu=outrh.*outvs.^2;
outla=outrh.*outvp.^2-outmu*2;

% find grids with la < 0
indxla=find(outla<0);
fprintf('cells with lambda < 0 (set to 0): %d\n', numel(indxla));
outla(indxla)=0;

% smoothing --------------------------------------------------
do_extra = any(strcmp(ite_nm, extra_vsmooth_ites));
if apply_smoothing
    disp('smoothening the model ...')
    if smooth_preserve_fixed, fixed_sm = fixed; else, fixed_sm = false(size(fixed)); end

    report_smoothing(ZSIM{n}, 'base smoothing', smooth_type, smooth_box, smooth_sigma);
    outrh = smooth_masked(outrh, fixed_sm, smooth_type, smooth_box, smooth_sigma);
    outmu = smooth_masked(outmu, fixed_sm, smooth_type, smooth_box, smooth_sigma);
    outla = smooth_masked(outla, fixed_sm, smooth_type, smooth_box, smooth_sigma);

    if do_extra
        extra_box   = [1 1 extra_vsmooth_box];
        extra_sigma = [0 0 extra_vsmooth_sigma];
        fprintf('EXTRA vertical smoothing applied for %s\n', ite_nm);
        report_smoothing(ZSIM{n}, 'extra vertical smoothing', extra_vsmooth_type, extra_box, extra_sigma);
        outrh = smooth_masked(outrh, fixed_sm, extra_vsmooth_type, extra_box, extra_sigma);
        outmu = smooth_masked(outmu, fixed_sm, extra_vsmooth_type, extra_box, extra_sigma);
        outla = smooth_masked(outla, fixed_sm, extra_vsmooth_type, extra_box, extra_sigma);
    end
else
    disp('no smoothing applied');
end

[outx, outy, outz]=size(outvp);
ib=floor(outx/MPII); jb=floor(outy/MPIJ); kb=floor(outz/MPIK);

% FINAL STABILITY CHECK, on the model that will be written (after smoothing)
vp_f = sqrt((outla + 2*outmu)./outrh);
vs_f = sqrt(outmu./outrh);
fprintf('Final Vp: %.0f - %.0f m/s, Vs: %.0f - %.0f m/s, Rho: %.1f - %.1f\n', ...
    min(vp_f(:)), max(vp_f(:)), min(vs_f(:)), max(vs_f(:)), min(outrh(:)), max(outrh(:)));
fprintf('Final Vp/Vs: %.3f - %.3f (max Vp sets the time step)\n', min(vp_f(:)./max(vs_f(:),eps)), max(vp_f(:)./max(vs_f(:),eps)));

if any(isnan(outla(:))) || any(isnan(outmu(:))) || any(isnan(outrh(:))) || ...
        any(outrh(:) <= 0) || any(outmu(:) <= 0) || any(vp_f(:) <= 0)
    error('NaN, zero, or negative values detected in output arrays! Check interpolation, smoothing, or density relation.');
end

% settings log ----------------------------------------------
if apply_maxdv_clip, clip_word = 'clipped'; else, clip_word = 'not clipped'; end
logtext = {
    sprintf('update_model_smth run: %s', datestr(now,'yyyy-mm-dd HH:MM:SS'));
    sprintf('iteration: %s', ite_nm);
    sprintf('inversion result: damp %g, smooth %g, file %s', sel_damp, sel_smooth, fnm_try);
    sprintf('dv_scale %g; maxdv clip %d (maxdv %g); drhdvp %g; gaussian weight %d (weight_sigma %g m)', dv_scale, apply_maxdv_clip, maxdv, drhdvp, apply_gaussian_weight, weight_sigma);
    sprintf('min_vp_cut %g, min_vs_cut %g (compared with the updated Vs)', min_vp_cut, min_vs_cut);
    sprintf('cells with |dv| > maxdv (%s): Vp %d, Vs %d; max |dv|: Vp %.4f, Vs %.4f; lambda<0: %d', clip_word, nexc_vp, nexc_vs, max(abs(itp_vp(:))), max(abs(itp_vs(:))), numel(indxla));
    sprintf('cells kept at reference: %d (reference Vp<=cut: %d; updated Vs<cut: %d, of which reference Vs already below the cut: %d)', nnz(fixed), nrev_vp, nrev_vs, nrev_vs_refbelow);
    sprintf('smoothing: apply %d, type %s, box [%s], sigma [%s], preserve fixed %d', apply_smoothing, smooth_type, num2str(smooth_box), num2str(smooth_sigma), smooth_preserve_fixed);
    sprintf('extra vertical smoothing applied: %d (type %s, sigma %g cells, box %d cells)', apply_smoothing && do_extra, extra_vsmooth_type, extra_vsmooth_sigma, extra_vsmooth_box);
    sprintf('final Vp %.0f-%.0f, Vs %.0f-%.0f, rho %.1f-%.1f', min(vp_f(:)), max(vp_f(:)), min(vs_f(:)), max(vs_f(:)), min(outrh(:)), max(outrh(:)))};
cellfun(@(s) fprintf('%s\n', s), logtext);

if write_model
%------------------- create output nc file -----------------------
% make sure the updated folder exists and holds media files to write into
has_media = ~isempty(dir(fullfile(dir_media_updated,'media_mpi*.nc')));
if has_media && backup_existing_updated
    dir_backup = [dir_media_updated '_backup_' datestr(now,'yyyymmdd_HHMMSS')];
    fprintf('Backing up existing updated media to %s ...\n', dir_backup);
    copyfile(dir_media_updated, dir_backup);
end
if ~has_media
    if ~exist(dir_media_updated,'dir'), mkdir(dir_media_updated); end
    fprintf('No media files in %s: copying the input media files there first ...\n', dir_media_updated);
    copyfile(fullfile(dir_media,'media_mpi*.nc'), dir_media_updated);
end

disp('write the updated model ...')
for i=0:MPII-1
    for j=0:MPIJ-1
        for k=0:MPIK-1
            fnm_out=fullfile(dir_media_updated, sprintf('media_mpi%02d%02d%02d.nc', i, j, k));
            % Write to NetCDF
            nc_varput(fnm_out,'mu',permute(outmu(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
            nc_varput(fnm_out,'lambda',permute(outla(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
            nc_varput(fnm_out,'rho',permute(outrh(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
        end
    end
end
fid = fopen(fullfile(dir_media_updated, ['update_log_' ite_nm '.txt']), 'w');
fprintf(fid, '%s\n', logtext{:});
fclose(fid);
disp('finished updating the model media')
else
    disp('write_model = 0: nothing written')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions

function out = smooth_masked(F, fixed, type, box, sigma)
% Smoothing in grid cells that ignores the cells flagged in 'fixed' (normalized convolution: smoothed
% values divided by smoothed weights) and puts those cells back exactly as they were. With no fixed cells
% the box option is identical to smooth3(F,'box',box), which pads the edges by replication. The Gaussian
% branch zero-pads and renormalizes at the grid edges.
    w = double(~fixed);
    if strcmp(type, 'gaussian')
        num = gsmooth(F.*w, sigma);
        den = gsmooth(w,    sigma);
    else
        num = smooth3(F.*w, 'box', box);
        den = smooth3(w,    'box', box);
    end
    out = num ./ max(den, eps);
    out(fixed) = F(fixed);
end

function out = gsmooth(F, sigma)
% Separable Gaussian filter, sigma = [s1 s2 s3] in cells (0 = no smoothing along that dimension).
% Core MATLAB only (convn). smooth3 'gaussian' takes a single sigma for all dimensions, so it cannot
% give different lateral and vertical widths.
    out = double(F);
    for d = 1:3
        if sigma(d) > 0
            n = ceil(3*sigma(d));
            g = exp(-0.5*((-n:n)/sigma(d)).^2);  g = g/sum(g);
            shp = [1 1 1];  shp(d) = numel(g);
            out = convn(out, reshape(g, shp), 'same');
        end
    end
end

function report_smoothing(ZS, label, type, box, sigma)
% print the physical vertical size of the smoothing at a few depths, using a column near the grid center.
% Smoothing in grid cells follows the resolution loss with depth only if the vertical cell size grows with depth.
% Depth is measured from the top of the column, assumed to be k = NZ (*).
    [NX, NY, NZ] = size(ZS);
    zc = double(squeeze(ZS(round(NX/2), round(NY/2), :)));
    depth_km = abs(zc(end) - zc)/1000;
    dz = abs(gradient(zc));
    if strcmp(type, 'gaussian')
        fprintf('%s: Gaussian sigma = [%g %g %g] cells\n', label, sigma);
    else
        fprintf('%s: box = [%d %d %d] cells\n', label, box);
    end
    for dd = [5 10 20 30 40]
        [~, k] = min(abs(depth_km - dd));
        if strcmp(type, 'gaussian')
            fprintf('   depth %2d km: dz = %.2f km, vertical FWHM = %.1f km\n', dd, dz(k)/1000, 2.355*sigma(3)*dz(k)/1000);
        else
            h = (box(3)-1)/2;  k1 = max(1, k-h);  k2 = min(NZ, k+h);
            fprintf('   depth %2d km: dz = %.2f km, vertical box = %.1f km\n', dd, dz(k)/1000, sum(dz(k1:k2))/1000);
        end
    end
end
