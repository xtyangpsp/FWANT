% matlab script to update the velocity model using inversion results
% last modified from, n.cascadia/update_model_smth.m, Y.Shen, 03-13-2011
% script to update the velocity model
% Restructured by Xiaotao Yang @UMASS for better user experiences.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up global parameters in this section.

close all; clear all
addpath(genpath('/depot/xtyang/data/codes/mexcdf'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise'))
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))
MFILEROOT='../mfiles';
path([MFILEROOT '/fun-spool'],path);

% current model iteration
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

sel_smooth=8; % model
sel_damp=10; %

% minimum Vp and Vs allowed. model values below these will be reverted back to the input model.
min_vp_cut = 4000; % in meters/second
min_vs_cut = 2500; % in m/s
maxdv = 0.1; %maximum allowed perturbation. 0.3=30% + or -.
sigma = 10000; % Transition zone in meters (adjust based on grid spacing)
apply_gaussian_weight = 0; %if 1, apply distance-based gaussian weight or smoothing using sigma

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
XYZ=zeros(SIZ,3);
ii=0;
for k=1:NZ
    for j=1:NY
        for i=1:NX
            ii=ii+1;
            XYZ(ii,1)=XSIM{n}(i,j,k);
            XYZ(ii,2)=YSIM{n}(i,j,k);
            XYZ(ii,3)=ZSIM{n}(i,j,k);
        end
    end
end

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

% put the inversion block locations in a three colume array
mzplus=mz+1; nb=mx*my*mzplus;
ii=0;
for k=1:mzplus
    for j=1:my
        for i=1:mx
            ii=ii+1;
            xyzb(ii,1)=xb(i,j,k);
            xyzb(ii,2)=yb(i,j,k);
            xyzb(ii,3)=zb(i,j,k);
            vp(ii)=rvp(i,j,k);
            vs(ii)=rvs(i,j,k);
        end
    end
end

% Note: TriScatteredInterp is replaced by scatteredInterpolant for stability
disp('Constructing interpolation function...');
xyzb=double(xyzb); vp=double(vp'); vs=double(vs');

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

% 2. Calculate Gaussian Weighting Mask for boundary stability
if apply_gaussian_weight
	disp('Calculating Gaussian weighting mask...');
% Find distance to nearest inversion point. 
% If knnsearch is missing, use: dist = pdist2(xyzb, XYZ, 'euclidean', 'Smallest', 1);
	[~, dist] = knnsearch(xyzb, XYZ, 'K', 1); 

	weight = exp(-(dist.^2) / (2 * sigma^2));
	weight(dist > 3*sigma) = 0; % Hard cut-off for far-field
	weight_3d = reshape(weight, NX, NY, NZ);

% 4. Apply Gaussian Weighting and Safety Clip
% This ensures a smooth transition and no "cliff" edges
	itp_vp = max(min(itp_vp_3d .* weight_3d, maxdv), -1.0*maxdv);
	itp_vs = max(min(itp_vs_3d .* weight_3d, maxdv), -1.0*maxdv);
else
	%only apply safety clip
	itp_vp = max(min(itp_vp_3d, maxdv), -1.0*maxdv);
	itp_vs = max(min(itp_vs_3d, maxdv), -1.0*maxdv);
end
% 5. Update the model
outvp = mvp .* (1 + itp_vp);
outvs = mvs .* (1 + itp_vs);
drhdvp = 0.297; 
outrh = mrh .* (1 + drhdvp * itp_vp);

% --- CRITICAL STABILITY TRAP ---
% This will stop the script BEFORE it writes the bad files to Gautschi
if any(outvp(:) <= 0) || any(outrh(:) <= 0) || any(isnan(outvp(:)))
    fprintf('!!! STABILITY ERROR DETECTED !!!\n');
    fprintf('Zero or NaN at i=%d, j=%d, k=%d\n', iv, jv, kv);
    error('Aborting: The updated model contains values that will cause dt=0.');
end

clear indx
% finally identify the water layer in the reference model and fix the velocity
indx = find(mvp <= min_vp_cut); 
outvp(indx)=mvp(indx);
outvs(indx)=mvs(indx); 
outrh(indx)=mrh(indx);
clear indx

indx = find(outvs < min_vs_cut); 
outvp(indx)=mvp(indx);
outvs(indx)=mvs(indx);
outrh(indx)=mrh(indx);

% convert vp,vs to mu,la
outmu=outrh.*outvs.^2;
outla=outrh.*outvp.^2-outmu*2;

% find grids with la < 0
indxla=find(outla<0);
outla(indxla)=0; 

if 1
    disp('smoothening the model ...')
    outrhtmp=smooth3(outrh(:,:,:),'box', [7 7 1]);
    outmutmp=smooth3(outmu(:,:,:),'box', [7 7 1]);
    outlatmp=smooth3(outla(:,:,:),'box', [7 7 1]);
    outrh=outrhtmp; outmu=outmutmp; outla=outlatmp;
end

[outx, outy, outz]=size(outvp);
ib=floor(outx/MPII); jb=floor(outy/MPIJ); kb=floor(outz/MPIK);

% FINAL STABILITY CHECK
fprintf('Min Vp: %f, Max Vp: %f\n', min(outvp(:)), max(outvp(:)));
fprintf('Min Rho: %f, Max Rho: %f\n', min(outrh(:)), max(outrh(:)));

if any(isnan(outla(:))) || any(isnan(outmu(:))) || any(isnan(outrh(:)))
    error('NaN detected in output arrays! Check interpolation or density relation.');
end

if 1
%------------------- create output nc file -----------------------
disp('write the updated model ...')
for i=0:MPII-1
    for j=0:MPIJ-1
        for k=0:MPIK-1
            if i<10   
                fnm_out=[dir_media_updated '/media_mpi0',num2str(i),'0',num2str(j),'0',num2str(k),'.nc'];
            else
                fnm_out=[dir_media_updated '/media_mpi',num2str(i),'0',num2str(j),'0',num2str(k),'.nc'];
            end
            % Write to NetCDF
            nc_varput(fnm_out,'mu',permute(outmu(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
            nc_varput(fnm_out,'lambda',permute(outla(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
            nc_varput(fnm_out,'rho',permute(outrh(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
        end
    end
end
disp('finished updating the model media')
end

