% matlab script to update the velocity model using inversion results
% last modified from, n.cascadia/update_model_smth.m, Y.Shen, 03-13-2011
% script to update the velocity model
% Restructured by Xiaotao Yang @UMASS for better user experiences.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up global parameters in this section.

close all; clear all
MFILEROOT='../mfiles';
path([MFILEROOT '/fun-spool'],path);

% current model iteration
ite_nm = 'ite_0.05deg_05';
inv_homedir='inv.structure.kerVpVs';
inv_workdir='Freq12345.Z.model';
inv_modeldir='result.1th';
mx=69; my=159; mz=52; %blockfile dimension
sel_smooth=8;
sel_damp=8;


% set up dimensions for output
MPII=6; MPIJ=10; MPIK=2; % 

% read previous model -----------------------------------
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./input_' ite_nm];
dir_metric=['./input_' ite_nm];
dir_media=['./input_' ite_nm];

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

fnm_try=strcat(inv_modeldirfull,'/','try.damp',num2str(sel_damp),'.smot',...
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
%clear xb yb zb rvp rvs

% set up the functions for interpolation
% Note: TriScatteredInterp is unavailable in matlab versions prior to 2009
disp('construct interpolation function ...');
xyzb=double(xyzb);vp=double(vp');vs=double(vs');
FP=TriScatteredInterp(xyzb,vp);
FS=TriScatteredInterp(xyzb,vs);

% interpolate vp and vs at simulation grids and replace resulting NaN with 0
disp('interpolate velocities at simulation grids ...');
itp_vp=FP(XYZ);itp_vpnan=isnan(itp_vp);
INDNaN=find(itp_vpnan==1);itp_vp(INDNaN)=0;

itp_vs=FS(XYZ);itp_vsnan=isnan(itp_vs);
clear INDNaN; 
INDNaN=find(itp_vsnan==1);itp_vs(INDNaN)=0;

% reshape the interpolated wave speed perturbations in to NX by NY by NZ matrix
itp_vp=reshape(itp_vp,NX,NY,NZ);
itp_vs=reshape(itp_vs,NX,NY,NZ);

% update the model, assuming a linear density-velocity relation (Christensen & Mooney, ?)
outvp=mvp.*(1+itp_vp);
outvs=mvs.*(1+itp_vs);
drhdvp=0.297; % density in kg/m^3, wave speed in m/s)
outrh=mrh.*(1+drhdvp*itp_vp);

% find grids that are unlikely for earth materials
indx=find(outvp.^2 <= 2.*outvs.^2); % lamda > 0 or Vp > ~1.4*Vs (for crust Vp =~1.7-1.8Vs)
alpha=0.3;
outvp(indx)=mvp(indx).*(1+alpha*itp_vp(indx));
outvs(indx)=mvs(indx).*(1+alpha*itp_vs(indx));
outrh(indx)=mrh(indx).*(1+drhdvp*itp_vp(indx));
clear indx
indx=find(outvp.^2 <= 2.*outvs.^2); % if the criterion is still unsatisfied 
outvp(indx)=mvp(indx);
outvs(indx)=mvs(indx);
outrh(indx)=mrh(indx);
clear indx
vpvs=2.2; % upper Vp/Vs limit
indx=find(outvp > vpvs*outvs); 
outvp(indx)=(outvp(indx)+vpvs*outvs(indx))/2;
outvs(indx)=(outvp(indx)/vpvs+outvs(indx))/2;
clear indx
vpvs=1.5; % lower Vp/Vs limit
indx=find(outvp < vpvs*outvs);
outvp(indx)=(outvp(indx)+vpvs*outvs(indx))/2;
outvs(indx)=(outvp(indx)/vpvs+outvs(indx))/2;
clear indx
% finally identify the water layer in the reference model and fix the velocity in the water layer
indx = find(mvp <= 1500); % find water in the initial model
outvp(indx)=mvp(indx);
outvs(indx)=mvs(indx); 
outrh(indx)=mrh(indx);
clear indx
indx = find(outvp < 1450); % no vp less than 1450 m/s in the output
outvp(indx)=mvp(indx);
outvs(indx)=mvs(indx);
outrh(indx)=mrh(indx);

% convert vp,vs to mu,la
outmu=outrh.*outvs.^2;
outla=outrh.*outvp.^2-outmu*2;
% find grids with la < 0
indxla=find(outla<0);
outla(indxla)=0; %clear indxla
smoothscale=3;
if 0
    disp('smoothening the model ...')
    outrhtmp=smooth3(outrh(:,:,:),'box',[smoothscale smoothscale 1]);
    outmutmp=smooth3(outmu(:,:,:),'box',[smoothscale smoothscale 1]);
    outlatmp=smooth3(outla(:,:,:),'box',[smoothscale smoothscale 1]);
    outrh=outrhtmp; outmu=outmutmp; outla=outlatmp;
end

[outx, outy, outz]=size(outvp);
ib=floor(outx/MPII);jb=floor(outy/MPIJ);kb=floor(outz/MPIK);
%%
if 1
%------------------- create output nc file -----------------------
disp('write the updated model ...')
for i=0:MPII-1
    for j=0:MPIJ-1
        for k=0:MPIK-1
            if i<10   
                fnm_out=['./updated_input_' ite_nm '/media_mpi0',num2str(i),'0',num2str(j),'0',num2str(k),'.nc']
            else
                fnm_out=['./updated_input_' ite_nm '/media_mpi',num2str(i),'0',num2str(j),'0',num2str(k),'.nc']
            end
            nc_varput(fnm_out,'mu',permute(outmu(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
            nc_varput(fnm_out,'lambda',permute(outla(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
            nc_varput(fnm_out,'rho',permute(outrh(i*ib+1:(i+1)*ib,j*jb+1:(j+1)*jb,k*kb+1:(k+1)*kb),[3 2 1]),[3 3 3],[kb jb ib],[1 1 1]);
        end
    end
end

disp('finished creating')

end


