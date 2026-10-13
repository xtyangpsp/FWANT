% script to generate the file containing the dimensions of inversion blocks and indexes
% modified from /net/fs01/data/ren/Ambient_Noise_New/inv.structure/inv_make_block_stride_fix.m 
% (also see /net/fs01/data/tibet/ite_02/inv.structure.64x84x29.kerVpVs/inv_make_block_stride_gridcent.m)

clear all, close all

set_netcdf

MFILE_ROOT='../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

fnm_conf=['../sim.station/skel/fx/SeisFD3D.conf']; %linked to sim.station/skel/fx/SeisFD3D.conf
pnm_metric='../sim.input'; %linked to ../sim.input

%%% plot the original 3D node space
id=1; %snap_id=1 is sgt volume; snap_id=2 is velocity on surface 
subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
% get coord data
[X,Y,Z]=gather_coord(snapinfo,'coorddir',pnm_metric);

% read kernel values ONLY where they are calculated (in case kernels are not calculated at 
% every (saved) grids
sub_start = [ 1 1 1];
sub_count = [64 58 48];
sub_stride =[ 1 1 1];

% number of kernel grids per block, and total number of blocks
ii=1; jj=1; kk=1; 
mx=floor(sub_count(1)/ii);
my=floor(sub_count(2)/jj);
mz=floor(sub_count(3)/kk);
volfct=sub_stride(1)*sub_stride(2)*sub_stride(3);

%file naming: block.block_sizes.stride.kernel-per-block
fnm_inv_blk=strcat('block.',num2str(sub_count(1)),'x',num2str(sub_count(2)),...
    'x',num2str(sub_count(3)),'.1x1x1.1x1x1.nc')
%%
% ---------------- read in coordinate -------------------------
% config
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
% get coord data
[X,Y,Z]=gather_coord(snapinfo,'coorddir',pnm_metric);
nx=size(X,1); if nx ~= sub_count(1); disp('nx & sub_count inconsistent');end
ny=size(X,2); if ny ~= sub_count(2); disp('ny & sub_count inconsistent');end
nz=size(X,3); if nz ~= sub_count(3); disp('nz & sub_count inconsistent');end

% ---------------- volume of each cell -------------------------
dz(:,:,2:nz-1)= (Z(:,:,2:nz-1)-Z(:,:,1:nz-2))/2 ...
               +(Z(:,:,3:nz  )-Z(:,:,2:nz-1))/2;
dz(:,:,1)=(Z(:,:,2)-Z(:,:,1))/2; %half height since the next upper grid accounts for half
%half height for the grid on the free surface because the next lower grid accounts for that space
dz(:,:,nz)=(Z(:,:,nz)-Z(:,:,nz-1))/2; 

dx(2:nx-1,:,:)= (X(2:nx-1,:,:)-X(1:nx-2,:,:))/2 ...
                     +(X(3:nx  ,:,:)-X(2:nx-1,:,:))/2;
dx(1,:,:)=(X(2,:,:)-X(1,:,:))/2;
dx(nx,:,:)=(X(nx,:,:)-X(nx-1,:,:))/2;

dy(:,2:ny-1,:)= (Y(:,2:ny-1,:)-Y(:,1:ny-2,:))/2 ...
                     +(Y(:,3:ny  ,:)-Y(:,2:ny-1,:))/2;
dy(:,1,:)=(Y(:,2,:)-Y(:,1,:))/2;
dy(:,ny,:)=(Y(:,ny,:)-Y(:,ny-1,:))/2;

%----- in Cartesian coordinate
%V=dz.*dy.*dx;

%----- in spherical coordinate (CLAT=X; LO=Y; R=Z)
V=Z.^2.*sin(X).*dz.*dy.*dx;

V=V*volfct;

% ---------------- inversion block ----------------------------
nblk=mx*my*mz;
v1=ones(1,ii);
v2=ones(1,jj);
v3=ones(1,kk);
[V1,V2,V3]=meshgrid(v1,v2,v3);
V1=permute(V1,[2 1 3]); V2=permute(V2,[2 1 3]); V3=permute(V3,[2 1 3]);
Vm=V1.*V2.*V3;

% ---------------- grids in block ----------------------------
max_pix=ii*jj*kk;

npix_blk=ones(1,nblk)*max_pix;

indx_blk=zeros(nblk,3);
volm_blk=zeros(nblk,max_pix);
pixs_blk=zeros(nblk,max_pix,3);

x_blk=zeros(1,nblk);
y_blk=zeros(1,nblk);
z_blk=zeros(1,nblk);

% - grid index for one block cell-
I1=[1:sub_stride(1):ii*sub_stride(1)];
J1=[1:sub_stride(2):jj*sub_stride(2)];
K1=[1:sub_stride(3):kk*sub_stride(3)];

[I,J,K]=meshgrid(I1,J1,K1);
I=permute(I,[2 1 3]); J=permute(J,[2 1 3]); K=permute(K,[2 1 3]);
IMAT=[reshape(I,[1,numel(I)]);
      reshape(J,[1,numel(J)]);
      reshape(K,[1,numel(K)])]';

n=0;
for k=1:mz
    k
for j=1:my
for i=1:mx

    n=n+1;
    % index of this block
    indx_blk(n,1:3)=[i,j,k];
    % kernel grid indexes (1D) contained in this block
    I2=I1+(i-1)*ii*sub_stride(1)+sub_start(1)-1;
    J2=J1+(j-1)*jj*sub_stride(2)+sub_start(2)-1;
    K2=K1+(k-1)*kk*sub_stride(3)+sub_start(3)-1;

    pixs_blk(n,1:npix_blk(n),1)=IMAT(:,1)+(i-1)*ii*sub_stride(1)+sub_start(1)-1;
    pixs_blk(n,1:npix_blk(n),2)=IMAT(:,2)+(j-1)*jj*sub_stride(2)+sub_start(2)-1;
    pixs_blk(n,1:npix_blk(n),3)=IMAT(:,3)+(k-1)*kk*sub_stride(3)+sub_start(3)-1;

    volm_blk(n,1:npix_blk(n))=reshape(V(I2,J2,K2).*Vm,[1,npix_blk(n)]);
    x_blk(n)=mean(mean(mean(X(I2,J2,K2))));
    y_blk(n)=mean(mean(mean(Y(I2,J2,K2))));
    z_blk(n)=mean(mean(mean(Z(I2,J2,K2))));
end
end
end

% ---------------- create nc file ----------------------------
if 1
%my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
%nc_create_empty ( fnm_inv_blk, my_mode );
nc_create_empty(fnm_inv_blk);
nc_add_dimension(fnm_inv_blk,'block',nblk);
nc_add_dimension(fnm_inv_blk,'max_pix',max_pix);
nc_add_dimension(fnm_inv_blk,'geometry',3);
nc_add_dimension(fnm_inv_blk,'mx',mx);
nc_add_dimension(fnm_inv_blk,'my',my);
nc_add_dimension(fnm_inv_blk,'mz',mz);

nc_attput(fnm_inv_blk,nc_global,'mx',mx);
nc_attput(fnm_inv_blk,nc_global,'my',my);
nc_attput(fnm_inv_blk,nc_global,'mz',mz);

var.Nctype='float';var.Attribute=[];
var.Name='y';var.Dimension={'block'};nc_addvar(fnm_inv_blk,var);
var.Name='x';var.Dimension={'block'};nc_addvar(fnm_inv_blk,var);
var.Name='z';var.Dimension={'block'};nc_addvar(fnm_inv_blk,var);
var.Name='volume'; var.Dimension={'block','max_pix'};nc_addvar(fnm_inv_blk,var);

var.Nctype='int';var.Attribute=[];
var.Name='num_of_pix';var.Dimension={'block'};nc_addvar(fnm_inv_blk,var);
var.Name='indx_of_pix'; var.Dimension={'block','max_pix','geometry'};nc_addvar(fnm_inv_blk,var);
var.Name='indx_of_block'; var.Dimension={'block','geometry'};nc_addvar(fnm_inv_blk,var);

nc_varput(fnm_inv_blk,'y',y_blk,[0],[nblk],[1]);
nc_varput(fnm_inv_blk,'x',x_blk,[0],[nblk],[1]);
nc_varput(fnm_inv_blk,'z',z_blk,[0],[nblk],[1]);
nc_varput(fnm_inv_blk,'volume',volm_blk);
nc_varput(fnm_inv_blk,'num_of_pix',npix_blk,[0],[nblk],[1]);
nc_varput(fnm_inv_blk,'indx_of_pix',pixs_blk);
nc_varput(fnm_inv_blk,'indx_of_block',indx_blk);

disp(['finished exporting ' fnm_inv_blk])
end
