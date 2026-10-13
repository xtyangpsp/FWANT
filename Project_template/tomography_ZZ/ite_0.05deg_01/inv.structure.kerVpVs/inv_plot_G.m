%check if the assembled kernel for the inversion blocks is OK
%modified from the same name file in /net/fs01/tibet/ite_02/inv.structure
clear all; close all

MFILE_ROOT='../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

%fnm_blk='block.130x170x60.1x1x1.1x1x1.fixvolm.nc';
%fnm_blk='block.64x84x29.stride.1x1x1.cell.2x2x2.nc';
%fnm_blk='block.49x45x26.1x1x1.2x2x1.nc';
fnm_blk='block.50x58x20.1x1x1.1x1x1.nc';
%mx=130; my=170; mz=60;
%mx=49; my=45; mz=26;
mx=50; my=58; mz=20;
% see inv_make_block_stride_xygridcent.m
%CLA=nc_varget(fnm_blk,'colatitude');
%LO=nc_varget(fnm_blk,'longitude');
%R=nc_varget(fnm_blk,'radius');
x=nc_varget(fnm_blk,'x');
y=nc_varget(fnm_blk,'y');
z=nc_varget(fnm_blk,'z');

%----- in Cartesian coordinate
%X=reshape(x,[mx,my,mz])/1.e3; %from m to km
%Y=reshape(y,[mx,my,mz])/1.e3;
%Z=reshape(z,[mx,my,mz])/1.e3;

%----- in spherical coordinate (X=lat; Y=lon; D=depth)
X=90-reshape(x,[mx,my,mz])/pi*180;
Y=reshape(y,[mx,my,mz])/pi*180;
Z=(6371-reshape(z,[mx,my,mz])/1e3);

size(X)
size(Y)
size(Z)

%fnm_G='G_spool/IC.KMI_2002.180.06.54.42_BHT_1_T3T4.P1.Kap.unformat'
%fnm_G='G_spool/IC.KMI_2002.180.06.54.42_BHT_1_T3T4.P1.Kbp.unformat'
%fnm_G='G_spool/XE.ES06_2003.220.09.59.24_BHZ_1_T3T4.P1.Kbp.unformat'
%fnm_G='G_spool/inv_G_sum.Kap.unformat'
%fnm_G='./G_spool/YA.MC04_2004.087.18.47.29_BHT_2_T1T2.P1.Kbp.unformat'
%fnm_G='./G_spool/XE.ES36_2003.289.12.28.11_BHR_2_T3T4.P1.Kbp.unformat'
%fnm_G='./G_spool/YA.MC07_2004.087.18.47.29_BHZ_2_T5T6.P1.Kbp.unformat'
%fnm_G='./G_spool/YA.MC07_2004.087.18.47.29_BHZ_2_T3T4.P1.Kap.unformat'
%fnm_G='G_spool/IC.KMI_2003.202.15.16.37_BHZ_2_T3T4.P1.Kap.unformat'
% sensitivities to Vp
%fnm_G='./G_spool/TA.E07A_TA.B04A_BHZ_2_T1T2.P2.Kap.unformat'
% sensitivities to Vs
%fnm_G='./G_spool/US.DUG_CN.BMSB_BHZ_2_T1T2.P2.Kbp.unformat'
%fnm_G='./G_spool/US.HAWA_CN.BMSB_BHZ_2_T1T2.P2.Kbp.unformat'
fnm_G='./G_spool/TA.ELFS_CN.EDB_BHZ_2_T1T2.P2.Kbp.unformat'
%fnm_G='./G_spool/TA.F06A_CN.EDB_BHZ_2_T1T2.P2.Kbp.unformat'

% ----------------------------------
fid=fopen(fnm_G,'r')

% n and v0
t=fread(fid,1,'float32');
num_blk=fread(fid,1,'int');
kv0=fread(fid,1,'float32');
t=fread(fid,1,'float32');

t=fread(fid,1,'float32');
npix=fread(fid,1,'int');
t=fread(fid,1,'float32');

% indx and kernel
t=fread(fid,1,'float32');
pos=ftell(fid);
indx=fread(fid,npix,'int',4);
stat=fseek(fid,pos+4,'bof');
ker=fread(fid,npix,'float32',4);

%for n=1:npix
%    indx(n)=fread(fid,1,'int');
%    ker(n)=fread(fid,1,'float32');
%end

fclose(fid);
% ----------------------------------

V=zeros(mx,my,mz);
V1=zeros(1,mx*my*mz);
for n=1:npix
%i=mod(indx(n)-1,mx)+1;
%j=ceil((mod(indx(n)-1,mx*my)+1)/mx);
%k=ceil(indx(n)/(mx*my));
%V(i,j,k)=ker(n);
V1(indx(n))=ker(n);
%pause
end
V=reshape(V1,[mx,my,mz]);
%V=V*1.0e-13;

cmax=5e-2; % n times the cut off value
%cmax=1e-4;
figure('Position',[100 200 600 1100]),

	%cmax=20e-17
	% x slice
	nxslice=10;
	if 1
	subplot(3,1,1), hold on, box on
	set(gcf,'renderer','zbuffer');
	W=squeeze(V(nxslice,:,:))'; %W(end+1,:)=W(end,:);
	YS=squeeze(Y(nxslice,:,:))';
	ZS=squeeze(Z(nxslice,:,:))';
	%pcolor(squeeze(V(25,:,:))');
	%pcolor(W);
	pcolor(YS,ZS,W);
	colormap('jetwr');colorbar;shading flat;
	caxis([-cmax cmax]);
%caxis([-1 1])

	%%V1=V; V1(:,:,46:end)=2.0*V1(:,:,46:end);
	%%W=squeeze(V(25,:,:))'; W(end+1,:)=W(end,:);
	%W(100:end,:)=2.0*W(100:end,:);
	%W(end-1,:)=2.0*W(end-1,:);
	%figure;set(gcf,'renderer','zbuffer');
	%%pcolor(squeeze(V1(25,:,:))');
	%pcolor(W);
	end

	% y slice
	if 1
	clear W ZS YS
	nyslice=10;
	%nsy=22;
	%nsy=45;
	%nsy=100;
	%nsy=19;

	subplot(3,1,2), hold on, box on
	set(gcf,'renderer','zbuffer');
	W=squeeze(V(:,nyslice,:))'; %W(end+1,:)=W(end,:);
	XS=squeeze(X(:,nyslice,:))'; %XS(end+1,:)=XS(end,:);
	ZS=squeeze(Z(:,nyslice,:))'; %ZS(end+1,:)=ZS(end,:);
	%imagesc(W);
	pcolor(XS,ZS,W);
	colorbar
	colormap('jetwr');shading flat
	%set(gca,'ydir','reverse');
	caxis([-cmax cmax])
%caxis([-1 1])

	%W(46:end,:)=2.0*W(46:end,:);
	%W(end-1,:)=2.0*W(end-1,:);
	%figure;set(gcf,'renderer','zbuffer');
	%pcolor(W);
	end

	% z slice
	if 1
	clear W XS
	nzslice=15;
	subplot(3,1,3), hold on, box on
	set(gcf,'renderer','zbuffer');
	%W=squeeze(V(:,:,end));
	W=squeeze(V(:,:,nzslice));
	XS=squeeze(X(:,:,nzslice));
	YS=squeeze(Y(:,:,nzslice));
	%W=squeeze(V(:,:,end-4));
	%W=squeeze(V(:,:,end-10));
	%W=squeeze(V(:,:,1));
	%imagesc(W);
	%pcolor(XS,YS,W);
	pcolor(YS,XS,W);
	colormap('jetwr');shading flat;
	caxis([-cmax cmax])
%caxis([-1 1])
	colorbar
	title(['Vs, at ' num2str(squeeze(Z(1,1,nzslice))) ' km depth'] )
	daspect([1 0.75 1])
    axis([233 251 39 51])

	% station location
	%vslat=47.923;vslon=239.1056;
	%reclat=44.6653; reclon=242.3357;
	%plot(reclon,reclat,'ks','MarkerSize',7,'LineWidth',2,'MarkerFaceColor','w');
	%plot(vslon,vslat,'ks','MarkerSize',8,'LineWidth',2,'MarkerFaceColor','m');
	end

	%Vd=load([fnm_G '.txt']);
	%Vd3=reshape(Vd',[mx my mz]);
	%figure;set(gcf,'renderer','zbuffer');
	%W=squeeze(Vd3(:,:,end));
	%imagesc(W);
	%
	%Vk=load([fnm_G '.ker']);
	%Vk3=reshape(Vk',[mx my mz]);
	%figure;set(gcf,'renderer','zbuffer');
	%W=squeeze(Vk3(:,:,end));
	%imagesc(W);
%caxis([-4 4]*1e-15*1e7)
	%
	%Vv=load([fnm_G '.vol']);
	%Vv3=reshape(Vv',[mx my mz]);
	%figure;set(gcf,'renderer','zbuffer');
	%W=squeeze(Vv3(:,:,end));
	%imagesc(W);
	%
	%Vkv3=Vk3.*Vv3;
	%figure;set(gcf,'renderer','zbuffer');
	%W=squeeze(Vkv3(:,:,end));
	%imagesc(W);


