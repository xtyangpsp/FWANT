% Matlab script to construct a composite medium with CRUST2.0 for the crust 
% and the uppermost mantle, and  AK135 for deeper (>100 km) mantle.
% Surface topography is ignored.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: This code does not work under matlab 2009, but works fine with matlab 2007.
% The reason is unclear.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

pnm_map='./ref_models/crust2/crust2-maps/';
% crust2 coordinate
vlat=1:2:179; nlat=length(vlat);
vlon=1:2:359; nlon=length(vlon);

% read interfaces
% t0: top of water
% t1-t7: bottom of (water, ice, soft sed, hard sed, upper c., mid c., lower c.=Moho)
fid=fopen([pnm_map '/' 'map2.t0'],'r');
t0=fscanf(fid,'%f',[nlon,nlat]); t0=t0';
fclose(fid);

nlayer=7;
for n=1:nlayer
    idt=fopen([pnm_map '/' 'map2.t' num2str(n)],'r');
    Lb{n}=fscanf(idt,'%f',[nlon,nlat]); Lb{n}=Lb{n}'; % boundary topography
    if n==1, Lh{n}=t0-Lb{n}; else Lh{n}=Lb{n-1}-Lb{n}; end; % layer thickness
    fclose(idt);
end

%%% to account for topography, we adjust the depths of all interfaces
%%% (pushing topography to sea level)
for n=1:nlayer
Lb{n}=Lb{n}-t0;
end
%%%

% read the velocities and densities: 7 crustal layers + mantle beneath Moho
for n=1:nlayer+1
    idvp=fopen([pnm_map '/' 'map2.vp' num2str(n)],'r');
    idvs=fopen([pnm_map '/' 'map2.vs' num2str(n)],'r');
    iddp=fopen([pnm_map '/' 'map2.rho' num2str(n)],'r');
    VpC{n}=fscanf(idvp,'%f',[nlon,nlat]); VpC{n}=VpC{n}';
    VsC{n}=fscanf(idvs,'%f',[nlon,nlat]); VsC{n}=VsC{n}';
    DpC{n}=fscanf(iddp,'%f',[nlon,nlat]); DpC{n}=DpC{n}';
    fclose(idvp);fclose(idvs);fclose(iddp);
end

% define study area: lat1 and lat2 are colatitudes (90-lat)
% define_latlon.m is a matlab snippet to define the output area.
% define_latlon;

% cascadia 
minlon=232.0; maxlon=250.0; minlat=36.0; maxlat=53.0;

lat1=90-maxlat;lat2=90-minlat;lon1=minlon;lon2=maxlon;

% remove values outside the area and assign medium parameters to
% the top and bottom of the layers
j1=find(vlat>=lat1); j1=j1(1); j2=find(vlat<=lat2); j2=j2(end);
i1=find(vlon>=lon1); i1=i1(1); i2=find(vlon<=lon2); i2=i2(end);

vlat([1:j1-1,j2+1:end])=[];
vlon([1:i1-1,i2+1:end])=[];
mlat=length(vlat);
mlon=length(vlon);
for n=1:nlayer
    Lb{n}([1:j1-1,j2+1:end],:)=[]; Lb{n}(:,[1:i1-1,i2+1:end])=[];
    Lh{n}([1:j1-1,j2+1:end],:)=[]; Lh{n}(:,[1:i1-1,i2+1:end])=[];
    Lb{n}=-Lb{n};
    H(n,:,:)=Lh{n};
end
for n=1:nlayer+1
    VpC{n}([1:j1-1,j2+1:end],:)=[]; VpC{n}(:,[1:i1-1,i2+1:end])=[];
    VsC{n}([1:j1-1,j2+1:end],:)=[]; VsC{n}(:,[1:i1-1,i2+1:end])=[];
    DpC{n}([1:j1-1,j2+1:end],:)=[]; DpC{n}(:,[1:i1-1,i2+1:end])=[];
    Vp(n,:,:,1)=VpC{n}; % top of the layer
    Vs(n,:,:,1)=VsC{n};
    Dp(n,:,:,1)=DpC{n};
    Vp(n,:,:,2)=VpC{n}; % bottom of the layer (uniform crustal layers)
    Vs(n,:,:,2)=VsC{n};
    Dp(n,:,:,2)=DpC{n};
end

Vp(n,:,:,1)=8.04;
Vs(n,:,:,1)=4.48;

% test: find and replace Vs==0 with a small positive value
jw=find(Vs==0);Vs(jw)=0.01;

% ------------------------------------------------------------
%% fill the deep part with ak135
ak135=load('./ref_models/ak135.tvel');
[mak nak]=size(ak135);
for i=1:mak
% the bottom of the uppermost mantle layer is the depth of the deepest Moho
%  if(ak135(i,1)>max(max(Lb{nlayer}))),bmoho=i;break;end;
% the bottom of the uppermost mantle layer is at 100 km depth
if(ak135(i,1)>60),bmoho=i;break;end;
end
H(nlayer+1,:,:)=ak135(bmoho,1)-sum(H,1);
Vp(nlayer+1,:,:,2)=ak135(bmoho,2); % replace vel. & den. at bottom of the first mantle layer
Vs(nlayer+1,:,:,2)=ak135(bmoho,3);
Dp(nlayer+1,:,:,2)=ak135(bmoho,4);

% i=23 corresponds to 760 km in ak135.tvel
% i=28 corresponds to 1007 km
% i=33 corresponds to 1255 km
% i=45 corresponds to 1849 km
il=0;
for i=bmoho:28
  if(ak135(i,1)==ak135(i+1,1)),continue;
  else
    il=il+1;
    Vp(nlayer+1+il,:,:,1)=ak135(i,2);
    Vs(nlayer+1+il,:,:,1)=ak135(i,3);
    Dp(nlayer+1+il,:,:,1)=ak135(i,4);
    Vp(nlayer+1+il,:,:,2)=ak135(i+1,2);
    Vs(nlayer+1+il,:,:,2)=ak135(i+1,3);
    Dp(nlayer+1+il,:,:,2)=ak135(i+1,4);
    H(nlayer+1+il,:,:)=ak135(i+1,1)-ak135(i,1);
  end  
end
Vp=permute(Vp,[1 3 2 4]);
Vs=permute(Vs,[1 3 2 4]);
Dp=permute(Dp,[1 3 2 4]);
H=permute(H,[1 3 2]);
Vp=Vp*1000;
Vs=Vs*1000;
Dp=Dp*1000;
H=H*1000;

% Define attenuation parameters
% if Qs > QsINF, then no attenuation effect (see subroutine atten_graves in 
% code/srcF/mod_macdrp.F90)
% QsF0: reference frequency (see Graves, 1996)
% WITHQs
%QsF0=0.1; QsINF=1000.0;
%Qs=ones(size(Vp))*200; % some arbitrary values to test the effect of attenuation

Vp_poly_d =ones(1,nlayer+1+il);
Vs_poly_d =ones(1,nlayer+1+il);
rho_poly_d=ones(1,nlayer+1+il);
% WITHQs
%Qs_poly_d =ones(1,nlayer+1+il);

%%---------------------------- create nc file -------------------------------------
%%fnm_out='SeisMedia.sichuan.sed1.crust2.composite.nc';
fnm_out='SeisMedia.composite.cascadia.crust2.d1000.nc';
if 1
%my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
%nc_create_empty ( fnm_out, my_mode );
nc_create_empty(fnm_out);
nc_add_dimension(fnm_out,'theta',mlat);
nc_add_dimension(fnm_out,'phi',mlon);
nc_add_dimension(fnm_out,'layer',nlayer+1+il);
nc_add_dimension(fnm_out,'side',2);
% WITHQs
%nc_attput(fnm_out,nc_global,'QsF0',QsF0);
%nc_attput(fnm_out,nc_global,'QsINF',QsINF);

var.Nctype='float';var.Attribute=[];
var.Name='theta';var.Dimension={'theta'};nc_addvar(fnm_out,var);
var.Name='phi';var.Dimension={'phi'};nc_addvar(fnm_out,var);
var.Name='Vp_poly_d';var.Dimension={'layer'};nc_addvar(fnm_out,var);
var.Name='Vs_poly_d';var.Dimension={'layer'};nc_addvar(fnm_out,var);
var.Name='rho_poly_d';var.Dimension={'layer'};nc_addvar(fnm_out,var);
% WITHQs
%var.Name='Qs_poly_d';var.Dimension={'layer'};nc_addvar(fnm_out,var);
var.Name='thickness';var.Dimension={'layer','phi','theta'};nc_addvar(fnm_out,var);
var.Name='Vp';  var.Dimension={'layer','phi','theta','side'};nc_addvar(fnm_out,var);
var.Name='Vs';  var.Dimension={'layer','phi','theta','side'};nc_addvar(fnm_out,var);
var.Name='rho'; var.Dimension={'layer','phi','theta','side'};nc_addvar(fnm_out,var);
% WITHQs
%var.Name='Qs'; var.Dimension={'layer','phi','theta','side'};nc_addvar(fnm_out,var);

nc_varput(fnm_out,'theta',vlat);
nc_varput(fnm_out,'phi',vlon);
nc_varput(fnm_out,'Vp_poly_d',Vp_poly_d);
nc_varput(fnm_out,'Vs_poly_d',Vs_poly_d);
nc_varput(fnm_out,'rho_poly_d',rho_poly_d);
% WITHQs
%nc_varput(fnm_out,'Qs_poly_d',Qs_poly_d);

nc_varput(fnm_out,'thickness',H);
nc_varput(fnm_out,'Vp',Vp);
nc_varput(fnm_out,'Vs',Vs);
nc_varput(fnm_out,'rho',Dp);
% WITHQs
%nc_varput(fnm_out,'Qs',Qs);

disp('finished creating')
end
%
