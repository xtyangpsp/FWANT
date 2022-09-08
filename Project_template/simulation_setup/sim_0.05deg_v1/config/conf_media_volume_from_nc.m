% Script combining CUB2.0 (0-396 km depth) and AK135 at great depth to generate a 
% media volume for finite-difference simulation in the spherical coordinate system
% The lat and long dimensions are defined in read_CUB.m
% Version: modified from the script in ./SE.Tibet/config/structure, Y.S. 03/2010
%
% Xiaotao Yang: extend the model to the depth of 1057 from ak135.
% added reference for the conversion from Vp to density for the crust.
% reference for mantle conversion is still not clear.
clear all; close all;
%------------------- read AK135     -----------
ak135=load('./vmodels/ak135.tvel');
[mak, nak]=size(ak135);

modelfile='vmodels/Chen2016/oiink-VSVPRHO-model.nc';
%
define_latlon;

[lat, lon, dep, vs]=read_netCDF_model3d(modelfile,'vs');
[~, ~, ~, vp]=read_netCDF_model3d(modelfile,'vp');
[~, ~, ~, rho]=read_netCDF_model3d(modelfile,'rho');
figure;
imagesc(lon,lat,vs(:,:,18));set(gca,'ydir','normal');
%%
colat=90-lat;
nlat=numel(lat);
nlon=numel(lon);
ndep=numel(dep);
% 
% 
% ndep=ndep+1;
% dep(ndep)=409;
% vp(:,:,ndep)=9.0302;
% vs(:,:,ndep)=4.8702;
% rho(:,:,ndep)=3.5068;
% 
% ndep=ndep+1;
% dep(ndep)=410;
% vp(:,:,ndep)=9.3601;
% vs(:,:,ndep)=5.0806;
% rho(:,:,ndep)=3.9317;

ndep=ndep+1;
dep(ndep)=460;
vp(:,:,ndep)=9.5208;
vs(:,:,ndep)=5.1864;
rho(:,:,ndep)=3.9273;

ndep=ndep+1;
dep(ndep)=510;
vp(:,:,ndep)=9.6962;
vs(:,:,ndep)=5.2922;
rho(:,:,ndep)=3.9233;

ndep=ndep+1;
dep(ndep)=560;
vp(:,:,ndep)=9.8640;
vs(:,:,ndep)=5.3989;
rho(:,:,ndep)=3.9218;

ndep=ndep+1;
dep(ndep)=610;
vp(:,:,ndep)=10.0320;
vs(:,:,ndep)=5.5047;
rho(:,:,ndep)=3.9206;

ndep=ndep+1;
dep(ndep)=659;
vp(:,:,ndep)=10.2000;
vs(:,:,ndep)=5.6104;
rho(:,:,ndep)=3.9201;

ndep=ndep+1;
dep(ndep)=660;
vp(:,:,ndep)=10.7909;
vs(:,:,ndep)=5.9607;
rho(:,:,ndep)=4.2387;

ndep=ndep+1;
dep(ndep)=710;
vp(:,:,ndep)=10.9222;
vs(:,:,ndep)=6.0898;
rho(:,:,ndep)=4.2986;

ndep=ndep+1;
dep(ndep)=760;
vp(:,:,ndep)=11.0558;
vs(:,:,ndep)=6.2095;
rho(:,:,ndep)=4.4305;

ndep=ndep+1;
dep(ndep)=809.5;
vp(:,:,ndep)=11.1353;
vs(:,:,ndep)=6.2426;
rho(:,:,ndep)=4.4596;

ndep=ndep+1;
dep(ndep)=859;
vp(:,:,ndep)=11.2221;
vs(:,:,ndep)=6.2798;
rho(:,:,ndep)=4.4885;

ndep=ndep+1;
dep(ndep)=908.5;
vp(:,:,ndep)=11.3068;
vs(:,:,ndep)=6.3160;
rho(:,:,ndep)=4.5173;

ndep=ndep+1;
dep(ndep)=958;
vp(:,:,ndep)=11.3896;
vs(:,:,ndep)=6.3512;
rho(:,:,ndep)=4.5459;

ndep=ndep+1;
dep(ndep)=1007.5;
vp(:,:,ndep)=11.4705;
vs(:,:,ndep)=6.3854;
rho(:,:,ndep)=4.5744;

ndep=ndep+1;
dep(ndep)=1057;
vp(:,:,ndep)=11.5495;
vs(:,:,ndep)=6.4187;
rho(:,:,ndep)=4.6028;

% convert to SI units
vp=vp*1000;
vs=vs*1000;
rho=rho*1000;
dep=dep*1000;

%%
vp=permute(vp,[3 1 2]); % depth lon lat
vp=flip(vp,3);
vs=permute(vs,[3 1 2]); % depth lon lat
vs=flip(vs,3);
rho=permute(rho,[3 1 2]); % depth lon lat
rho=flip(rho,3);
colat=fliplr(colat);
%%
figure
imagesc(lon,lat,squeeze(vs(18,:,:)));set(gca,'ydir','normal');
%% ------------------- create output nc file -----------------------
fnm_out='SeisMedia.volume.CUS_Chen2016_AK135.nc';
if 1
    disp('saving to file')
    %my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
    %nc_create_empty ( fnm_out, my_mode );
    nc_create_empty(fnm_out);
    nc_add_dimension(fnm_out,'phi',nlon);
    nc_add_dimension(fnm_out,'theta',nlat);
    nc_add_dimension(fnm_out,'depth',ndep);

    nc_attput(fnm_out,nc_global,'sealevel',6371*1e3);

    var.Nctype='float';var.Attribute=[];
    var.Name='theta';var.Dimension={'theta'};nc_addvar(fnm_out,var);
    var.Name='phi';var.Dimension={'phi'};nc_addvar(fnm_out,var);
    var.Name='depth';var.Dimension={'depth'};nc_addvar(fnm_out,var);
    var.Name='depth2sealevel';var.Dimension={'depth'};nc_addvar(fnm_out,var);

    var.Name='Vp'; var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);
    var.Name='Vs'; var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);
    var.Name='rho';var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);
%     var.Name='Vp'; var.Dimension={'theta','phi','depth'};nc_addvar(fnm_out,var);
%     var.Name='Vs'; var.Dimension={'theta','phi','depth'};nc_addvar(fnm_out,var);
%     var.Name='rho';var.Dimension={'theta','phi','depth'};nc_addvar(fnm_out,var);
    
    nc_varput(fnm_out,'depth2sealevel',dep);
    nc_varput(fnm_out,'phi', lon);
    nc_varput(fnm_out,'theta',colat);
    

    nc_varput(fnm_out,'Vp',vp);
    nc_varput(fnm_out,'Vs',vs);
    nc_varput(fnm_out,'rho',rho);

    disp('finished creating')
end

%%
vs2=ncread(fnm_out,'Vs');
size(vs2)
figure
imagesc(lon,lat,squeeze(vs2(:,:,18))');set(gca,'ydir','normal');