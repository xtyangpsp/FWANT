%This script converts velocity model from text file to nc file used as
%input of the simulation program.
%FORMAT OF INPUT FILE:
%lat lon depth Vp Vs Density
% 48.000 283.500 0.00 2.19 1.26 1.63 
% 48.000 283.500 2.75 4.96 2.85 2.43 
% 48.000 283.500 5.50 6.31 3.63 2.81 
% 48.000 283.500 8.25 6.31 3.63 2.81 

%Xiaotao Yang @ UMASS Amherst
%11/23/2016
%

clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%set up global parameters. %%%%%%%%%%%%%%
modelfile='./ref_models/NewEngland_VpVs_Gao_ite_0.025deg_02.dat';
fnm_out='../SeisMedia.volume.newengland.Gao.d500_0.15.nc';

define_latlon;
gridextend=1.5; 
%extend range outword for 'gridextend' degrees in both x and y directions.
%modified by Xiaotao Yang
minlon=minlon-gridextend; maxlon=maxlon+gridextend; minlat=minlat-gridextend; maxlat=maxlat+gridextend;

%minlon=360-79.5; maxlon=360-68.5; minlat=35; maxlat=46;

origgrid=0.15;
%%%%%%%%%%%%%%%%%%%%%%%% END of setting up global parameters. %%%%%%%%%%%%%%

latgrid=minlat:origgrid:maxlat;
longrid=minlon:origgrid:maxlon;
[LAT, LON]=meshgrid(latgrid, longrid);

% run ~/FWT/ANT/Proj/cascadia/model_updates/set_netcdf;

vmodel=load(modelfile);
%FORMAT OF INPUT FILE:
%lat lon depth Vp Vs Density

% lat0=vmodel(:,1);lon0=vmodel(:,2);dep0=vmodel(:,3);
% vp0=vmodel(:,4);vs0=vmodel(:,5);rho0=vmodel(:,6);

%sort and group by depth first, then lat, then lon
vmodelsort=sortrows(vmodel,[3,1,2]);
[nr,nc]=size(vmodel);
ngroup=0;
deptemp=-99999.0;
j=1;
for i=1:nr
    if(vmodelsort(i,3)==deptemp)
        if(ngroup==1)
            lat0(j)=vmodelsort(i,1);
            lon0(j)=vmodelsort(i,2);
        end
        dep0(ngroup)=vmodelsort(i,3);
        vp0(j,ngroup)=vmodelsort(i,4);
        vs0(j,ngroup)=vmodelsort(i,5);
        rho0(j,ngroup)=vmodelsort(i,6);
        j=j+1;  
        if(i==nr), disp(['layer#, npoints: ',num2str(ngroup),', ',num2str(j-1)]),end
    else
        if(ngroup>=1), disp(['layer#, npoints: ',num2str(ngroup),', ',num2str(j-1)]),end
        ngroup=ngroup+1;
        j=1;
    end
    deptemp=vmodelsort(i,3);
end        

for ii=1:length(dep0)
    
    disp(['Layer: ',num2str(ii)]);
    FFvp=scatteredInterpolant(lat0',lon0',vp0(:,ii));
    vp(:,:,ii)=FFvp(LAT, LON);
    
    FFvs=scatteredInterpolant(lat0',lon0',vs0(:,ii));
    vs(:,:,ii)=FFvs(LAT, LON);
    
    FFrho=scatteredInterpolant(lat0',lon0',rho0(:,ii));
    rho(:,:,ii)=FFrho(LAT, LON);

end

% convert to SI units
vp=vp*1000;
vs=vs*1000;
rho=rho*1000;
dep=dep0*1000;

vp=permute(vp,[3 1 2]); % depth lon lat
vp=flipdim(vp,3);
vs=permute(vs,[3 1 2]); % depth lon lat
vs=flipdim(vs,3);
rho=permute(rho,[3 1 2]); % depth lon lat
rho=flipdim(rho,3);
lat=latgrid;
lon=longrid;

colat=fliplr(90.0-lat);
%%
%------------------- create output nc file -----------------------

if 1
    nc_create_empty(fnm_out);
    nc_add_dimension(fnm_out,'phi',length(lon));
    nc_add_dimension(fnm_out,'theta',length(lat));
    nc_add_dimension(fnm_out,'depth',length(dep));

    nc_attput(fnm_out,nc_global,'sealevel',6371*1e3);

    var.Nctype='float';var.Attribute=[];
    var.Name='theta';var.Dimension={'theta'};nc_addvar(fnm_out,var);
    var.Name='phi';var.Dimension={'phi'};nc_addvar(fnm_out,var);
    var.Name='depth';var.Dimension={'depth'};nc_addvar(fnm_out,var);
    var.Name='depth2sealevel';var.Dimension={'depth'};nc_addvar(fnm_out,var);

    var.Name='Vp'; var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);
    var.Name='Vs'; var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);
    var.Name='rho';var.Dimension={'depth','phi','theta'};nc_addvar(fnm_out,var);

    nc_varput(fnm_out,'theta',colat);
    nc_varput(fnm_out,'phi', lon);
    nc_varput(fnm_out,'depth2sealevel',dep);

    nc_varput(fnm_out,'Vp',vp);
    nc_varput(fnm_out,'Vs',vs);
    nc_varput(fnm_out,'rho',rho);

    disp('finished creating')
end
