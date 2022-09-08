% Script combining Vs_Yang_2008.dat (0-200 km depth) and AK135 at great depth to generate a 
% media volume for finite-difference simulation in the spherical coordinate system
% for the study area out of Yang's model, we use CUB2 model;
% The lat and long dimensions are defined in read_CUB_yang2008.m

clear all; close all;
define_latlon;
%extend range outword for 3 degrees in both x and y directions.
%modified by Xiaotao Yang
minlon=minlon-3; maxlon=maxlon+3; minlat=minlat-3; maxlat=maxlat+3;

%minlon=360-79.5; maxlon=360-68.5; minlat=35; maxlat=46;

origgrid=0.15;

latgrid=[minlat:origgrid:maxlat];
longrid=[minlon:origgrid:maxlon];
[LAT, LON]=meshgrid(latgrid, longrid);

% run ~/FWT/ANT/Proj/cascadia/model_updates/set_netcdf;


%%
ak135=load('./ref_models/ak135.tvel');
[mak nak]=size(ak135);

% CUB from 0-70 km
[lat1 lon1 dep1 vs1]=read_CUB_revised(minlon,maxlon,minlat,maxlat);

%%
for ii=1:length(dep1)
vs(:,:,ii)=griddata(lat1,lon1,squeeze(vs1(:,ii)),LAT,LON);

clear x1 y1 ii0
[x1 y1]=find(isnan(vs(:,:,ii)));
ii0=find(~isnan(vs(:,:,ii)));
if ~isempty(x1) & ~isempty(y1) & ~isempty(ii0) 
  for kk=1:length(x1)
     vs(x1(kk),y1(kk),ii)=mean(vs(ii0));
  end
end
    
end


%%
modelname='./ref_models/NA07_kmps.nc'

lat2=nc_varget(modelname,'latitude');
lon2=360+nc_varget(modelname,'longitude');
dep2=nc_varget(modelname,'depth');
vs2=nc_varget(modelname,'vs');

fid=fopen('NA07_Vs.dat','w');
fprintf(fid,'%s %s %s %s\n','Lat','Lon','Depth','Vs');
for ii=1:size(vs2,3)
    for jj=1:size(vs2,2)
        for kk=1:size(vs2,1)-1
        if lat2(ii)>=minlat & lat2(ii)<=maxlat & lon2(jj)>=minlon & lon2(jj)<=maxlon    
        fprintf(fid,'%f %f %f %f\n', lat2(ii),lon2(jj),dep2(kk),vs2(kk));
        end
        end
    end
end
fclose(fid);

% NA07 from 70-670 km
for ii=1+length(dep1):length(dep1)+length(dep2)-1
vs(:,:,ii)=griddata(lat2,lon2,squeeze(vs2(ii-length(dep1),:,:)),LAT,LON);

clear x1 y1 ii0
[x1 y1]=find(isnan(vs(:,:,ii)));
ii0=find(~isnan(vs(:,:,ii)));
if ~isempty(x1) & ~isempty(y1) & ~isempty(ii0) 
  for kk=1:length(x1)
     vs(x1(kk),y1(kk),ii)=mean(vs(ii0));
  end
end

end

dep = [dep1'; dep2(1:end-1)];
    
lat = latgrid;
lon = longrid;
colat=90-lat;
nlat=numel(lat);
nlon=numel(lon);
ndep=numel(dep);

moho=zeros(nlon,nlat);
%%
%------------------- convert Vs to Vp  -----------
for i=1:nlon
  for j=1:nlat
    for k=1:ndep
      %depth=(k-1)*4;
      depth=dep(k);
      if (vs(i,j,k)<4.0);  % proxy for crust
        vp(i,j,k)=vs(i,j,k)*1.74; 
        moho(i,j)=k;
      else
	% use vp/vs of ak135
        for m=1:mak
          if(ak135(m,1)>depth), break; end;
        end
        ratu=ak135(m-1,2)/ak135(m-1,3);
        ratd=ak135(m,2)/ak135(m,3);
        vp(i,j,k)=vs(i,j,k)*(ratu+(ratd-ratu)/(ak135(m,1)-ak135(m-1,1))*(depth-ak135(m-1,1)));
      end
      % water layer
      if(vs(i,j,k)==0)
	  vp(i,j,k)=1.5;
      end
    end
  end
end

%-------------------calculate density using vp-density relation  -----------
for i=1:nlon
  for j=1:nlat
    % for the mantle (reference?)
    for k=ndep:-1:moho(i,j)+1
      if(k==ndep)
        rho(i,j,k)=3.5068+(vp(i,j,k)-9.0302)*0.384;
      else
        rho(i,j,k)=rho(i,j,k+1)+(vp(i,j,k)-vp(i,j,k+1))*0.384;
      end
    end
    % for the crust (reference?)
    for k=moho(i,j):-1:1
      if(dep(k)<=10),rho(i,j,k)=0.9893+0.2891*vp(i,j,k);
      elseif(dep(k)>10 && dep(k)<=20),rho(i,j,k)=0.9473+0.2966*vp(i,j,k);
      elseif(dep(k)>20 && dep(k)<=30),rho(i,j,k)=0.9466+0.2997*vp(i,j,k);
      elseif(dep(k)>30 && dep(k)<=40),rho(i,j,k)=0.9645+0.3005*vp(i,j,k);
      % ? wrong unit in SE.Tibet?
      %elseif(dep(k)>40),rho(i,j,k)=1078.3+299.0*vp(i,j,k); 
      elseif(dep(k)>40),rho(i,j,k)=1.0783+0.2990*vp(i,j,k); 
      end;
      % water layer
      if(vs(i,j,k)==0); rho(i,j,k)=1.0; end
    end
  end
end

% tmp = find(ak135(:,1)>670); index670=tmp(1);
% tmp = find(ak135(:,1)>1000); index1000=tmp(1);
% 
% for ii=1:(index1000-index670+1)
% ndep=ndep+1;
% dep(ndep)=ak135(index670+ii,1);
% vp(:,:,ndep)=ak135(index670+ii,2);
% vs(:,:,ndep)=ak135(index670+ii,3);
% rho(:,:,ndep)=ak135(index670+ii,4);
% end
%%
% convert to SI units
vp=vp*1000;
vs=vs*1000;
rho=rho*1000;
dep=dep*1000;

vp=permute(vp,[3 1 2]); % depth lon lat
vp=flipdim(vp,3);
vs=permute(vs,[3 1 2]); % depth lon lat
vs=flipdim(vs,3);
rho=permute(rho,[3 1 2]); % depth lon lat
rho=flipdim(rho,3);
colat=fliplr(colat);
%------------------- create output nc file -----------------------

% fnm_out='SeisMedia.volume.eUS.NA07.d500.nc';
fnm_out='../SeisMedia.volume.eUS.NA07.d500.nc';
if 1

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

nc_varput(fnm_out,'theta',colat);
nc_varput(fnm_out,'phi', lon);
nc_varput(fnm_out,'depth2sealevel',dep);

nc_varput(fnm_out,'Vp',vp);
nc_varput(fnm_out,'Vs',vs);
nc_varput(fnm_out,'rho',rho);

disp('finished creating')
end
