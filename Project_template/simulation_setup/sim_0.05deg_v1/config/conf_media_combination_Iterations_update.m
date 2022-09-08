% Script combining Vs_Yang_2008.dat (0-200 km depth) and AK135 at great depth to generate a 
% media volume for finite-difference simulation in the spherical coordinate system
% for the study area out of Yang's model, we use CUB2 model;
% The lat and long dimensions are defined in read_CUB_yang2008.m

clear all; close all;

warning off;

run ~/FWT/ANT/Proj/cascadia/model_updates/set_netcdf.m

define_latlon

origgridlat = 0.05;
origgridlon = origgridlat/sind((minlat+maxlat)/2);

lat0=[minlat:origgridlat:maxlat];
lon0=[minlon:origgridlon:maxlon];
[LATgrid, LONgrid]=meshgrid(lat0, lon0);

% load velocity model from the ANT inversion

[lat1 lon1 depth1 vp1 vs1 rho1]=textread('./ref_models/Cascades_VpVs_Gao_ite_0.05deg_06.dat','%f %f %f %f %f %f','headerlines',1);
ndep1 = 132; % 0-1000 km
npoints1 = length(depth1)/ndep1;

intpdepth1 = depth1(1:ndep1);

for ii = 1:ndep1, 
    jj = ii:ndep1:length(depth1);
    VSgrid0(:,:,ii)=griddata(lat1(jj),lon1(jj),vs1(jj),LATgrid,LONgrid);
    VPgrid0(:,:,ii)=griddata(lat1(jj),lon1(jj),vp1(jj),LATgrid,LONgrid);
    
    clear x1 y1 ii0
    [x1 y1]=find(isnan(VSgrid0(:,:,ii)));
    ii0=find(~isnan(VSgrid0(:,:,ii)));
    if ~isempty(x1) & ~isempty(y1) & ~isempty(ii0) 
      for kk=1:length(x1)
      VSgrid0(x1(kk),y1(kk),ii)=mean(VSgrid0(ii0));
      end
    end
    
    clear x1 y1 ii0
    [x1 y1]=find(isnan(VPgrid0(:,:,ii)));
    ii0=find(~isnan(VPgrid0(:,:,ii)));
    if ~isempty(x1) & ~isempty(y1) & ~isempty(ii0) 
      for kk=1:length(x1)
      VPgrid0(x1(kk),y1(kk),ii)=mean(VPgrid0(ii0));
      end
    end
     
end


%%
disp('assign Vs value point by point ...');
for ii=1:length(lat0),
    for jj=1:length(lon0)
        clear idx1 idx2 kk
        idx2=find( abs(LATgrid(1,:)-lat0(ii))==min(abs(LATgrid(1,:)-lat0(ii))) );
        idx1=find( abs(LONgrid(:,1)-lon0(jj))==min(abs(LONgrid(:,1)-lon0(jj))) );
        vs(ii,jj,:)=VSgrid0(idx1,idx2,:);
        vp(ii,jj,:)=VPgrid0(idx1,idx2,:);     
    end
end

disp('Done ...')
%%

dep=intpdepth1;

colat=90-lat0;
nlat=numel(lat0);
nlon=numel(lon0);
ndep=numel(dep);

moho=zeros(nlat,nlon);
%------------------- convert Vs to Vp  -----------
for i=1:nlat
  for j=1:nlon
    for k=1:ndep
      depth=dep(k);
      if (vs(i,j,k)<4.0);  % proxy for crust
        moho(i,j)=k;        
      end
    end
  end
end

%-------------------calculate density using vp-density relation  -----------
for i=1:nlat
  for j=1:nlon
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
      elseif(dep(k)>40),rho(i,j,k)=1.0783+0.2990*vp(i,j,k); 
      end;
      % water layer
      if(vs(i,j,k)==0); rho(i,j,k)=1.0; end
    end
  end
end

% convert to SI units
vp=vp*1000;
vs=vs*1000;
rho=rho*1000;
dep=dep*1000;

vp=permute(vp,[3 2 1]); % depth lon lat
vp=flipdim(vp,3);
vs=permute(vs,[3 2 1]); % depth lon lat
vs=flipdim(vs,3);
rho=permute(rho,[3 2 1]); % depth lon lat
rho=flipdim(rho,3);
colat=fliplr(colat);
%------------------- create output nc file -----------------------

fnm_out='SeisMedia.volume.cascadia.d1000.nc';
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
nc_varput(fnm_out,'phi', lon0);
nc_varput(fnm_out,'depth2sealevel',dep);

nc_varput(fnm_out,'Vp',vp);
nc_varput(fnm_out,'Vs',vs);
nc_varput(fnm_out,'rho',rho);

disp('finished creating')
end


