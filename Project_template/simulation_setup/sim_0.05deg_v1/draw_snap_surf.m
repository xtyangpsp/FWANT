% draw_snap_surf: Draw wavefield snapshot by using surf.

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Initial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date: 2008-04-27 17:31:28 -0400 (Sun, 27 Apr 2008) $
% $Revision: 469 $
% $LastChangedBy: zhangw $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

set_mfiles_path
[fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
  dir_station,dir_out]=get_simul_path;

% ----------------------- parameter -----------------------
flag_overlap = 0;

flag_km=0;
flag_emlast=1;
flag_jetwr=1;
flag_light=1;

flag_print = 0;
flag_avi = 1;

% define output component
varnm='Vz'; 
taut=0.5;
scl_daspect=[1 1 1];
%scl_daspect=[10 10 1];
scl_caxis=[-0.05 0.05];

id = 2; 
%id = 1;
subs=[1,1,1];subc=[-1,-1,1];subt=[1,1,1];
%id = 2; subs=[1,1,1];subc=[-1,1,-1];subt=[1,1,1];
%id = 3; subs=[1,1,1];subc=[1,-1,-1];subt=[1,1,1];
%id = 4; subs=[1,1,165];subc=[-1,-1,1];subt=[1,1,1];
%n1=100; n2=5000; dn=100;
%n1=10; n2=n1; dn=1;
n1=1;n2=400;dn=5;% saved every 4 simulation time steps

% -------------------- load data --------------------------

[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
%-- get coord data
%[x,y,z]=gather_coord(snapinfo,'coorddir',dir_coord);
%nx=size(x,1);ny=size(x,2);nz=size(x,3);

    [CLAT,LON,R]=gather_coord(snapinfo,'coorddir',dir_coord);
    nx=size(CLAT,1);ny=size(CLAT,2);nz=size(CLAT,3);
 
    clat=CLAT*180/pi;
    lon=LON*180/pi;
 
    lat=90-clat;
 
%    [y{n},x{n},z{n}]=sph2cart(LON{n},pi/2-CLAT{n},R{n});
  
    y=lon;
    x=lat;
    z=-6371+R/1e3;

str_unit='m';
if flag_km
   x=x/1e3;y=y/1e3;z=z/1e3; str_unit='km';
end

% ----------------------- plot figure -----------------------------------
% -- create new window --
if flag_overlap==1
   hold on
else
   hid=figure;
   set(hid,'BackingStore','on');
   set(hid,'renderer','zbuffer');
   %set(hid,'menubar','none');
   set(gcf, 'PaperPositionMode', 'manual');
   set(gcf,'PaperUnits','points')
   set(gcf,'PaperPosition',[0 0 1024 768])
end

if flag_avi
   aviid = VideoWriter(['snap_' num2str(id,'%3.3i') '_' varnm '.avi']);
   open(aviid);
end

% -- time loop --
lon1d=lon(1,:);
lat1d=clat(:,1);

for nlayer=n1:dn:n2

[v,t]=gather_snap(snapinfo,id,nlayer,varnm,'outdir',dir_out);

v=double(v);

disp([ '  draw ' num2str(nlayer) 'th layer(t=' num2str(t) ')']);

if flag_emlast
%   sid=surf(squeeze(permute(x,[2 1 3])), ...
%            squeeze(permute(y,[2 1 3])), ...
%            squeeze(permute(z,[2 1 3])), ...
%            squeeze(permute(v,[2 1 3])));
%   sid=surf(squeeze(permute(x,[2 1 3])), ...
%            squeeze(permute(y,[2 1 3])), ...
%            squeeze(permute(v,[2 1 3])));
   sid=surf(squeeze(y), ...
            squeeze(x), ...
            squeeze(v));
%    sid=surf(lon1d,lat1d',v);
else
   sid=surf(flipdim(squeeze(permute(x,[2 1 3])),3), ...
            flipdim(squeeze(permute(y,[2 1 3])),3), ...
            flipdim(squeeze(permute(z,[2 1 3])),3), ...
            flipdim(squeeze(permute(v,[2 1 3])),3));
end

%set(sid,'DiffuseStrength',1.0,'SpecularStrength',0.2, ...
%    'SpecularExponent',50, ...
%    'SpecularColorReflectance',0.1)

%axis image
%shading interp;
shading flat;
if exist('scl_caxis'); caxis(scl_caxis); end
% if exist('scl_daspect'); daspect(scl_daspect); end
daspect([1 cosd(mean([x(1,1), x(end,1)])) 1]);
if flag_jetwr; colormap(jetwr); end
if flag_light
   view(0,90)
   set(gca,'box','on');
   camlight(0,10,'local');
   lighting phong
   %camlight
end

colorbar('vert')
%colorbar('location','manual','position',[0.85 0.4 0.068 0.08],'plotboxaspectratiomode','manual','plotboxaspectratio',[1 5 1]...) 
%colorbar('location','manual','position',[0.85 0.1 0.03 0.2])

titlestr=['Snapshot of ' varnm ' at ' ...
          '{\fontsize{16}{\bf ' ...
          num2str(double(t),'%07.3f') ...
          '}}s'];
title(titlestr)

axis tight
%axis([min(x(:,1)) max(x(:,1)) min(y(1,:)) max(y(1,:)) min(min(z)) max(max(z))]);
%axis([min(y(1,:)) max(y(1,:)) min(x(:,1)) max(x(:,1))]);

drawnow
pause(taut);

if flag_print==1
   fnm_out=[varnm '_ndim',num2str(nlayer,'%5.5i')];
   set(gca,'FontName','FixedWidth');
   print(gcf,'-dpng',[fnm_out '.png']);
end

if flag_avi==1
   F = getframe(gca);
    writeVideo(aviid,F);
end

end

if flag_avi==1
   close(aviid);
end

