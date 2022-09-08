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

clear all; close all;

set_mfiles_path
% [fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
%   dir_station,dir_out]=get_simul_path
fnm_conf='./SeisFD3D.conf';

dir_coord='./input';  
dir_metric='./input'; 
dir_media='./input'; 
dir_source='./input'; 
dir_station='./input';
dir_out='./output';
% ----------------------- parameter -----------------------
flag_overlap = 0;

flag_km=1;
flag_emlast=1;
flag_jetwr=1;
flag_light=1;

flag_print = 0;
flag_avi = 1;
% define output component
varnm='Vz'; 
taut=0.25; %pause between each snapshot
scl_daspect=[1 1 1];
scl_caxis=[-.1 .1];

id = 2; % snapshot id (see SeisFD3d.conf)
%surface snapshot
subs=[1,1,1];subc=[256,232,1];subt=[1,1,1];
%e-w crosssection
n1=20;n2=300;dn=1;% saved every 4 simulation time steps

% -------------------- load data --------------------------

[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

 [CLAT,LON,R]=gather_coord(snapinfo,'coorddir',dir_coord);
 nx=size(CLAT,1);ny=size(CLAT,2);nz=size(CLAT,3);

%    [y{n},x{n},z{n}]=sph2cart(LON{n},pi/2-CLAT{n},R{n});
    [y,x,z]=sph2cart(LON,pi/2-CLAT,R);
  
str_unit='m';
if flag_km
   x=x/1e3;y=y/1e3;z=z/1e3; str_unit='km';
end

% ----------------------- plot figure -----------------------------------
% -- create new window --
if flag_overlap==1
   hold on
else
   hid=figure('Position',[600 400 900 900]);
   set(hid,'BackingStore','on');
   set(hid,'renderer','zbuffer');
   %set(hid,'menubar','none');
   set(gcf, 'PaperPositionMode', 'manual');
   set(gcf,'PaperUnits','points')
   set(gcf,'PaperPosition',[0 0 1024 768])
end

% -- time loop --
aviid = VideoWriter(['snap_' num2str(id,'%3.3i') '_' varnm '.avi']);
open(aviid);
for nlayer=n1:dn:n2

[v,t]=gather_snap(snapinfo,id,nlayer,varnm,'outdir',dir_out);

v=double(v);

disp([ '  draw ' num2str(nlayer) 'th layer(t=' num2str(t) ')']);

if flag_emlast
   sid=surf(squeeze(permute(y,[2 1 3])), ...
            squeeze(permute(x,[2 1 3])), ...
            squeeze(permute(z,[2 1 3])), ...
            squeeze(permute(v,[2 1 3])));
else
   sid=surf(flipdim(squeeze(permute(x,[2 1 3])),3), ...
            flipdim(squeeze(permute(y,[2 1 3])),3), ...
            flipdim(squeeze(permute(z,[2 1 3])),3), ...
            flipdim(squeeze(permute(v,[2 1 3])),3));
end

set(sid,'AmbientStrength',1,'DiffuseStrength',0.9,'SpecularStrength',0.0, ...
    'SpecularExponent',1, ...
    'SpecularColorReflectance',0.0);

%axis image
shading interp;
% shading flat;
if exist('scl_caxis'); caxis(scl_caxis); end
if exist('scl_daspect'); daspect(scl_daspect); end

if flag_jetwr; colormap(jetwr); end
if flag_light
    view(0,90)
   set(gca,'box','off');
   %camlight(-80,0,'local');
   camlight(0,-45,'local');
%   camlight(0,0,'infinite');
   lighting phong
%   lighting gouraud
   material dull
end

colorbar('vert')


titlestr=['Snapshot of ' varnm ' at ' ...
          '{\fontsize{16}{\bf ' ...
          num2str(double(t),'%07.3f') ...
          '}}s'];
title(titlestr);


drawnow
pause(taut);

if flag_print==1
   fnm_out=[varnm '_ndim',num2str(nlayer,'%5.5i')];
   set(gca,'FontName','FixedWidth');
   print(gcf,'-dpng',[fnm_out '.png']);
end

F = getframe(gca);
writeVideo(aviid,F);


end
close(aviid);
