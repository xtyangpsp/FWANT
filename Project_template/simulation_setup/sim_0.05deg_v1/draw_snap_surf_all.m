% draw_snap_surf_all: Draw wavefield snapshot on muliti cross-sections by using surf.

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

flag_km=0;
flag_emlast=0;
flag_jetwr=1;
flag_light=1;

flag_print = 0;
flag_avi = 0;

varnm='Vz'; 
!cat SeisFD3D.conf | grep time_interval_in_s | awk '{print $3}' > taut
load taut;
!/bin/rm taut

scl_daspect=[1 1 150];
scl_caxis=[-0.1 0.1];

!cat SeisFD3D.conf | grep total_grids_in_x | awk '{print $5}' > mx
!cat SeisFD3D.conf | grep total_grids_in_y | awk '{print $5}' > my
!cat SeisFD3D.conf | grep total_grids_in_z | awk '{print $5}' > mz
load mx; load my; load mz;
!/bin/rm mx my mz

% id = 1; % volume output
id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];
id{end+1} = 2; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,1];subt{end+1}=[1,1,1];
indxem{end+1}=[];
indxkp{end+1}=[1,276,1,172,1,1];
%id{end+1} = 1; subs{end+1}=[round(mx/2),1,round(0.5*mz)];subc{end+1}=[1,-1,-1];subt{end+1}=[1,1,1];
%indxem{end+1}=[];
%indxkp{end+1}=[];
n1=1;n2=300;dn=5; % saved every 4 simulation time steps; total 800 saved time steps

nsnap=numel(id);
% -------------------- load data --------------------------
for n=1:nsnap
    [snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
    [CLAT{n},LON{n},R{n}]=gather_coord(snapinfo{n},'coorddir',dir_coord);
    nx{n}=size(CLAT{n},1);ny{n}=size(CLAT{n},2);nz{n}=size(CLAT{n},3);

    clat{n}=CLAT{n}*180/pi;
    lon{n}=LON{n}*180/pi;

    lat{n}=90-clat{n};

    [y{n},x{n},z{n}]=sph2cart(LON{n},pi/2-CLAT{n},R{n});
  
    str_unit='m';
    if flag_km
       x{n}=x{n}/1e3;y{n}=y{n}/1e3;z{n}=z{n}/1e3; str_unit='km';
    end
end

% ----------------------- plot figure -----------------------------------
% -- create new window --
   hid=figure;
   set(hid,'BackingStore','on');
   set(hid,'renderer','zbuffer');
   %set(hid,'menubar','none');
   set(gcf, 'PaperPositionMode', 'manual');
   set(gcf,'PaperUnits','points')
   set(gcf,'PaperPosition',[0 0 1024 768])

if flag_avi
   aviid = avifile(['snap_' num2str(id,'%3.3i') '_' varnm '.avi']);
end

% -- time loop --
for nlayer=n1:dn:n2

for n=1:nsnap
    [v{n},t]=gather_snap(snapinfo{n},id{n},nlayer,varnm,'outdir',dir_out);
    if ~ isempty(indxem{n})
       i1=indxem{n}(1);i2=indxem{n}(2);
       j1=indxem{n}(3);j2=indxem{n}(4);
       k1=indxem{n}(5);k2=indxem{n}(6);
       v{n}(i1:i2,j1:j2,k1:k2)=NaN;
    end
    %debug
    %debug
    if ~ isempty(indxkp{n})
       i1=indxkp{n}(1);i2=indxkp{n}(2);
       j1=indxkp{n}(3);j2=indxkp{n}(4);
       k1=indxkp{n}(5);k2=indxkp{n}(6);
       %debug
       vtmp=v{n}(i1:i2,j1:j2,k1:k2);
       %debug
       v{n}(:,:,:)=NaN;
       %debug
       v{n}(i1:i2,j1:j2,k1:k2)=vtmp;
    end
    disp([ '  draw ' num2str(nlayer) 'th layer(t=' num2str(t) ')']);
    
    % matlab 2009b
    v{n}=double(v{n});

    if flag_emlast
       sid{n}=surf(squeeze(y{n}), ...
		squeeze(x{n}), ...   
                squeeze(z{n}), ...
                   squeeze(v{n}));
    else
       sid{n}=surf(flipdim(squeeze(y{n}),3), ...
                   flipdim(squeeze(x{n}),3), ...
                   flipdim(squeeze(z{n}),3), ...
                   flipdim(squeeze(v{n}),3));
    end
    hold on
    if n>1
    set(sid{n},'DiffuseStrength',1.0,'SpecularStrength',0.2, ...
        'SpecularExponent',50, ...
        'SpecularColorReflectance',0.1)
    end
    hold on
end
hold off

%axis image
%shading interp;
shading flat;
if exist('scl_caxis'); caxis(scl_caxis); end
if exist('scl_daspect'); daspect(scl_daspect); end

if flag_jetwr; colormap(jetwr); end
if flag_light
   view(-30,35)
%   view(80,-60)
   set(gca,'box','off');
   camlight(0,10,'local');
   lighting phong
end
axis tight

%colorbar('vert')

titlestr=['Snapshot of ' varnm ' at ' ...
          '{\fontsize{16}{\bf ' ...
          num2str(double(t),'%07.3f') ...
          '}}s'];
title(titlestr)

drawnow
pause(taut);

if flag_print==1
   fnm_out=[varnm '_ndim',num2str(nlayer,'%5.5i')];
   set(gca,'FontName','FixedWidth');
   print(gcf,'-dpng',[fnm_out '.png']);
end

if flag_avi==1
   F = getframe(gca);
   aviid = addframe(aviid,F);
end

end

if flag_avi==1
   aviid = close(aviid);
end

