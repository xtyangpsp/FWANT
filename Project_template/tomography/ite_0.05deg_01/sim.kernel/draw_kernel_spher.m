% draw_kernel
%
% $Date: 2008-03-03 09:16:59 +0800 (Mon, 03 Mar 2008) $
% $Revision: 2 $
% $LastChangedBy: zhangw $

close all

MFILEROOT='../../mfiles';
path([MFILEROOT '/fun-spool'],path);

% ----------------------- parameter -----------------------
flag_surf=0;
flag_pcolor=1;
flag_overlap = 1;

flag_subplot = 1;
flag_print = 1;
flag_jetwr = 1;
flag_km=1;

% topo
%load /net/fs01/data/yang/easthemi/misc/topogrid.mat
% load ../../misc/topogrid.mat
% 
% topo=qt;
% 
% % coastline
% !cp /home/yang/Proj/Shared/coastline/coastline_1to5m.dat coast.dat
% load coast.dat

% load('us_states.mat');
RUN_ROOT='../sim.station/skel/fx/';
fnm_conf   =[RUN_ROOT '/' 'SeisFD3D.conf'];
dir_coord  =['./input'];
dir_metric = dir_coord;

id = 2; % snapshot id (see SeisFD3d.conf) 
subs=[1,1,1];subc=[-1,-1,1];subt=[1,1,1];
[snapinfosurf]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

[CLATSURF,LONSURF,RSURF]=gather_coord(snapinfosurf,'coorddir',dir_coord);
nxsurf=size(CLATSURF,1);nysurf=size(CLATSURF,2);nzsurf=size(CLATSURF,3);

LATSURF=pi/2-CLATSURF;
clatsurf=CLATSURF*180/pi;
lonsurf=LONSURF*180/pi;
latsurf=90-clatsurf;

%Zin=topo;
Zlat=LATSURF;Zlon=LONSURF;
% define area to be plotted
%flatmin=-45*pi/180;
flatmin=Zlat(end,1);flatmax=Zlat(1,1);
%flonmin=-30*pi/180;
flonmin=Zlon(1,1);flonmax=Zlon(1,end);
nindlat=find(Zlat<flatmin);
nindlon=find(Zlon<flonmin);
Zin(nindlat)=NaN; Zin(nindlon)=NaN;
phi_min=flatmin; phi_max=flatmax;
theta_min=flonmin; theta_max=flonmax;

pband=[7.5 15; 10 20; 15 30; 20 40; 30 60; 40 75; 60 100; 75 125];
fband=flip(flip(1./pband),2);

freq_id=5;
freq_key=['f' num2str(freq_id)];
source='TA.Q38A';
receiver='TA.C35A';
pnm_out=[source,'/',receiver,'/BHZ/',num2str(freq_id),'/T1T2.P2/'];
period_tag=[num2str(1/fband(freq_id,2)),' - ',num2str(1/fband(freq_id,1)),' s'];

% for vertical-vertical cross correlation
% wavename = 'Rayleigh wave';

id=1;
% [x,y,z] order
%subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
subs=[ 1 1 41]; subc=[ -1 -1 1]; subt=[ 1   1  1 ];

var_list=[]; scl_caxis_list=[];
var_list{end+1}='phase_Vp'; scl_caxis_list{end+1}=[-1.5e-14,1.5e-14];
var_list{end+1}='phase_Vs'; scl_caxis_list{end+1}=[-1.5e-13,1.5e-13];
%var_list{end+1}='Kbpz'; scl_caxis_list{end+1}=[-1e-11,1e-11];
%var_list{end+1}='Kbqz'; scl_caxis_list{end+1}=[-2e-11,1e-21];

% -------------------- load data --------------------------

[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

%-- get coord data
[X,Y,Z]=gather_coord(snapinfo,'coorddir',dir_metric);
%X=90-X/pi*180; Y=Y/pi*180; Z=6371e3-Z; Z=Z/1e3;
R=Z; CLAT=X; LON=Y;
nx=size(X,1);ny=size(X,2);nz=size(X,3);

clat=CLAT*180/pi;
lon=LON*180/pi;
lat=90-clat;

[y,x,z]=sph2cart(LON,pi/2-CLAT,R);

str_unit='m';
if flag_km
   x=x/1e3;y=y/1e3;z=z/1e3; str_unit='km';
end

%X=permute(X,[2,1,3]); Y=permute(Y,[2,1,3]); Z=permute(Z,[2,1,3]);
%[x,y,z]=sph2cart(X,Y,Z);
%[x,y,z]=sph2cart(Y,X-pi/2,Z);
%%
% ----------------------- plot kernel -----------------------------------
nvar=length(var_list);
figure('position',[400 400 1000 400])
for n=1:nvar

    [V,varnm]=gather_dist(snapinfo,id,var_list{n},'outdir',pnm_out);
    %v=permute(V,[2,1,3]);

    V=double(V);

    %plot3(cy,cx,cz,'.','Color',[0.2 0.2 0.2],'LineWidth',0.5);
    hold on

    if flag_overlap==1
       hold on
    elseif flag_subplot==0 || n==1
       hid=figure;set(hid,'renderer','zbuffer');
    end
    if flag_subplot==1
       subplot(1,2,n)
    elseif flag_subplot==2
       subplot(4,1,n)
    end

    if flag_surf==1
       sid=surf(squeeze(Y),squeeze(X),squeeze(Z),squeeze(V));
       %View = [-37.5 70];
    else
%        if nx==1
%     %      hid=pcolor(squeeze(Y),squeeze(Z),squeeze(V));
%     %      set(gca,'ydir','reverse');
%           %hid=pcolor(squeeze(Y/pi*180),squeeze(Z),squeeze(v));
%        sid=surf(squeeze(permute(y,[2 1 3])), ...
%                 squeeze(permute(x,[2 1 3])), ...
%                 squeeze(permute(z,[2 1 3])), ...
%                 squeeze(permute(V,[2 1 3])));
%        elseif ny==1
%     %      hid=pcolor(squeeze(X),squeeze(Z),squeeze(V));
%     %      set(gca,'ydir','reverse');
%           %hid=pcolor(squeeze(X/pi*180-90),squeeze(Z),squeeze(v));
%        sid=surf(squeeze(permute(y,[2 1 3])), ...
%                 squeeze(permute(x,[2 1 3])), ...
%                 squeeze(permute(z,[2 1 3])), ...
%                 squeeze(permute(V,[2 1 3])));
%        else
%     %      hid=pcolor(squeeze(Y),squeeze(X),squeeze(V));
%           %hid=pcolor(squeeze(Y/pi*180),squeeze(X/pi*180-90),squeeze(v));
       sid=surf(squeeze(permute(y,[2 1 3])), ...
                squeeze(permute(x,[2 1 3])), ...
                squeeze(permute(z,[2 1 3])), ...
                squeeze(permute(V,[2 1 3])));
%        end
%        view(-90,90)
    end % flag_surf

    %axis image
    shading flat;

    if flag_jetwr==1
       colormap('jetwr');
    end
    %if exist('scl_caxis_list') && ~ isempty(scl_caxis_list)
        caxis(scl_caxis_list{n});
    %end
    xlabel('X');
    view(0,0)
    set(gca,'box','off');
    %camlight(-80,0,'local');
    camlight(0,0,'local');
    %  camlight(0,0,'infinite');
    lighting phong
    %   lighting gouraud
    material dull
    c=colorbar;
    c.Location='eastoutside';
    cpos=c.Position;
    cpos(2)=1.1*cpos(2);
    cpos(3:4)=0.5*cpos(3:4);
    c.Position=cpos;
    %colorbar('location','manual','position',[0.9 0.15 0.02 0.15])

    % depth of horizontal surface in km
    dep=(6371000-R(1,1))/1.e3;
    if strcmp(var_list{n}, 'phase_Vs' )
        phase = 'Vs';
    elseif strcmp(var_list{n}, 'phase_Vp' )
        phase = 'Vp';
    end

    title([phase ', ' num2str(round(dep)) ' km, ' period_tag],'FontSize',14);
    cR=6371*1.e3;
    %sphere3d(Zin,theta_min,theta_max,phi_min,phi_max,cR,4.0*pi/180,'contour','spline',0);
    %Zin=double(Zin);
    theta_min=double(theta_min);theta_max=double(theta_max);
    phi_min=double(phi_min);phi_max=double(phi_max);
    %sphere3d1c(Zin,theta_min,theta_max,phi_min,phi_max,cR,1,'contour','nearest',0);
    %plot3(cy,cx,cz,'--','Color',[0 0 0],'LineWidth',1);
    %view(-10,5);
    %view(0,0);
    %axis off; axis tight;
    grid off;

    drawnow


    
    hold off
end % nvar

if flag_print==1
   print(gcf,'-dpng',['kernel_',source,'_',receiver, '_',freq_key,'_gaussian_a2t6.png']);
end
