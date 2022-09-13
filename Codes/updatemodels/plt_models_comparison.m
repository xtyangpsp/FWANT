%
clear all; close all;
% /opt/matlab/2007b/bin/matlab

set_mfiles_path
set_netcdf

ite_nm = ['ite_0.05deg_03'];

idx = 0; % 0, absolute velocity; 1, velocity perturbation

%read previous model
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./updated_input_' ite_nm];
dir_media=['./updated_input_' ite_nm];
disp(['Read current model... ' ite_nm]);

id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];
id{end+1} = 0; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,-1];subt{end+1}=[1,1,1];
               indxem{end+1}=[];
               indxkp{end+1}=[];
n=1;
[snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
[XSIM{n},YSIM{n},ZSIM{n}]=gather_coord(snapinfo{n},'coorddir',dir_coord);
% convert from radian to degrees
XSIM{n}=90-XSIM{n}*180/pi; %latitude
YSIM{n}=YSIM{n}*180/pi;

%define the area of plot (exclude pmls)
npml=12; %number of pml layers
minlat=XSIM{1}(end-npml,1,end);maxlat=XSIM{1}(1+npml,1,end);
minlon=YSIM{1}(1,1+npml,end);maxlon=YSIM{1}(1,end-npml,end);

id{end+1} = 0; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,-1];subt{end+1}=[1,1,1];
               indxem{end+1}=[];
               indxkp{end+1}=[];

mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
mvs{n}=(mmu{n}./mrh{n}).^0.5;
mvs{n}=smooth3(mvs{n},'box',[9 9 1]);

minlon=360-167.5; maxlon=360-129.5; minlat=53.5;maxlat=67.5;

% define vertical grid in the inversion model
nzbgrid=[37 36 35 34 33 32 31 30 29 28 27 26];
% define the corresponding simulation grid
nzsimgrid=17+(nzbgrid-1)*3; %ite_0.025deg


state=[];
load ../../../cascadia/misc/us_states;

% plot body-wave tomography result from Brandon Schmandt
clear LON LAT depint

orggrid = 0.2;   % grid spacing
[LON,LAT] = meshgrid(minlon:orggrid:maxlon, minlat:orggrid:maxlat);

% load result from Brandon Schmandt
ll = 64; % # of grid points at each depth

[lat,lon,Depth, Vs] = textread('US_CrustVs_SLK_GRL_2015.txt','%f %f %f %f');
nlayer = length(Depth)/ll;
lon = lon+360;
Vs = double(Vs);
lat = double(lat);
lon = double(lon);

for zz=1:length(nzbgrid)
    
figure('Position',[400 400 1100 800])

dep=6371-abs(ZSIM{1}(npml,npml,nzsimgrid(zz))/1000); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1),  hold on, box on, axis on
v=squeeze(mvs{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz)));
v=double(v)/1000;
vmeans=mean(mean(v));

cabsmaxs=vmeans*1.1; cabsmins=vmeans*0.9; %vs

pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz))),squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz))),v);
caxis([cabsmins cabsmaxs]);

shading interp;
colormap('jetwr'); 
ch=colorbar('location','southoutside','plotboxaspectratio',[10 0.5 1]);

title(['Vs at ' num2str(round(dep)) ' km'],'FontSize',20);
axis([minlon maxlon minlat maxlat]);
daspect([1 cosd((minlat+maxlat)/2) 1]);
set(gca,'FontSize',12)

for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

set(gca,'XTick',[minlon:2:maxlon],'XTickLabel',[minlon:2:maxlon])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ipoint = find( abs(Depth-dep)==min(abs(Depth-dep)) );

VSgrid=[];
VSgrid = griddata(lon(ipoint),lat(ipoint),Vs(ipoint),LON,LAT,'nearest');
VSgrid=smoothing_hori(VSgrid,3);
%%%

subplot(1,2,2), hold on, box on, axis on
pcolor(LON,LAT,VSgrid);

caxis([cabsmins cabsmaxs]);

shading interp;
colormap('jetwr'); 
ch=colorbar('location','southoutside','plotboxaspectratio',[10 0.5 1]);

title(['Vs at ' num2str(round(dep)) ' km'],'FontSize',20);
axis([minlon maxlon minlat maxlat]);
daspect([1 cosd((minlat+maxlat)/2) 1]);
set(gca,'FontSize',12)


for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

set(gca,'XTick',[minlon:2:maxlon],'XTickLabel',[minlon:2:maxlon]);

drawnow
pause
close all

end

%%%% save figure
% set(gcf,'PaperPositionMode','auto');
% switch idx
%     case 0
% 	figname = ['VelModel_' num2str(dep,3) 'km.eps'];
%     case 1
% 	figname = ['VelModel_Percent_' num2str(dep,3) 'km.eps'];
% end
% eval(['print -depsc ' figname])

