%
clear all; close all;
% /opt/matlab/2007b/bin/matlab

set_mfiles_path

%run ~/FWT/ANT/Proj/cascadia/model_updates/set_netcdf.m

ite_nm = ['sim_0.025deg'];

fnm_conf=['SeisFD3D.conf'];

idx = 0; % 0, absolute velocity; 1, velocity perturbation

id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];

%read updated model (same conf and coord)
dir_media=['input' ];
dir_coord=['input' ];

disp(['Read model... ']);

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

mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
mvs{n}=(mmu{n}./mrh{n}).^0.5;

figure('Position',[200 200 1200 800]);
%figure('Position',[400 400 850 900])

% minlon=235; maxlon=245;minlat=38; maxlat=49;

figwid = 0.3; figheight=0.35;

% define vertical grid in the inversion model
% nzbgrid=[35 25 20 15 10 1]; 

% define the corresponding simulation grid
%nzsimgrid=40+(nzbgrid-1)*3; %OBS_0.05deg_03
nzsimgrid=[134 100 50 1]

for zz=1:length(nzsimgrid)

dep=6371-abs(ZSIM{1}(npml,npml,nzsimgrid(zz))/1000); 

%%%%%%%%% Vs %%%%%%%%%%%%%%%%


subplot(2,2,zz),
hold on, box on, axis on
v=squeeze(mvs{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz)));
v=double(v)/1000;
vmeans=mean(mean(v));
cabsmaxs=vmeans*1.15; cabsmins=vmeans*0.85; %vs

pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz))),squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz))),v);
%caxis([cabsmins cabsmaxs]);
%caxis([5.2 5.6]);
shading interp;
colormap('jetwr'); 
%ch=colorbar('location','eastoutside','plotboxaspectratio',[0.43 6 1]);

title([num2str(round(dep)) ' km'],'FontSize',20);
axis([minlon maxlon minlat maxlat]);
daspect([1 cosd((minlat+maxlat)/2) 1]);
set(gca,'FontSize',14)

set(gca,'XTick',[minlon:3:maxlon],'XTickLabel',[minlon:3:maxlon])
set(gca,'YTick',[minlat:2:maxlat],'YTickLabel',[minlat:2:maxlat])
set(gca,'TickDir','out');

end    



