%
clear all; close all;
% /opt/matlab/2007b/bin/matlab

set_mfiles_path
set_netcdf

ite_nm = ['ite_0.025deg_02'];

fnm_conf=['./SeisFD3D.conf_' ite_nm];


id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];

%read updated model (same conf and coord)
dir_media=['./updated_input_' ite_nm ];
dir_coord=['./updated_input_' ite_nm ];

disp(['Read updated model... ']);

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
mvs{n}=smooth3(mvs{n},'box',[9 9 1]);

figure('Position',[200 200 1500 1400])

minlon=284; maxlon=292; minlat=39;maxlat=47;

figwid = 0.3; figheight=0.35;

% define vertical grid in the inversion model
nzbgrid=[37 36 35 34 33 32 31 30 29 28 27 26];

% define the corresponding simulation grid
nzsimgrid=17+(nzbgrid-1)*3; %ite_0.025deg

figlabel={'(a)','(b)','(c)','(d)','(e)','(f)'};
for zz=1:length(nzbgrid)

dep=6371-abs(ZSIM{1}(npml,npml,nzsimgrid(zz))/1000); 

subplot(3,4,zz),
hold on, box on, axis on
v=squeeze(mvs{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz)));
v=double(v)/1000;
vmeans=mean(mean(v));

if zz<=6
cabsmaxs=vmeans*1.1; cabsmins=vmeans*0.9; %vs
elseif zz>6
cabsmaxs=vmeans*1.05; cabsmins=vmeans*0.95; %vs
end

pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz))),squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(zz))),v);
caxis([cabsmins cabsmaxs]);

shading interp;
colormap('jetwr'); 
ch=colorbar('location','southoutside','plotboxaspectratio',[10 0.5 1]);

%title([char(figlabel(zz)) ' Vs at ' num2str(round(dep)) ' km'],'FontSize',20);
title(['Vs at ' num2str(round(dep)) ' km'],'FontSize',20);
axis([minlon maxlon minlat maxlat]);
daspect([1 cosd((minlat+maxlat)/2) 1]);
set(gca,'FontSize',12)

state=[];
load ../../../cascadia/misc/us_states;
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

set(gca,'XTick',[minlon:2:maxlon],'XTickLabel',[minlon:2:maxlon])
set(gca,'YTick',[minlat:2:maxlat],'YTickLabel',[minlat:2:maxlat])
set(gca,'TickDir','out');

end    

%%%% save figure
set(gcf,'PaperPositionMode','auto');   
figname = ['VelModel_Vs_' ite_nm '.eps'];
eval(['print -depsc ' figname]);
 
 

