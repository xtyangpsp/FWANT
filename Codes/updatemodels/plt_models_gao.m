%
clear all; close all;

set_mfiles_path
set_netcdf

ite_nm = ['ite_0.025deg_03'];

idx = 0; % 0, absolute velocity; 1, velocity perturbation

%read previous model
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./input_' ite_nm];
dir_media=['./input_' ite_nm];
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
minlat=XSIM{1}(end,1,end);maxlat=XSIM{1}(1,1,end);
minlon=YSIM{1}(1,1,end);maxlon=YSIM{1}(1,end,end);

mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
mvs{n}=(mmu{n}./mrh{n}).^0.5;

%read updated model (same conf and coord)
dir_media=['./updated_input_' ite_nm ];
disp(['Read updated model... ']);

id{end+1} = 0; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,-1];subt{end+1}=[1,1,1];
               indxem{end+1}=[];
               indxkp{end+1}=[];
n=2;
[snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
mvs{n}=(mmu{n}./mrh{n}).^0.5;
%mvs{n}=smooth3(mvs{n},'box',[19 19 1]);

% difference between the two models
dmvp=mvp{2}-mvp{1};dmvp=dmvp./mvp{1};
dmvs=mvs{2}-mvs{1};dmvs=dmvs./mvs{1};

% define vertical grid in the inversion model
nzbgrid=31;
NZ=size(ZSIM{1},3); 

% define the corresponding simulation grid
% starting vertical grid + saved every other grid
nzsimgrid=17+(nzbgrid-1)*3;

figure('Position',[200 200 1600 1200])

mv=mvs;
vmean=mean(mean(mv{1}(:,:,nzsimgrid)));
cabsmax=vmean*1.15; cabsmin=vmean*0.85; %vs
crelmin=-12;crelmax=12;

subplot(1,3,1); hold on, box on, axis on

dep=6371-abs(ZSIM{1}(npml,npml,nzsimgrid)/1000); 
v=squeeze(mv{1}(npml:end-npml,npml:end-npml,nzsimgrid));
v=double(v);

switch idx
    case 0
      pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),v);
      caxis([cabsmin cabsmax]);
    case 1
      vpercent = (v-mean(mean(v)))/mean(mean(v))*100;
      pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),vpercent);
      caxis([crelmin crelmax]);
end
shading interp;

colormap('jetwr');colorbar('plotboxaspectratio',[0.5 9 1]);
title(['Initial Model at Depth = ' num2str(dep,3) ' km']);

axis([minlon maxlon minlat maxlat]);
daspect([1 cosd((minlat+maxlat)/2) 1]);

state=[];
load ../../../cascadia/misc/us_states;
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

drawnow
clear v

subplot(1,3,2),  hold on, box on, axis on
v=squeeze(mv{2}(npml:end-npml,npml:end-npml,nzsimgrid));
v=double(v);
%v=smoothing_hori(v,3);
switch idx
    case 0
      pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),v);
      caxis([cabsmin cabsmax]);
    case 1
      vpercent = (v-mean(mean(v)))/mean(mean(v))*100;
      pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),vpercent);
      caxis([crelmin crelmax]);
end
shading interp;
colormap('jetwr'); colorbar('plotboxaspectratio',[0.5 9 1]);
title(['Updated Model at Depth = ' num2str(dep,3) ' km']);
axis([minlon maxlon minlat maxlat]);
daspect([1 cosd((minlat+maxlat)/2) 1]);

state=[];
load ../../../cascadia/misc/us_states;
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end

dmv=dmvs;
subplot(1,3,3), hold on, box on, axis on
v=100*squeeze(dmv(npml:end-npml,npml:end-npml,nzsimgrid));
v=double(v);

pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid)),v);hold on
shading interp;
colormap('jetwr');colorbar('plotboxaspectratio',[0.5 9 1]);
title(['Velocity Anomaly at Depth = ' num2str(dep,3) ' km']);
caxis([crelmin crelmax]);
axis([minlon maxlon minlat maxlat]);
daspect([1 cosd((minlat+maxlat)/2) 1]);

state=[];
load ../../../cascadia/misc/us_states;
for sb=1:length(state)
    plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
end


%%%% save figure
% set(gcf,'PaperPositionMode','auto');   
% switch idx
%     case 0
%         figname = ['VelModel_' num2str(dep,3) 'km.eps'];
%     case 1
%         figname = ['VelModel_Percent_' num2str(dep,3) 'km.eps'];
% end
% eval(['print -depsc ' figname])
%} 
 

