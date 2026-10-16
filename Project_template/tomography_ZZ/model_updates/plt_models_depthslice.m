%
clear all; close all;
MFILEROOT='../../../mfiles';
path([MFILEROOT '/fun-spool'],path);
addpath(genpath('/depot/xtyang/data/codes/mexcdf'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise'))
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))

ite_nm = 'ite_0.05deg_05';
velocitytag='S'; % 'P' for P velocities
idx = 0; % 0, absolute velocity; 1, velocity perturbation
savefigtag=1;
visible='on';
%read previous model
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./input_' ite_nm];
dir_media=['./input_' ite_nm];

%read-in key configuration parameters
id_stress = 1; %the snap id for the stress tensor.
confinfo=read_fdconf(fnm_conf,id_stress); 

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
minlon=YSIM{1}(1,1,end)-360;maxlon=YSIM{1}(1,end,end)-360;

mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
mvs{n}=(mmu{n}./mrh{n}).^0.5;

%read updated model (same conf and coord)
dir_media=['./updated_input_' ite_nm ];
disp('Read updated model... ');

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

%% plotting
% minlat=53.5;maxlat=67.5;
% define vertical grid in the inversion model
NZ=size(ZSIM{1},3); 
nzbgrid=[47 46 44 42 40 39 38 37 36 34 33 32 31]; % for Vs
% define the corresponding simulation grid
% starting vertical grid + saved every other grid
nzsimgrid=confinfo.snap_subs(3)+(nzbgrid)*confinfo.snap_subt(3);
if strcmp(velocitytag,'P')
    mv=mvp;
    dmv=dmvp;
elseif strcmp(velocitytag,'S')
    mv=mvs;
    dmv=dmvs;
end
mv{1}=mv{1}/1000; mv{2}=mv{2}/1000;
caxisdataabs=[3.1 3.6; 3.4 4.2;3.5 4.3;3.8 4.2;3.8 4.3;4.1 4.6;...
    4.5 4.8;4.5 4.8;4.6 4.9; 4.6 4.9; 4.6 4.9; 4.6 4.9;...
    4.6 4.9; 4.6 4.9; 4.6 4.9; 4.6 4.9];
% caxisdataabs=[3 5;4.5 6; 5.5 7; 5.5 7;5.5 7; 5.5 7;5.8 7.2;]; %for Vp.
caxisdatarel=[-8 8;-8 8;-8 8;-8 8;-8 8;-8 8;-8 8;-8 8; -8 8;-8 8];
crelmin_diff=-8;
crelmax_diff=8;
% %get Alaska state border
% [alat, alon]=borders('alaska');
state=[];
load us_states;
for iz=1:length(nzsimgrid)
    vmean=mean(mean(mv{1}(:,:,nzsimgrid(iz))));

    figure('Position',[200 200 1250 350],'visible',visible)
    subplot(1,3,1); hold on, box on, axis on

    dep=6371-abs(ZSIM{1}(npml,npml,nzsimgrid(iz))/1000); 
    
    disp(['working on: ', num2str(dep)]);
    
    v=squeeze(mv{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)));
    v=double(v);

    switch idx
        case 0
            cabsmin=caxisdataabs(iz,1);cabsmax=caxisdataabs(iz,2);
            pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)))-360,...
                squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),v);
            caxis([cabsmin cabsmax]);
        case 1
            crelmin=caxisdatarel(iz,1);crelmax=caxisdatarel(iz,2);
            vpercent = (v-mean(mean(v)))/mean(mean(v))*100;
            pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)))-360,...
                squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),vpercent);
            caxis([crelmin crelmax]);
    end
    shading flat;
    colormap('jetwr'); colorbar %('plotboxaspectratio',[0.5 9 1]);
    title(['Initial ',velocitytag, ', ', num2str(dep,3) ' km']);

    axis([minlon maxlon minlat maxlat]);
    daspect([1 cosd((minlat+maxlat)/2) 1]);
%     plot(alon,alat,'k-');
    
    for sb=1:length(state)
        plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[0.3 0.3 0.3],'LineWidth',1);
    end
    set(gca,'TickDir','out')
    hold off;

    drawnow
    clear v

    subplot(1,3,2),  hold on, box on, axis on
    v=squeeze(mv{2}(npml:end-npml,npml:end-npml,nzsimgrid(iz)));
    v=double(v);
    %v=smoothing_hori(v,3);
    switch idx
        case 0
          pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)))-360,...
              squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),v);
          caxis([cabsmin cabsmax]);
        case 1
          vpercent = (v-mean(mean(v)))/mean(mean(v))*100;
          pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)))-360,...
              squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),vpercent);
          caxis([crelmin crelmax]);
    end
    shading flat;
    colormap('jetwr'); colorbar %('plotboxaspectratio',[0.5 9 1]);
    title(['Updated ',velocitytag, ', ',num2str(dep,3) ' km']);
    axis([minlon maxlon minlat maxlat]);
    daspect([1 cosd((minlat+maxlat)/2) 1]);
%     plot(alon,alat,'k-');
    
    
    for sb=1:length(state)
        plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[0.3 0.3 0.3],'LineWidth',1);
    end
    set(gca,'TickDir','out')
    hold off;

    
    subplot(1,3,3), hold on, box on, axis on
    v=100*squeeze(dmv(npml:end-npml,npml:end-npml,nzsimgrid(iz)));
    v=double(v);

    pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)))-360,...
        squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),v);
    shading flat;
    colormap('jetwr');colorbar %('plotboxaspectratio',[0.5 9 1]);
    title([velocitytag, ' anomaly (%), ',num2str(dep,3) ' km']);
    caxis([crelmin_diff crelmax_diff]);
    axis([minlon maxlon minlat maxlat]);
    daspect([1 cosd((minlat+maxlat)/2) 1]);
%     plot(alon,alat,'k-');
    for sb=1:length(state)
        plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[0.3 0.3 0.3],'LineWidth',1);
    end
    set(gca,'TickDir','out')
    hold off;

    %%%% save figure
    set(gcf,'PaperPositionMode','auto');   
    switch idx
        case 0
            figname = [ite_nm '_VelModel_' num2str(dep,3) 'km.png'];
        case 1
            figname = [ite_nm '_VelModel_Percent_' num2str(dep,3) 'km.png'];
    end
    eval(['print -dpng ' figname])
%     pause
%     if savefigtag
%         fignm=strcat('vmodel_compare_',velocitytag,'_',num2str(int16(dep)));
%         saveas(gca,strcat(fignm,'.png'),'png');
%     end
end 

