%
clear all; close all;
MFILEROOT='../mfiles';
path([MFILEROOT '/fun-spool'],path);

ite_nm = 'ite_0.05deg_05';
velocitytag='S'; % 'P' for P velocities
idx = 0; % 0, absolute velocity; 1, velocity perturbation
savefigtag=1;

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
minlat=53.5;maxlat=67.5;
% define vertical grid in the inversion model
NZ=size(ZSIM{1},3); 
% nzbgrid=[28 25 23 21 19 16 12 10 8];
nzbgrid=[51 50 47 44 42 40 38 35 29]; % for Vs
% define the corresponding simulation grid
% starting vertical grid + saved every other grid
nzsimgrid=34+(nzbgrid-1)*2;
if strcmp(velocitytag,'P')
    mv=mvp;
    dmv=dmvp;
elseif strcmp(velocitytag,'S')
    mv=mvs;
    dmv=dmvs;
end
mv{1}=mv{1}/1000; mv{2}=mv{2}/1000;
caxisdataabs=[2.5 3.5;2.5 4; 3.4 4.1;3.75 4.3;4 4.65;4.1 4.7; 4. 4.7;4.1 4.7; 4.3 4.6];
% caxisdataabs=[3 5;4.5 6; 5.5 7; 5.5 7;5.5 7; 5.5 7;5.8 7.2;]; %for Vp.
caxisdatarel=[-8 8;-8 8;-8 8;-8 8;-8 8;-8 8;-8 8;-8 8; -8 8;-8 8];

%get Alaska state border
[alat, alon]=borders('alaska');
for iz=3:length(nzsimgrid)
    vmean=mean(mean(mv{1}(:,:,nzsimgrid(iz))));
%     cabsmax=vmean*1.15; 
%     cabsmin=vmean*0.85; %vs
%     crelmin=-0.01;crelmax=0.01;
    cabsmin=caxisdataabs(iz,1);cabsmax=caxisdataabs(iz,2);
    crelmin=caxisdatarel(iz,1);crelmax=caxisdatarel(iz,2);
    
    figure('Position',[200 200 1250 350])
    subplot(1,3,1); hold on, box on, axis on

    dep=6371-abs(ZSIM{1}(npml,npml,nzsimgrid(iz))/1000); 
    
    disp(['working on: ', num2str(dep)]);
    
    v=squeeze(mv{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)));
    v=double(v);

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
    title(['Initial ',velocitytag, ', ', num2str(dep,3) ' km']);

    axis([minlon maxlon minlat maxlat]);
    daspect([1 cosd((minlat+maxlat)/2) 1]);
    plot(alon,alat,'k-');
    hold off;
    % state=[];
    % load ../../../cascadia/misc/us_states;
    % for sb=1:length(state)
    %     plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
    % end

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
    plot(alon,alat,'k-');
    hold off;
    % state=[];
    % load ../../../cascadia/misc/us_states;
    % for sb=1:length(state)
    %     plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
    % end

    
    subplot(1,3,3), hold on, box on, axis on
    v=100*squeeze(dmv(npml:end-npml,npml:end-npml,nzsimgrid(iz)));
    v=double(v);

    pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)))-360,...
        squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),v);
    shading flat;
    colormap('jetwr');colorbar %('plotboxaspectratio',[0.5 9 1]);
    title([velocitytag, ' anomaly (%), ',num2str(dep,3) ' km']);
    caxis([crelmin crelmax]);
    axis([minlon maxlon minlat maxlat]);
    daspect([1 cosd((minlat+maxlat)/2) 1]);
    plot(alon,alat,'k-');
    hold off;
    % state=[];
    % load ../../../cascadia/misc/us_states;
    % for sb=1:length(state)
    %     plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
    % end


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
    pause
%     if savefigtag
%         fignm=strcat('vmodel_compare_',velocitytag,'_',num2str(int16(dep)));
%         saveas(gca,strcat(fignm,'.png'),'png');
%     end
end 

