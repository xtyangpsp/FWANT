%
clear all; close all;
% MFILEROOT='../mfiles';
% path([MFILEROOT '/fun-spool'],path);

ite_nm = 'ite_0.015deg_01';
velocitytag='S'; % 'P' for P velocities
idx = 0; % 0, absolute velocity; 1, velocity perturbation
savefigtag=1;

%read previous model
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./input_' ite_nm];
dir_media=['./updated_input_' ite_nm];
disp(['Read updated model... ' ite_nm]);

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
%{
%read updated model (same conf and coord)
% dir_media=['./updated_input_' ite_nm ];
% disp('Read updated model... ');
% 
% id{end+1} = 0; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,-1];subt{end+1}=[1,1,1];
%                indxem{end+1}=[];
%                indxkp{end+1}=[];
% n=2;
% [snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
% mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
% mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
% mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
% mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
% mvs{n}=(mmu{n}./mrh{n}).^0.5;
% %mvs{n}=smooth3(mvs{n},'box',[19 19 1]);
% 
% % difference between the two models
% dmvp=mvp{2}-mvp{1};dmvp=dmvp./mvp{1};
% dmvs=mvs{2}-mvs{1};dmvs=dmvs./mvs{1};
%}
%% plotting
% define vertical grid in the inversion model
NZ=size(ZSIM{1},3); 
% nzbgrid=[25 23 21 19 16 12];
nzbgrid=[25 23 19 14 9];
% define the corresponding simulation grid
% starting vertical grid + saved every other grid
nzsimgrid=73+(nzbgrid-1)*3;
if strcmp(velocitytag,'P')
    mv=mvp;
    %dmv=dmvp;
    figtitletag='Vp';
elseif strcmp(velocitytag,'S')
    mv=mvs;
    figtitletag='Vs';
    %dmv=dmvs;
end
mv{1}=mv{1}/1000; %mv{2}=mv{2}/1000;
%caxisdataabs=[3.4 4; 3.6 4.3;3.75 4.3;3.8 4.65;4.3 4.9; 4.3 4.9];
caxisdataabs=[3.4 4; 3.6 4.3;3.8 4.5;4.4 4.9;4.3 4.8];
caxisdatarel=[-10 10;-10 10;-10 10;-10 10;-10 10; -10 -10];
figure('Position',[400 400 1200 400]);
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) '};
ptLon1 = [284.24, 284.05, 284.16, 289.95];
ptLon2 = [291.34, 289.13, 289.22, 285.03];
ptLat1 = [46.98, 45.83, 43,47.26];
ptLat2 = [44.56, 43.02, 41.9,39.96];

plongridnum = length(ptLon1);

proflabel={'A-A'' ','B-B'' ','C-C'' ','D-D'' '};

for iz=1:length(nzsimgrid)
    vmean=mean(mean(mv{1}(:,:,nzsimgrid(iz))));
%     cabsmax=vmean*1.15; 
%     cabsmin=vmean*0.85; %vs
%     crelmin=-0.01;crelmax=0.01;
    cabsmin=caxisdataabs(iz,1);cabsmax=caxisdataabs(iz,2);
    crelmin=caxisdatarel(iz,1);crelmax=caxisdatarel(iz,2);
    
    subplot(1,5,iz); hold on, box on, axis on

    dep=6371-abs(ZSIM{1}(npml,npml,nzsimgrid(iz))/1000); 
    
    disp(['working on: ', num2str(dep)]);
    
    v=squeeze(mv{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz)));
    v=double(v);
    %v=smoothing_hori(v,3);
    switch idx
        case 0
          pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),...
              squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),v);
          caxis([cabsmin cabsmax]);
        case 1
          vpercent = (v-mean(mean(v)))/mean(mean(v))*100;
          pcolor(squeeze(YSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),...
              squeeze(XSIM{1}(npml:end-npml,npml:end-npml,nzsimgrid(iz))),vpercent);
          caxis([crelmin crelmax]);
    end
    shading interp;
    colormap('jetwr'); 
    c=colorbar('southoutside');
    ax = gca;
    axpos = ax.Position;
    cpos = c.Position;
    cpos(2) = 1.2*cpos(2);
    cpos(4) = 0.5*cpos(4);
    c.Position = cpos;
    ax.Position = axpos;
    title([figlabel{iz},figtitletag, ' at ',num2str(dep,3) ' km']);
    axis([minlon maxlon minlat maxlat]);
    daspect([1 cosd((minlat+maxlat)/2) 1]);
    
    
    if iz==length(nzsimgrid)
        for np=1:plongridnum,
            plot([ptLon1(np) ptLon2(np)],[ptLat1(np) ptLat2(np)],'k-','linewidth',2);
            text(ptLon1(np), ptLat1(np),proflabel{np},'fontsize',14);
        end
    end
    hold off;
    drawnow;
%     if savefigtag
%         fignm=strcat('vmodel_compare_',velocitytag,'_',num2str(int16(dep)));
%         saveas(gca,strcat(fignm,'.png'),'png');
%     end
end 

