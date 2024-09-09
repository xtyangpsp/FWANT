%
% clear all; 
close all;
% MFILEROOT='../mfiles';
% path([MFILEROOT '/fun-spool'],path);

ite_nm = 'ite_0.05deg_11';
velocitytag='S'; % 'P' for P velocities, 'PR' for Poisson's ratios, 'RHO' for density.
iplot = 0; % 0, absolute velocity; 1, velocity perturbation
% savefigtag=1;

%read previous model
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./updated_input_' ite_nm];
dir_media=['./updated_input_' ite_nm];
disp(['Read model... ' dir_media]);

id = 0; subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
[XSIM,YSIM,ZSIM]=gather_coord(snapinfo,'coorddir',dir_coord);
% convert from radian to degrees
XSIM=90-XSIM*180/pi; %latitude
YSIM=YSIM*180/pi;

%define the area of plot (exclude pmls)
npml=12; %number of pml layers
% minlat=XSIM(end,1,end);maxlat=XSIM(1,1,end);
% minlon=YSIM(1,1,end);maxlon=YSIM(1,end,end);

mrh=gather_media(snapinfo,'rho','mediadir',dir_media);
mmu=gather_media(snapinfo,'mu','mediadir',dir_media);
mla=gather_media(snapinfo,'lambda','mediadir',dir_media);
mvp=((mla+2*mmu)./mrh).^0.5;
mvs=(mmu./mrh).^0.5;

%% 
if strcmp(velocitytag,'P')
    mv=mvp/1000;
    figtitletag='Vp';
%     mv=mv/1000;
elseif strcmp(velocitytag,'S')
    mv=mvs/1000;
    figtitletag='Vs';
elseif strcmp(velocitytag,'PR')
    mv=(0.5*(mvp./mvs).^2 - 1)./((mvp./mvs).^2 - 1);
    figtitletag='Poisson''s Ratio';
elseif strcmp(velocitytag,'RHO')
    mv=mrh/1000;
    figtitletag='Density';
end
 %mv{2}=mv{2}/1000;

load ../00masterdatafiles/us_states.mat;

%%
slicetype='grid'; 
raycutoff=num2str(10);
%grid OR layer;
% grid: extracts at the specified grid depth; you need to check the grid
% depths first to find the grid depth closest to your expectations.
% layer: velocity between two depth values (interpolates between
% them).
nlayers=0;
if strcmp(velocitytag,'S')
    if strcmp(slicetype,'layer')
        zgridtoplist=[50];
        zgridbotlist=[120];
        nlayers=length(zgridtoplist);
    elseif strcmp(slicetype,'grid')
        zgridlist=[13,22,31,40,51,63,74,97,109,134,165,199];
        nlayers=length(zgridlist);
    else
        error('Wrong layer type, has to be: grid OR average.');
    end
    % ray coverlist for each layer
    raycoverlist={...
        strcat('RayCoverOutline_',ite_nm,'_10-20s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_15-30s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_20-40s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_20-40s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_30-60s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_40-75s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_40-75s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_60-100s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_60-100s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_75-125s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_75-125s_cutoff',raycutoff,'.mat'),...
        strcat('RayCoverOutline_',ite_nm,'_75-125s_cutoff',raycutoff,'.mat')};
    smoothlistx=[9,11,11,13,13,15,15,17,17,19,21,23];
%     smoothlistx=[17,17,19,21,23];
    smoothlisty=smoothlistx; %[15,15,19,23,23];
    smoothscale=smoothlistx.*0.05*110; %[70,70,110,170,250, 320];
elseif strcmp(velocitytag,'P')
    if strcmp(slicetype,'layer')
        zgridtoplist=[5];
        zgridbotlist=[15];
        nlayers=length(zgridtoplist);
    elseif strcmp(slicetype,'grid')
        zgridlist=[9];
        nlayers=length(zgridlist);
    else
        error('Wrong layer type, has to be: grid OR average.');
    end
    % ray coverlist for each layer
    raycoverlist={'AlaskaRayCoverOutline_ite_0.05deg_05_15-35s_cutoff10.mat'};
    smoothlistx=[7];
    smoothlisty=[15];
    smoothscale=smoothlistx.*0.05*110; %[70,70,110,170,250, 320];
end

%%
mvsmoothed=cell(nlayers,1);
disp('Smoothing model ...');
for i=1:nlayers
    mvsmoothed{i}.v=smooth3(mv,'box',[smoothlistx(i) smoothlisty(i) 1]);
    mvsmoothed{i}.tag=num2str(smoothlistx(i));
end
%mv=smooth3(mv,'box',[31 31 1]);

%% get layer velocity
griddepth=squeeze(ZSIM(1,1,end)/1000-abs(ZSIM(npml,npml,:))/1000);
vout=cell(nlayers,1);
if strcmp(slicetype,'layer')
    dz=1.0;
    operation=@mean;
    clear Xgrid Ygrid Zgrid;
    [Ygrid,Xgrid,Zgrid]=meshgrid(squeeze(YSIM(1,:,1)),squeeze(XSIM(:,1,1)),griddepth);
    
    % layerv=cell(length(zgridtoplist),1);
    for iz=1:nlayers
        disp(['Extracting layer: ',num2str(zgridtoplist(iz)),' km to ',num2str(zgridbotlist(iz)),' km ...']);
        vout{iz}.data=velocity_layeroperator(Xgrid,Ygrid,Zgrid,mvsmoothed{iz}.v,operation,'flat',...
            zgridtoplist(iz):dz:zgridbotlist(iz));
        vout{iz}.tag=[num2str(zgridtoplist(iz)),' - ',num2str(zgridbotlist(iz)),' km'];
        vout{iz}.depth=[zgridtoplist(iz), zgridbotlist(iz)]; %depth range: [min, max].
    end
elseif strcmp(slicetype,'grid')
    for iz=1:nlayers
        clear idxv;
        [~,idxv]=nanmin(abs(griddepth - zgridlist(iz)));
        vout{iz}.data=squeeze(mvsmoothed{iz}.v(:,:,idxv));
        vout{iz}.tag=[num2str(int16(griddepth(idxv))),' km'];
        vout{iz}.depth=[int16(griddepth(idxv)),int16(griddepth(idxv))];
    end
else
    error('Wrong layer type, has to be: grid OR average.');
end
%% example of extracting layer between surfaces
% layergrid={zgrid10.data,zgridmoho.data,dz};
% vout_crust=velocity_layeroperator(Xgrid,Ygrid,Zgrid,mvsmoothed{iz}.v,@mean,'surface',...
%     layergrid);
% pcolor(YSIM(1,:,1)-360,XSIM(:,1,1),vout_test);shading flat;
%%
savefilenm=['LayeredVelocity_griddepth_' velocitytag '_' ite_nm '.mat'];
if strcmp(slicetype,'layer')
    save(savefilenm,'zgridtoplist', 'zgridbotlist','vout','smoothlistx','smoothlisty','XSIM','YSIM','ZSIM');
elseif strcmp(slicetype,'grid')
    save(savefilenm,'zgridlist','vout','smoothlistx','smoothlisty','XSIM','YSIM','ZSIM');
else
    error('Wrong layer type, has to be: grid OR average.');
end    
%% load, use this when restart after MATLAB crashed during plotting.
% load(savefilenm);
%% plotting
maparea.lon=[-98,-79]; % temporary 
maparea.lat=[34.4, 48];
gmtmaparea=[num2str(maparea.lon(1)) '/' num2str(maparea.lon(2)) '/' ...
              num2str(maparea.lat(1)) '/' num2str(maparea.lat(2))];
gmtmapsize=[num2str(mean(maparea.lon)) '/' num2str(mean(maparea.lat)) '/' ...
    num2str(maparea.lat(1)) '/' num2str(maparea.lat(2)) '/2.5i'];
%
%caxisdataabs=[3.8 4.25; 3.8 4.2;3.8 4.35;3.85 4.15;4.3 4.75; 4.3 4.75]; %1/3 moho version
if strcmp(velocitytag,'S')
    %depth [15,22,31,40,51,63,74,97];
%     caxisdataabs=[4.3 4.8;];
    caxisdataabs=[...
        3.7 4.1; 3.7 4;3.7 4;4.0 4.4;4.45 4.85;...
        4.45 4.85;4.45 4.85; 4.45 4.85;4.45 4.85;...
        4.45 4.85;4.45 4.85;4.45 4.85;];
    ctickint=0.1;
elseif strcmp(velocitytag,'P')
%     caxisdataabs=[6 7.5;6 7.5;6 7.5;6 7.5;7.8 8.4;7.8 8.4];
    caxisdataabs=[5.7 6.5;6 7;6 7;6 7;7.8 8.4;7.8 8.4];
    ctickint=0.2;
elseif strcmp(velocitytag,'PR')
    caxisdataabs=[.1 .275;.1 .275;.1 .275;.1 .3;.225 .325;.225 .325];
    ctickint=0.025;
elseif strcmp(velocitytag,'RHO')
    caxisdataabs=[2.75 3;2.75 3;2.75 3;2.75 3;3.1 3.2;3.1 3.3];
    ctickint=0.05;
end
%caxisdataabs=[3.8 4.1; 3.6 4.3;3.8 4.5;4.4 4.9;4.3 4.8; 4.3 4.9];
caxisdatarel=[-8 8;-8 8;-8 8;-8 8;-8 8; -8 -8];

figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) ','(i) ','(j) '};

proflabel={'A','A'' ';'B','B'' ';'C','C'' ';'D','D'' ';'E','E'' ';'F','F'' '};

YS=squeeze(YSIM(npml:end-npml,npml:end-npml,1));
XS=squeeze(XSIM(npml:end-npml,npml:end-npml,1));

plotbedrock=0;
if plotbedrock
    bedrock=load('NewEnglandBedRock');
    bedrocklabelsize=14;
end
plotstations=0;
if plotstations
    stations=load('../Stations/craton_station_withdata.txt');
end
%%% load structural outlines
plotoutline=1;
if plotoutline
    load('outlines_NA-MC.mat');
    noutline=length(lineinfo);
end
%velocity regions;
% load extracted_velocity_outline.mat
project_quakes=0;
% earthquakefile='AK_quakes_1970_gt3.5_NEIC_matlab.txt';
earthquakefile='AK_quakes_1990to2020_gt3_NEIC_selectedcolumns.txt';
% earthquakefile='southernAK_dd_ge3.xyd';
if project_quakes
    quakes=load(earthquakefile);
    quakeedgecolor=[.4 .4 .4];
    quakelinewidth=1;
    quakesizescale=0.15;
    quakeserror=2;
end
plotgmtfig=0;
figure('Position',[400 400 1750 930]);
for iz=1:size(vout,1) 
    %
    clear amask raycover;
    load(raycoverlist{iz});
    amask=zeros(length(squeeze(XSIM(:,1,1))),length(squeeze(YSIM(1,:,1))));
    for i=1:size(amask,1)
            clear id00;
            id00=inpolygon(YSIM(1,:,1)-360,XSIM(i,1,1)*ones(size(amask,2),1),raycover.data(:,1),...
                    raycover.data(:,2));
            amask(i,id00)=1;
    end
    amaskcut=amask(npml:end-npml,npml:end-npml);
    disp(['Plotting ',vout{iz}.tag,' ...']);
    clear vplot;
    clear amaskcuttemp;
    amaskcuttemp=amaskcut;
    vplot=vout{iz}.data(npml:end-npml,npml:end-npml);
    vplot(amaskcuttemp==0)=nan;
    amaskcuttemp(isnan(vplot))=0;
    
    if project_quakes
        clear qidx;
        qidx=find(quakes(:,3) >= mean(vout{iz}.depth)-quakeserror ...
            & quakes(:,3) <=mean(vout{iz}.depth)+quakeserror);
        if plotgmtfig
            quakestmpfile='quakestmp.txt';
            fidqks=fopen(quakestmpfile,'w');
            for qq=1:length(qidx)
               fprintf(fidqks,'%g  %g  %g %g\n',quakes(qidx(qq),1),quakes(qidx(qq),2),...
                   quakes(qidx(qq),3),quakes(qidx(qq),4));
            end
            fclose(fidqks);
        end
    end
    vmean=nanmean(nanmean(vplot));
    subplot(3,4,iz); 
    
    hold on; box on; axis on;
    switch iplot
        case 0
          cabsmin=caxisdataabs(iz,1)-0.1;cabsmax=caxisdataabs(iz,2)+0.1;
          
          if plotgmtfig
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % plot using GMT matlab API to single PS file for each layer.
              % To create a Grid structure from a 2-D Z array and a 1x9 header vector:
              % G = gmt ('wrapgrid', Z, head)
              % header is a vector with [x_min x_max, y_min y_max z_min z_max reg x_inc y_inc]
              % reg: registration. we use 0.
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              clear vG vcpt gheader psfilenm;
              gheader=[min(YS(1,:)-360) max(YS(1,:)-360) min(XS(:,1)) max(XS(:,1)) ...
                  min(min(vplot)) max(max(vplot)) 0 abs(mean(diff(YS(1,:)))) abs(mean(diff(XS(:,1))))];
              % create cpt file
              vcpt=['myjetwr_' num2str(cabsmin) '_' num2str(cabsmax) '.cpt'];
              gmt('gmtset MAP_FRAME_TYPE = plain');
              gmt('gmtset MAP_FRAME_WIDTH 0.1c');
              gmt('gmtset FONT_TITLE = 14p,Helvetica,black');
              gmt('gmtset MAP_TITLE_OFFSET = 0p');
              gmt(['makecpt -Cjetwr -T' num2str(cabsmin) '/' num2str(cabsmax) ...
                  '/.01 -D > ' vcpt]);
              %gmt('destroy');
              psfilenm=[ite_nm,'_', velocitytag,'_', num2str(vout{iz}.depth(1)),...
                  '_to_',num2str(vout{iz}.depth(1)),'.ps'];
              vG=gmt('wrapgrid',flip(vplot),gheader);
              %gmt('destroy');
              gmt(['grdimage -R' gmtmaparea ' -JL' gmtmapsize ' -C' vcpt ' -Bf3a6/f1a2:."' ...
                  vout{iz}.tag '"'...
                  ':SEWn -X4 -Y13 -P -K > ' psfilenm],vG);
              %gmt('destroy');
              gmt(['pscoast -R -JL -Di -N2/.5p,150 -N1/2p,150 -W.5p,150 -A50+r -O -K >> ' psfilenm]);
              %gmt('destroy');
              %plot profile locations.

%               gmt(['psxy -R -JL DenaliFault.txt ',...
%                   '  -W1p,255 -O -K >> ', psfilenm]);
%               gmt(['grdcontour -JL -R ',slabgridfile, ' -C+',num2str(mean(vout{iz}.depth)),...
%                   ' -W.5p,0,solid -K  -O >> ', psfilenm]);
              %gmt('destroy');
              % dark purpler:-G200/0/200
              gmt(['psxy -JL -R -W1p,0,dashed -K -O YakutatOutline.txt >> ', psfilenm]);

%               gmt(['psxy -R -JL AKvolclatlong_ready_matlab.txt -St0.2 ',...
%                   '  -W.5p,255 -Gred -O -K >> ', psfilenm]);
              if project_quakes
                 gmt(['psxy -R -JL ',quakestmpfile,' -Sc0.08 ',...
                    '  -W.2p,255 -G50 -O -K >> ', psfilenm]); 
              end
              if iz==3 || iz==4 || iz==5
              gmt(['psxy -R -JL vertical_profile_locations4gmt.txt -O -K -A -W2p,0 >> ' psfilenm]);
              %plot profile labels.
              gmt(['pstext -R -JL vertical_profile_labels4gmt_new.txt -F+f+a+j -O -K >> ' psfilenm]);
              end
              system(['echo -145 58 12 0 CB mean: ',num2str(vmean,2),'>stmp']);
              gmt(['pstext -R -JL stmp -F+f+a+j -O -K >> ' psfilenm]);
              %gmt('destroy');
              gmt(['psscale -D3.2/-1/4.7/0.25h -C' vcpt ' -B0.2 -O >> ' psfilenm]);
              gmt('destroy');
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %%%%%%%%%%%%% end of plotting by GMT API %%%%%%%%%%%%%%%%%%%%%%
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          end
          h=image(YS(1,:)-360,XS(:,1),vplot);
%           [CC,Ch]=contour(YS(1,:)-360,XS(:,1),vplot,4.3:0.2:4.7,'k',...
%               'Color',[0.6 0.6 0.6],'LineWidth',1);
          %caxis([cabsmin cabsmax]);
          amaskcuttemp(vplot >4.4 & vplot<4.7)=0;
          h.CDataMapping='scaled';
          h.AlphaData=amaskcuttemp;
          colormap('jetwr'); 
            hcbar=colorbar('eastoutside');
          set(gca,'CLim',[cabsmin cabsmax]);
          set(hcbar,'TickDirection','out','Ticks',cabsmin:ctickint:cabsmax);
          
          if strcmp(velocitytag,'P')
               hcbar.Label.String='V_P (km/s)';
          elseif strcmp(velocitytag,'S')
               hcbar.Label.String='V_S (km/s)';
          elseif strcmp(velocitytag,'PR')
              hcbar.Label.String='Poisson''s Ratio';
          elseif strcmp(velocitytag,'RHO')
              hcbar.Label.String='Density (kg/m^3)';
          end
%           hcbar.Label.String='V_S (km/s)';
          hcbar.FontSize=12;
        case 1
          crelmin=caxisdatarel(iz,1);crelmax=caxisdatarel(iz,2);
          clear vpercent;
          
          vpercent = 100*(vplot-vmean)/vmean;
          h=image(YS(1,:)-360,XS(:,1),vpercent);
          h.CDataMapping='scaled';
          h.AlphaData=amaskcuttemp;
          %caxis([crelmin crelmax]);
          set(gca,'CLim',[crelmin crelmax]);
          colormap('jetwr'); 
            hcbar=colorbar('eastoutside');
          set(hcbar,'TickDirection','out','Ticks',crelmin:2:crelmax);
          if strcmp(velocitytag,'P')
               hcbar.Label.String='dV_P (km/s)';
          elseif strcmp(velocitytag,'S')
               hcbar.Label.String='dV_S (km/s)';
          elseif strcmp(velocitytag,'PR')
              hcbar.Label.String='\delta\sigma';
          end
          hcbar.FontSize=12;
    end
    
    
    %title([figlabel{iz},' ',vout{iz}.tag],'fontsize',14);
    title([vout{iz}.tag],'fontsize',14);
    axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
    daspect([1 cosd(mean(maparea.lat)) 1]);

    for sb=1:length(state)
        plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[.5 .5 .5],'LineWidth',.75);
    end
 
    %
    if plotoutline
        for no = 1:noutline
            tag=lineinfo{no}.tag;
            disp(tag);
            if ~strcmp(tag,'NA-MC')
                if contains(tag,'Basin') || contains(tag,'Rift') || contains(tag,'MCR')
                    plot(lineinfo{no}.data(:,1),lineinfo{no}.data(:,2),'m-','linewidth',1);
                else
                    plot(lineinfo{no}.data(:,1),lineinfo{no}.data(:,2),'m--','linewidth',1);
                end
            end
        end
    end
    %plot bedrockoutline
    if plotbedrock %&& (iz<=4)
        for k=1:length(bedrock.boundaries) %1:4,7 is specially selected to include the Adirondacks
            if strcmp(bedrock.boundaries{k}.tag,'Adirondacks')
                plot(bedrock.boundaries{k}.data(:,1), bedrock.boundaries{k}.data(:,2),...
                    'k-','LineWidth',1.5); 
            elseif strcmp(bedrock.boundaries{k}.tag,'Laurentia-TB')
                 plot(bedrock.boundaries{k}.data(:,1), bedrock.boundaries{k}.data(:,2),...
                    'w-','LineWidth',2);
            end
        end
    end
    if plotstations && (iz==1)
        plot(stations(:,1)-360,stations(:,2),'k^','markersize',3,'linewidth',.5,'markeredgecolor',[.2 .2 .2]);
        plot(288.4727-360, 45.282700,'k^','markersize',7,'linewidth',.5,'markerfacecolor','k')
        plot(-74.2228, 43.9734,'ks','markersize',10,'linewidth',.5,'markerfacecolor','k','markeredgecolor','w')
    end
    if project_quakes
       plot(quakes(qidx,1),quakes(qidx,2),'k.');
    end
    %text of smooth scale in km
%     text(mean(maparea.lon)-3, maparea.lat(1)+1, ['smooth=' num2str(smoothscale(iz)) ' km'],'FontSize',14);
    set(gca,'TickDir','out');
    
%     if iz==length(zgridtoplist)
%         for np=1:4%plongridnum,
%             slat = ptLat1(np); 
%             slon = ptLon1(np)-360; 
%             elat = ptLat2(np); 
%             elon = ptLon2(np)-360;
% 
%             plot([slon elon],[slat elat],'k-','linewidth',2,'color',[0 .7 0]);
%             text(slon,slat+0.1,['\bf ' proflabel{np,1}],'fontsize',14,'color','k');
%             text(elon,elat+0.1,['\bf ' proflabel{np,2}],'fontsize',14,'color','k');
%         end
%     end
    hold off;
    set(gca,'FontSize',14);
    %box on;
    drawnow;
end