%plot isosurface
close all;clear all;
% MFILEROOT='../mfiles';
% path([MFILEROOT '/fun-spool'],path);

%% load master data
load ../00masterdatafiles/us_states.mat;
canada_provs=shaperead('PROVINCE.SHP','UseGeoCoords',true);
% load('outlines_NA-MC.mat');
% noutline=length(lineinfo);

%%
ite_nm = 'ite_0.05deg_07';
velocitytag='S'; % 'P' for P velocities, 'PR' for Poisson's ratios.
% idx = 0; % 0, absolute velocity; 1, velocity perturbation
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
minlat=XSIM(end,1,end);maxlat=XSIM(1,1,end);
minlon=YSIM(1,1,end);maxlon=YSIM(1,end,end);

mrh=gather_media(snapinfo,'rho','mediadir',dir_media);
mmu=gather_media(snapinfo,'mu','mediadir',dir_media);
mla=gather_media(snapinfo,'lambda','mediadir',dir_media);
mvp=((mla+2*mmu)./mrh).^0.5;
mvs=(mmu./mrh).^0.5;

x=YSIM(1,:,1)-360;
y=XSIM(:,1,1);
z=6371-squeeze(ZSIM(1,1,:))/1000;
%%
master_pars;
maparea.lat(1) = 37; 
% maparea.lon(1) = -96;
%
lonidx=find(x>=maparea.lon(1) & x<=maparea.lon(2));
latidx=find(y>=maparea.lat(1) & y<=maparea.lat(2));
%% MLD grid

dinfile='Liu2018_US_NVD_depth.txt';
clear datat;
datat=load(dinfile);
datax=datat(:,2);
datay=datat(:,1);
grid_dx=0.5;
grid_dy=grid_dx;
gridx=maparea.lon(1):grid_dx:maparea.lon(2);
gridy=maparea.lat(1):grid_dy:maparea.lat(2);

gridd=griddata(datax,datay,datat(:,3),gridx,gridy');

%%
gridplot=smooth2a(gridd,1);
% h=pcolor(gridx,gridy,gridplot);
% % h.FaceAlpha=0.75;
% h.CDataMapping='scaled';
% colormap(flip(jet(20)));  %(length(cabsmin:ctickint/10:cabsmax)-1)
% hcbar=colorbar('eastoutside');
% set(gca,'CLim',[70,95]);
% set(hcbar,'TickDirection','out');
% shading interp;
% hcbar.Label.String='Depth (km)';
%            hcbar.Label.String='V_S (km/s)';
% hcbar.FontSize=12;
% axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
% daspect([1 cosd(mean(maparea.lat)) 1]);
%% low-velocity volumes in the mantle lithosphere
myzlimUM=[60 200];
smoothscale=[17 17 3];
vplot_all=smooth3(mvs/1000,'box',smoothscale); %
amask=nan(length(squeeze(XSIM(:,1,1))),length(squeeze(YSIM(1,:,1))));
load('RayCoverOutline_ite_0.05deg_07_40-75s_cutoff10.mat');
for i=1:size(amask,1)
    clear id00;
    id00=inpolygon(YSIM(1,:,1)-360,XSIM(i,1,1)*ones(size(amask,2),1),raycover.data(:,1),...
            raycover.data(:,2));
    amask(i,id00)=1;
end
amask3d=nan(size(vplot_all));
for k=1:length(z)
    amask3d(:,:,k)=amask;
end
vplot=vplot_all.*amask3d;

depthidx2=find(z>=myzlimUM(1) & z<=myzlimUM(2));
%
clear isov isocap_low;
zlimindexmin=min(depthidx2); zlimindexmax=max(depthidx2);
isovalue1=4.4;
zlimindexmax_slow=zlimindexmax;
isovUM_low=isosurface(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax_slow),...
    vplot(latidx,lonidx,zlimindexmin:zlimindexmax_slow),isovalue1);
isocap_low=isocaps(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax_slow),...
    vplot(latidx,lonidx,zlimindexmin:zlimindexmax_slow),isovalue1,'below');

isovalue2=4.85;
zlimindexmax_slow=zlimindexmax;
isovUM_high=isosurface(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax_slow),...
    vplot(latidx,lonidx,zlimindexmin:zlimindexmax_slow),isovalue2);
isocap_high=isocaps(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax_slow),...
    vplot(latidx,lonidx,zlimindexmin:zlimindexmax_slow),isovalue2,'above');

%% plot 3d isosurface.
bgcolor='w'; %'w' for print; 'k' for screen play
figure('position',[1400 400 1280 720],'color',bgcolor);
surfacealpha=1;
view(3); axis tight

hold on;

pUM_low=patch(isovUM_low,'FaceColor',[1 .8 0],...
   'EdgeColor','none','FaceAlpha',surfacealpha);
% alpha(0.7);
isonormals(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax),...
    vplot(latidx,lonidx,zlimindexmin:zlimindexmax),pUM_low);
patch(isocap_low,'FaceColor','interp','EdgeColor','none','FaceAlpha',surfacealpha)
% patch(isocap_low,'FaceColor',[0.8,.2,0],'EdgeColor','none','FaceAlpha',surfacealpha)
% 
pUM_high=patch(isovUM_high,'FaceColor',[0,0.2,0.9],...
   'EdgeColor','none','FaceAlpha',surfacealpha);
% alpha(0.7);
isonormals(x(lonidx),y(latidx),z(zlimindexmin:zlimindexmax),...
    vplot(latidx,lonidx,zlimindexmin:zlimindexmax),pUM_high);
patch(isocap_high,'FaceColor','interp','EdgeColor','none','FaceAlpha',surfacealpha)

% plot MLD
% gridplot=smooth2a(gridd,1);
% hs=surface(gridx,gridy,gridplot,'FaceColor','interp','EdgeColor','none');
% hs.FaceAlpha=0.75;
% hs.CDataMapping='scaled';
% colormap(flip(summer(10)));  %(length(cabsmin:ctickint/10:cabsmax)-1)
% hcbars=colorbar('eastoutside');
% % set(gca,'CLim',[70,95]);
% set(hcbars,'TickDirection','out');
% % shading interp;
% hcbars.Label.String='Depth (km)';
% hcbars.FontSize=12;
% caxis([70,95])

colormap(jetwr);
caxis([4. 5]);
view(25,51); %axis tight
camlight right
% camlight left
% lightangle(125,65);
% lightangle(90,65); 
lightangle(-50,-60);
light('Position',[0.1 0 -1])
lighting gouraud
set(gca,'ZDir','reverse');
hold on;

flist=dir('outline_*.txt');

if strcmp(bgcolor,'k')
    for sb=1:length(state)
        plot3(state(sb).polygon(:,1), state(sb).polygon(:,2),myzlimUM(1)*ones(length(state(sb).polygon(:,1)),1),...
            'color',[.9 .9 .9],'LineWidth',1);
    end
    outlinecolor='m';
    for k =1:length(flist)
        infile=flist(k).name;
        outline=load(infile);
        if contains(infile,'basin') || contains(infile,'rift') || ...
                contains(infile,'mcr') || contains(infile,'rr') || contains(infile,'dome')
            %close polygon
            outline(end+1,:)=outline(1,:);
        end
        if contains(infile,'basin')
            plot3(outline(:,1),outline(:,2),myzlimUM(1)*ones(length(outline(:,2)),1),...
                'g--','linewidth',1.5,'color',outlinecolor);
        elseif contains(infile,'rift') || contains(infile,'mcr') || contains(infile,'rr') 
            plot3(outline(:,1),outline(:,2),myzlimUM(1)*ones(length(outline(:,2)),1),...
                'g:','linewidth',1.5,'color',outlinecolor);
        elseif contains(infile,'dome')
            plot3(outline(:,1),outline(:,2),myzlimUM(1)*ones(length(outline(:,2)),1),...
                'g:','linewidth',2,'color',outlinecolor);
    %             else
    %                 plot(outline(:,1),outline(:,2),'k-','linewidth',1,'color',0.5*[1,1,1]);
        end
    end
else
    for sb=1:length(state)
        plot3(state(sb).polygon(:,1), state(sb).polygon(:,2),myzlimUM(1)*ones(length(state(sb).polygon(:,1)),1),...
            'color',[.4 .4 .4],'LineWidth',.5);
    end
    for sb=1:length(canada_provs)
        plot3(canada_provs(sb).Lon, canada_provs(sb).Lat,...
            myzlimUM(1)*ones(length(canada_provs(sb).Lon),1),...
            'color',[.4 .4 .4],'LineWidth',.5);
    end
    outlinecolor='m';
    for k =1:length(flist)
        infile=flist(k).name;
        outline=load(infile);
        if contains(infile,'basin') || contains(infile,'rift') || ...
                contains(infile,'mcr') || contains(infile,'rr') || contains(infile,'dome')
            %close polygon
            outline(end+1,:)=outline(1,:);
        end
        
        if contains(infile,'basin')
            fill3(outline(:,1),outline(:,2),myzlimUM(1)*ones(length(outline(:,2)),1),'w','facealpha',0.,...
            'edgealpha',1,'edgecolor','m','LineStyle','--','linewidth',1.5)
%             plot3(outline(:,1),outline(:,2),myzlimUM(1)*ones(length(outline(:,2)),1),...
%                 'g--','linewidth',1.5,'color',outlinecolor);
        elseif contains(infile,'rift') || contains(infile,'mcr') || contains(infile,'rr') 
            plot3(outline(:,1),outline(:,2),myzlimUM(1)*ones(length(outline(:,2)),1),...
                'g:','linewidth',1.5,'color',outlinecolor);
        elseif contains(infile,'dome')
            plot3(outline(:,1),outline(:,2),myzlimUM(1)*ones(length(outline(:,2)),1),...
                'g:','linewidth',2,'color',outlinecolor);
        elseif contains(infile,'mgborder')
            plot3(outline(:,1),outline(:,2),myzlimUM(1)*ones(length(outline(:,2)),1),...
                'g-','linewidth',3,'color',0.3*[1 1 1]);
    %             else
    %                 plot(outline(:,1),outline(:,2),'k-','linewidth',1,'color',0.5*[1,1,1]);
        end
    end

    hold off;
end

axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2) myzlimUM(1)-1 myzlimUM(2)]);

hold off

box on;
ax = gca;
ax.BoxStyle = 'full';
% grid on;
daspect([1 cosd(mean(maparea.lat)) 30]);
%
% hcbar=colorbar('eastoutside');
% set(hcbar,'TickDirection','out','Ticks',4:.2:6);
% 
% hcbar.Label.String='V_S (km/s)';

if strcmp(bgcolor,'k')
    xlabel('Longitude','fontsize',16,'color','w');
    ylabel('Latitude','fontsize',16,'color','w');
    zlabel('Depth (km)','fontsize',16,'color','w');
    set(gca,'Fontsize',16,'ZTick',myzlimUM(1):40:myzlimUM(2),'ZTickLabel',myzlimUM(1):40:myzlimUM(2),...
        'YTick',maparea.lat(1):2:maparea.lat(2),'XTick',round(maparea.lon(1):3:round(maparea.lon(2)))); 
    ax.XColor='w';
    ax.YColor='w';
    ax.ZColor='w';
    hcbar.Color='w';
    hcbar.Label.Color='w';
    set(gca,'Color','k')
else
    xlabel('Longitude','fontsize',16);
    ylabel('Latitude','fontsize',16);
    zlabel('Depth (km)','fontsize',16);
    set(gca,'Fontsize',16,'ZTick',myzlimUM(1):40:myzlimUM(2),'ZTickLabel',myzlimUM(1):40:myzlimUM(2),...
        'YTick',maparea.lat(1):2:maparea.lat(2),'XTick',round(maparea.lon(1):3:round(maparea.lon(2))));
%zlim([40 100]);
end