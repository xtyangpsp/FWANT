%
close all;
% /opt/matlab/2007b/bin/matlab
% MFILE_ROOT='../mfiles';
% path([MFILE_ROOT '/fun-spool'],path);
% path([MFILE_ROOT '/saclab'],path);

ite_nm = ['ite_0.05deg_09'];
fnm_conf=['./SeisFD3D.conf_' ite_nm];

iplot = 0; % 0, absolute velocity; 1, velocity perturbation

%read updated model (same conf and coord)
dir_media=['./updated_input_' ite_nm];
dir_coord=['./updated_input_' ite_nm];

disp(['Read model... ',dir_media]);

id = 0; subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
               indxem=[];
               indxkp=[];
n=1;
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
[XSIM,YSIM,ZSIM]=gather_coord(snapinfo,'coorddir',dir_coord);
% convert from radian to degrees
XSIM=90-XSIM*180/pi; %latitude
YSIM=YSIM*180/pi;
%define the area of plot (exclude pmls)
npml=12; %number of pml layers
minlat=XSIM(end-npml,1,end);maxlat=XSIM(1+npml,1,end);
minlon=YSIM(1,1+npml,end);maxlon=YSIM(1,end-npml,end);

mrh=gather_media(snapinfo,'rho','mediadir',dir_media);
mmu=gather_media(snapinfo,'mu','mediadir',dir_media);
mla=gather_media(snapinfo,'lambda','mediadir',dir_media);
mvp=((mla+2*mmu)./mrh).^0.5;
mvs=(mmu./mrh).^0.5;
%mvs=smooth3(mvs,'box',[3 3 1]);
%%
topofile='ETOPO1_Bed_g_gmt4.grd';
%
load ../00masterdatafiles/us_states.mat;

%plot topography and locations of the cross-sections
tlon=nc_varget(topofile,'x');tlat=nc_varget(topofile,'y');topogrid=nc_varget(topofile,'z');


%%
velocitytag='S';
mvbkp=nan(size(mvs));
if strcmp(velocitytag,'P')
    mvbkp=smooth3(mvp,'box',[3 3 1])/1000;
    %dmv=dmvp;
    figtitletag='Vp';
elseif strcmp(velocitytag,'S')
    mvbkp=smooth3(mvs,'box',[7 7 3])/1000;
    figtitletag='Vs';
end
%
mv=mvbkp;
if iplot==1
   clear vmean3d;
   vmean3d=squeeze(mean(mean(mvbkp)));
   for i=1:length(vmean3d)
       mv(:,:,i)=100*(mvbkp(:,:,i)-vmean3d(i))/vmean3d(i);
   end
end
%% prepare grid file for projecting lines on the profiles.
clear xyzgrid xyz_x xyz_y;
xyzfilelist4mask={'US_Crust_GRL_2015_CrustThickness.txt',...
    'CUS_Moho_latlon.txt','Liu2018_US_NVD_depth.txt'};
xyzfilelistsmooth={'US_Crust_GRL_2015_CrustThickness.txt',...
    'CUS_Moho_latlon.txt','Liu2018_US_NVD_depth.txt'};
grididxlist=[1];
grid_dx=0.5;grid_dy=grid_dx;
for k =1:length(grididxlist)
    grididx=grididxlist(k);
    %
    clear datat;
    datat=load(xyzfilelistsmooth{k});
    datax=datat(:,2);
    datay=datat(:,1);
    xyz_x{grididx}=min(datax):grid_dx:max(datax);
    xyz_y{grididx}=min(datay):grid_dy:max(datay);
    xyzgrid{grididx}=griddata(datax,datay,datat(:,3),xyz_x{grididx},xyz_y{grididx}');

    clear moho;
    moho=load(xyzfilelist4mask{k});
    mohoidx=boundary(moho(:,2),moho(:,1),1);
    mohooutline.lon=moho(mohoidx,2);
    mohooutline.lat=moho(mohoidx,1);
    %mask moho
    mohomask=nan(size(xyzgrid{grididx}));
    for i=1:size(mohomask,1)
            clear id00;
            id00=inpolygon(xyz_x{grididx},xyz_y{grididx}(i)*ones(size(mohomask,2),1),...
                mohooutline.lon,mohooutline.lat);
            if ~isempty(id00);mohomask(i,id00)=1;end
    end
    xyzgrid{grididx}=mohomask.*xyzgrid{grididx};
end

%
% save('surfacegrids_combined_E125_Zhang2019GRL_Moho_heatflow.mat','xyz_x','xyz_y','xyzgrid');
%%
master_pars;
% ptLon1 = [-91,   -85, -93,   -91.4,   -95]; ptLon1=ptLon1+360;
% ptLon2 = [-80,   -84.,   -84.5,   -85.7,    -83];ptLon2=ptLon2+360;
% ptLat1 = [43,    47,     37, 43.5,    39.5];
% ptLat2 = [43.8,   39.5,     42,  35.5,    39.5];
%for paper
ptLon1 = [-89.5,   -87.5,   -93,   -91.4]; ptLon1=ptLon1+360;
ptLon2 = [-79.5,   -82,   -84.,  -85.7];ptLon2=ptLon2+360;
ptLat1 = [41.3,  46.8,    37.5,  43.5];
ptLat2 = [45.4,   39.8,  41.5,  35.5];

%subplotpar=[3,1];
plongridnum = length(ptLon1);
platgridnum = length(ptLat1);

proflabel={'A','A'' ';'B','B'' ';'C','C'' ';'D','D'' ';'E','E'' ';'F','F'' ';'G','G'' ';...
    'H','H'' ';'I','I'' ';'J','J'' ';'K','K'' ';'L','L'' ';'M','M'' ';'N','N'' ';'O','O'' '};

%%% load structural outlines
% load('outlines_NA-MC.mat');
% noutline=length(lineinfo);

figure('Position',[100,800,500,450]); hold on;
h0=image(tlon,tlat,topogrid); hold on;
h0.CDataMapping='scaled';
demcmap([-200,1000]);
hcbar=colorbar('eastoutside');
set(hcbar,'YDir','normal','TickDirection','out');
hcbar.Label.String='Elevation (m)';
hcbar.FontSize=14;
set(gca,'TickDir','out','YDir','normal');
hl=[];

axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
daspect([1 cosd(mean(maparea.lat)) 1]);
    
for sb=1:length(state)
    plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',0.9+[.0 .0 .0],'LineWidth',1);
end
flist=dir('outline_*.txt');
outlinecolor=0.*[1,1,1];
for k =1:length(flist)
    infile=flist(k).name;
    outline=load(infile);
    if contains(infile,'basin') || contains(infile,'rift') || ...
            contains(infile,'mcr') || contains(infile,'rr') || contains(infile,'dome')
        %close polygon
        outline(end+1,:)=outline(1,:);
    end
    if contains(infile,'basin')
        plot(outline(:,1),outline(:,2),'g--','linewidth',1.5,'color',outlinecolor);
    elseif contains(infile,'rift') || contains(infile,'mcr') || contains(infile,'rr') 
        plot(outline(:,1),outline(:,2),'g-','linewidth',1.5,'color',outlinecolor);
    elseif contains(infile,'dome')
        plot(outline(:,1),outline(:,2),'g:','linewidth',2,'color',outlinecolor);
%             else
%                 plot(outline(:,1),outline(:,2),'k-','linewidth',1,'color',0.5*[1,1,1]);
    end
end
% save profile locations and labels for GMT plotting
fout_xsection_locations='vertical_profile_locations4gmt.txt';
fout_xsection_labels='vertical_profile_labels4gmt.txt';

fido_locations=fopen(fout_xsection_locations,'w');
fido_labels=fopen(fout_xsection_labels,'w');
for np=1:plongridnum
% start point (slat slon)  and end point (elat elon) for plotting vertical image profile
    slat = ptLat1(np); 
	slon = ptLon1(np); 
    elat = ptLat2(np); 
	elon = ptLon2(np);
				    
	proflength = deg2km(distance(slat, slon, elat, elon));
				    
	% number of pts on this profile
	profptnum = round((proflength/20));
	% final pts (flat flon) that are shown on the profile
    % 	flat = slat + (0:(profptnum-1))*(elat - slat)/(profptnum - 1);
    % 	flon = slon + (0:(profptnum-1))*(elon - slon)/(profptnum - 1);
    [flat,flon]=track2(slat,slon-360,elat,elon-360,wgs84Ellipsoid,'degrees',profptnum); %get great cirlce path
    flon=360+flon;
    % save points along the profile.
    fprintf(fido_locations,'> %s\n',char(proflabel{np}));
    for kk=1:length(flat)
        fprintf(fido_locations,'%g   %g\n',flon(kk)-360,flat(kk));
    end

    fprintf(fido_labels,'%g   %g   %d,Helvetica,black   %d   %s   %s\n',...
        slon-.5,slat-.5,10,0,'CB',proflabel{np,1});
    fprintf(fido_labels,'%g   %g   %d,Helvetica,black   %d   %s   %s\n',...
        elon+.5,elat,10,0,'CB',proflabel{np,2});
    % continue to plot the profiles.
    plot(flon-360,flat,'-','linewidth',3,'color','m');
    text(slon-360-0.8,slat-.3,proflabel{np,1},'fontsize',16,'color','m');
    text(elon-360+0.2,elat-.3,proflabel{np,2},'fontsize',16,'color','m');
%     text(elon-360,elat+.1,proflabel{np,2},'fontsize',16,'color','m');
end
fclose(fido_locations);
fclose(fido_labels);  
box on;
axis on;
% legend(hl,'Rifts/Basins','Domes/Arches','fontsize',12)
set(gca,'FontSize',14,'TickDir','out');
hold off;
set(gcf,'renderer','Painters','PaperPositionMode','auto'); 
eval(['print -dpdf -r300 -painters topography_profilelocations.pdf']);
% legend([hl1 hl2 hl3],'SL2014','Menke2016','This study');
% saveas(gca,'mohomapwithlocations.eps','epsc');

%% get topo data
vert_exaggeration=5;
mindepth=6;
mindepthaxis=-10;
depthlabel_interval=10;
maxdepth=60; %20 for P, 30 for S
depthinterval=1;
figheightscale=2.8; %2 for 120 depth
figwidthscale=1.7; %1.5 for 120 depth
xlabelunit='deg'; %deg (for degrees) or km (for km distance)
xlabeldegtype='auto'; %auto (automatically decides which of lat or lon is used;
                      %lon: use lontitude
                      %lat: use latitude.
nmarkpt=40; %degree to km coefficient.
plotcontour=0;
if plotcontour
    contouritv=[4.2,4.2];%4.4:0.4:4.8;
    clinecolor=0.6*[1 1 1];
    clinewidth=0.5;
    clabelcolor=0.6*[1 1 1];
    clabelspace=600;
    clabelsize=10;
end

imagetype='image'; %could be: contourf (FILLED CONTOUR) OR image
imagecontourlevels=3:0.01:5; %only effective when imagetype is contourf.
if iplot==0
    colorinterval=0.2;
    if strcmp(velocitytag,'P')
        colorrange=[6 7];
    elseif strcmp(velocitytag,'S')
        colorrange=[3.6,4.6];%[3.7 4.2] for 5-40;
    end
elseif iplot==1
    colorinterval=2;
    colorrange=[-10 10];
end
projectsurface=1;
projectsurface_line={'m--','m-','k:'};
projectsurface_linewidth=[1,1,1.5];

plottopo=1;
if plottopo
    topocolor={[0 0.8 .8],[.3 .3 .3]};
    toposhift=5; % scale and shift the topo data by (mindepth - toposhift)
    topomax=4; % scale facter for the maximum topo values. 
end
project_quakes=0;
if project_quakes
    quakedata=readtable("NACraton_quake_200001_202210.csv");
    quakedata_xyz=[quakedata{:,"longitude"},quakedata{:,"latitude"},quakedata{:,"depth"},...
        quakedata{:,"mag"}];
    quakeedgecolor=[.8 .8 .8];
    quakelinewidth=1;
    quakesizescale=2;
    quakes_proximity=5;
end

project_stations=0;
if project_stations
    stationfile='../00masterdatafiles/NACraton_station.txt';
    [stationlon,stationlat,stationele]=textread(stationfile,'%*s %*s %f %f %f\n');
    stationlon=stationlon-360;
    station_xyz=[stationlon,stationlat,stationele];
    station_proximity=2.5;
    stationsymblecolor='^r';
    stationsize=40;
end
for np=1:plongridnum
    disp(['Working on [ ',num2str(np),' / ',num2str(plongridnum),' ]']);
% start point (slat slon)  and end point (elat elon) for plotting vertical image profile
    slat = ptLat1(np); 
	slon = ptLon1(np); 
    elat = ptLat2(np); 
	elon = ptLon2(np);
				    
	[proflength,az]=distance(slat, slon, elat, elon,wgs84Ellipsoid);
    proflength=proflength/1000; %from meters to km.
    latdist=distance(slat, slon, elat, elon,wgs84Ellipsoid)/1000;
    %nmarkpt = round(proflength/deg2kmcoe);
	markptdeltadist = proflength/(nmarkpt-1);
	%markptdist = zeros(1, nmarkpt);
    markptdist = 0:markptdeltadist:proflength;
% 	markptlon = zeros(1, nmarkpt);
% 	markptlat = zeros(1, nmarkpt);
												    
	% number of pts on this profile
	profptnum = 40*round((proflength/5));
	% final pts (flat flon) that are shown on the profile
    [flat,flon]=track2(slat,slon-360,elat,elon-360,wgs84Ellipsoid,'degrees',profptnum);
    flon=360+flon;
    
	dist = zeros(1, profptnum);
    for i = 1:(profptnum-1)
        ptdist = deg2km(distance(flat(i), flon(i), flat(i+1), flon(i+1)));
        dist(i+1) = ptdist + dist(i);
    end
	markptlon = interp1(dist, flon, markptdist, 'pchip');
	markptlat = interp1(dist, flat, markptdist, 'pchip');
																												       
	intpdist = 0:(proflength/(profptnum-1)):proflength;
	proflength = dist(profptnum);
	depth=6371-abs(ZSIM(npml,npml,:)/1000);
	[distmesh,depthmesh] = meshgrid(dist, depth);
    idepth=mindepth:depthinterval:maxdepth;
	[intpdistmesh,intpdepthmesh] = meshgrid(intpdist, idepth);
																																	       
    intpvVert = nan(size(mv,3),profptnum);

	for k = 1:size(mv,3)
        v=squeeze(mv(:,:,k));
        intpvVert(k, 1:profptnum) = interp2(YSIM(1,:,end), XSIM(:,1,end), ...
            v, flon, flat, 'spline'); 
    end
    clear intpvVert2;
    intpvVert2 = interp2(distmesh,depthmesh, intpvVert, intpdistmesh,intpdepthmesh, 'spline');

    if plottopo
        clear tx ty tz td;
        tx=flon-360;ty=flat;td=intpdist;
        disp('  Getting topographic line ...');
%         [tx,ty,tz,td]=surface2profile(topofile,slon-360,slat,elon-360,elat,0.1,0.1,round(profptnum/2),'yes','other');
        tz=interp2(tlon,tlat,topogrid,flon-360,flat); %this is much faster.
        tz=tz/8000; %normalization with 8000, maximum amplitude of the topography/bathymetry in the region.
        %tz=-1.0*(mindepth-toposhift)*tz/1000+(mindepth-toposhift);
        tz(tz>0)=tz(tz>0)*1.2; %force to scale the land topography, otherwise the trench is too deep to make the
        %land topography visible.
        tz(tz<0)=tz(tz<0)*.9;
        tz=-topomax*abs(mindepthaxis-mindepth)*tz; %+(mindepth-toposhift);
        
    end
    
    figure('Position',[400 400 figwidthscale*(maxdepth-mindepth)*max(intpdist)/(vert_exaggeration*(maxdepth-mindepth)) ...
        figheightscale*(maxdepth-mindepth)]), hold on, box on, axis on;
    
    %
    if projectsurface
        clear mx my mz md ;
        mz=cell(length(xyzgrid));
        mx=flon-360;my=flat;md=intpdist;
        for pp=1:length(xyzgrid)
            disp(['  Getting xyz projection ... ',num2str(pp)]);
            mz{pp}=interp2(xyz_x{pp},xyz_y{pp},xyzgrid{pp},flon-360,flat);
            mz{pp}(mz{pp}< mindepth)=nan;
        end
    end

    if project_quakes
        clear qdist qlon qlat qdep qmag;
        disp('  Projecting earthquakes ...');
       [qdist,~,qlon,qlat,qdep,qmag]=project_points2profile(quakedata_xyz,...
           slon-360,slat,elon-360,elat,quakes_proximity,'yes','depth');
       qdep(qdep<mindepth)=nan;
    end

    if project_stations
        clear stdist stlon stlat stele;
        disp('  Projecting earthquakes ...');
       [stdist,~,stlon,stlat,stele]=project_points2profile(station_xyz,...
           slon-360,slat,elon-360,elat,station_proximity,'yes','elevation');
    end
    
    if strcmp(xlabelunit,'km')
        if strcmp(imagetype,'image')
            h2=image(intpdist, idepth, intpvVert2);
            h2.CDataMapping='scaled';
        elseif strcmp(imagetype,'contourf')
            contourf(intpdist, idepth, intpvVert2,imagecontourlevels,'linecolor','none');
        end
        plot([intpdist(1) intpdist(end)],[mindepth,mindepth],'k','linewidth',.5);
        
        daspect([1 1/vert_exaggeration 1]);

        if iplot==0
            if plotcontour
                [C, h]=contour(intpdist,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
            end
        elseif iplot==1
            if plotcontour
            [C, h]=contour(intpdist,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
            clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
            end
        end
        if project_quakes
           scatter(qdist,qdep,quakesizescale*qmag.^3,'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth);
        end
        
        if projectsurface
            for pp=1:length(xyzgrid)
                plot(md(:,1),mz{pp},projectsurface_line{pp},'linewidth',projectsurface_linewidth(pp)); 
            end
        end
        xlabeltag='Distance (km)';
        if plottopo
           %plot(td,tz,'k-','linewidth',1);
           seisplot(td,tz,mindepth-toposhift,'+-',topocolor);
%            area(td,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
           axis([min(qdist) max(qdist) mindepthaxis maxdepth]);
        else
            axis([min(qdist) max(qdist) mindepth maxdepth]);
        end
        if project_stations
           scatter(stdist,interp1(td,tz,stdist),stationsize,stationsymblecolor,'filled','linewidth',1); 
        end
    elseif strcmp(xlabelunit,'deg')
        if strcmp(xlabeldegtype,'auto')
            if (az>=45 && az<=135) || (az >=225 && az<=315 )%range(flon) >= range(flat)
                if strcmp(imagetype,'image')
                    h2=image(flon-360, idepth, intpvVert2);
                    h2.CDataMapping='scaled';
                elseif strcmp(imagetype,'contourf')
                    contourf(flon-360, idepth, intpvVert2,imagecontourlevels,'linecolor','none');
                end
                plot([flon(1)-360 flon(end)-360],[mindepth,mindepth],'k','linewidth',.5);
                
                %axis([min(flon-360) max(flon-360) -12 maxdepth]);
                daspect([1 (max(intpdist)/range(flon))/vert_exaggeration 1]);
                
                if iplot==0 && plotcontour
                    [C, h]=contour(flon-360,idepth,...
                        intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                    clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
                elseif iplot==1 && plotcontour
                    [C, h]=contour(flon-360,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                    clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
                end
                if project_quakes
                   scatter(qlon,qdep,quakesizescale*qmag.^3,...
                       'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth); 
                end
                
                if projectsurface
                    for pp=1:length(xyzgrid)
                        plot(mx,mz{pp},projectsurface_line{pp},'linewidth',projectsurface_linewidth(pp)); 
                    end
                end
                xlabeltag='Longitude (degree)';
                if plottopo
                   %plot(tx,tz,'k-','linewidth',1);
                   seisplot(tx,tz,mindepth-toposhift,'+-',topocolor);
%                    area(tx,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
                    axis([min(flon)-360 max(flon)-360 mindepthaxis maxdepth]);
                else
                    axis([min(flon)-360 max(flon)-360 mindepth maxdepth]);
                end
                if project_stations
                   scatter(stlon,interp1(tx,tz,stlon),stationsize,stationsymblecolor,'filled','linewidth',1); 
                end
                if flon(end) < flon(1)
                   set(gca,'XDir','reverse'); 
                end
            else
                if strcmp(imagetype,'image');
                    h2=image(flat, idepth, intpvVert2);
                    h2.CDataMapping='scaled';
                elseif strcmp(imagetype,'contourf')
                    contourf(flat, idepth, intpvVert2,imagecontourlevels,'linecolor','none');
                end
                plot([flat(1) flat(end)],[mindepth,mindepth],'k','linewidth',.5);
                %axis([min(flat) max(flat) 0 maxdepth]);
                daspect([1 (max(intpdist)/range(flat))/vert_exaggeration 1]);
                
                if iplot==0 && plotcontour
                    [C, h]=contour(flat,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                    clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
                elseif iplot==1 && plotcontour
                    [C, h]=contour(flat,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                    clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
                end
                if project_quakes
                   scatter(qlat,qdep,...
                       quakesizescale*qmag.^3,'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth); 
                end
                
                if projectsurface
                    for pp=1:length(xyzgrid)
                        plot(my,mz{pp},projectsurface_line{pp},'linewidth',projectsurface_linewidth(pp)); 
                    end
                end
                xlabeltag='Latitude (degree)';
                %clabel(C,h,'LabelSpacing',.2);
                if plottopo
%                    plot(ty,tz,'k-','linewidth',1);
                    seisplot(ty,tz,mindepth-toposhift,'+-',topocolor);
%                     area(ty,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
                    axis([min(flat) max(flat) mindepthaxis maxdepth]);
                else
                    axis([min(flat) max(flat) mindepth maxdepth]);
                end
                if project_stations
                   scatter(stlat,interp1(ty,tz,stlat),stationsize,stationsymblecolor,'filled','linewidth',1); 
                end
                if flat(end) < flat(1)
                   set(gca,'XDir','reverse'); 
                end
            end
        elseif strcmp(xlabeldegtype,'lon')
            if strcmp(imagetype,'image')
                h2=image(flon-360, idepth, intpvVert2);
                h2.CDataMapping='scaled';
            elseif strcmp(imagetype,'contourf')
                contourf(flon-360, idepth, intpvVert2,imagecontourlevels,'linecolor','none');
            end
            plot([flon(1)-360 flon(end)-360],[mindepth,mindepth],'k','linewidth',.5);
            %axis([min(flon-360) max(flon-360) 0 maxdepth]);
            daspect([1 (max(intpdist)/range(flon))/vert_exaggeration 1]);
            
            if iplot==0 && plotcontour
                [C, h]=contour(flon-360,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
            elseif iplot==1 && plotcontour
                [C, h]=contour(flon-360,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
            end
            if project_quakes
               scatter(qlon,qdep,...
                   quakesizescale*qmag.^3,'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth);
            end
            
            if projectsurface
                for pp=1:length(xyzgrid)
                    plot(mx,mz{pp},projectsurface_line{pp},'linewidth',projectsurface_linewidth(pp)); 
                end
            end
            xlabeltag='Longitude (degree)';
            if plottopo
               %plot(tx,tz,'k-','linewidth',1);
               seisplot(tx,tz,mindepth-toposhift,'+-',topocolor);
%                area(tx,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
                axis([min(flon)-360 max(flon)-360 mindepthaxis maxdepth]);
            else
                axis([min(flon)-360 max(flon)-360 mindepth maxdepth]);
            end
            if project_stations
               scatter(stlon,interp1(tx,tz,stlon),stationsize,stationsymblecolor,'filled','linewidth',1); 
            end
            if flon(end) < flon(1)
               set(gca,'XDir','reverse'); 
            end
        elseif strcmp(xlabeldegtype,'lat')
            if strcmp(imagetype,'image');
                h2=image(flat, idepth, intpvVert2);
                h2.CDataMapping='scaled';
            elseif strcmp(imagetype,'contourf')
                contourf(flat, idepth, intpvVert2,imagecontourlevels,'linecolor','none');
            end
            plot([flat(1) flat(end)],[mindepth,mindepth],'k','linewidth',.5);
            %axis([min(flat) max(flat) 0 maxdepth]);
            daspect([1 (max(intpdist)/range(flat))/vert_exaggeration 1]);
            
            if iplot==0 && plotcontour
                [C, h]=contour(flat,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
            elseif iplot==1 && plotcontour
                [C, h]=contour(flat,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                clabel(C,h,'FontSize',clabelsize,'Color',clabelcolor,'labelspacing', clabelspace);
            end
            if project_quakes
               scatter(qlat,qdep,...
                   quakesizescale*qmag.^3,'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth); 
            end
            
            if projectsurface
                for pp=1:length(xyzgrid)
                    plot(my,mz{pp},projectsurface_line{pp},'linewidth',projectsurface_linewidth(pp)); 
                end
            end
            xlabeltag='Latitude (degree)';
            if plottopo
               %plot(ty,tz,'k-','linewidth',1);
               seisplot(ty,tz,mindepth-toposhift,'+-',topocolor);
%                area(ty,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
                axis([min(flat) max(flat) mindepthaxis maxdepth]);
            else
                axis([min(flat) max(flat) mindepth maxdepth]);
            end
            if project_stations
               scatter(stlat,interp1(ty,tz,stlat),stationsize,stationsymblecolor,'filled','linewidth',1); 
            end
            if flat(end) < flat(1)
               set(gca,'XDir','reverse'); 
            end
        else
            error('***ERROR: wrong value for [xlabeldegtype]');
        end
    else
        error('***ERROR: wrong value for [xlabelunit]');
    end

%     title(strcat(proflabel{np,1}, ' - ',proflabel{np,2},' (',num2str(slon-360),', ',num2str(slat),') to (',num2str(elon-360),...
%         ', ',num2str(elat),')'),'FontSize',14);
    title(strcat(proflabel{np,1}, ' - ',proflabel{np,2}),'FontSize',13);
    %set(gca,'XTick',[0:100:max(intpdist)],'XTickLabel',[0:100:max(intpdist)]);
    set(gca,'YTick',10:depthlabel_interval:maxdepth,'YTickLabel',10:depthlabel_interval:maxdepth);
	set(gca,'TickDir','out');
    set(gca,'FontSize',11);
    set(gca,'YDir','reverse');
   
    xlabel(xlabeltag,'FontSize',12,'color','k');

    ylabel('Depth (km)','FontSize',12,'color','k');
%     shading flat;

    colormap(jetwr); 
    set(gca,'CLim',colorrange);
    hcbar=colorbar;
    hcbar.TickDirection='out';
    hcbar.Ticks=colorrange(1):colorinterval:colorrange(2);
    if iplot==0
        if strcmp(velocitytag,'P')
           hcbar.Label.String='V_P (km/s)';
        elseif strcmp(velocitytag,'S')
           hcbar.Label.String='V_S (km/s)';
        elseif strcmp(velocitytag,'PR')
          hcbar.Label.String='Poisson''s Ratio';
        elseif strcmp(velocitytag,'RHO')
          hcbar.Label.String='Density (kg/m^3)';
        end
    else
        hcbar.Label.String='anomaly (%)';
    end
    hcbar.Label.Color='k';

    if iplot==0
        vlabel='';
    elseif iplot==1
        vlabel='_dv';
    end
    orient('landscape')
    set(gcf,'renderer','Painters','PaperPositionMode','auto'); 
    eval(['print -dpdf -r300 -painters ',strcat('Profile_',proflabel{np,1},'_',figtitletag,'_', ite_nm,'_',...
            num2str(slon-360),'_',num2str(slat),'_to_',num2str(elon-360),...
            '_',num2str(elat)),'_',num2str(mindepth),'_',num2str(maxdepth),'km',vlabel,'.pdf']);

end  