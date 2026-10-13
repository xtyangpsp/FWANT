%
clear all;
%% close all;
mx=69;
my=159;
mz=52;
fnm_blk=strcat('../../00masterdatafiles/block.',num2str(mx),'x',num2str(my),'x',num2str(mz),'.1x1x1.1x1x1.nc');
X=nc_varget(fnm_blk,'x');
Y=nc_varget(fnm_blk,'y');
Z=nc_varget(fnm_blk,'z');

x=90-reshape(X,[mx,my,mz])/pi*180;
y=reshape(Y,[mx,my,mz])/pi*180;
R=reshape(Z,[mx,my,mz]);
z=(6371-reshape(Z,[mx,my,mz])/1e3);

X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
%%
topofile='etopo1_bedrock_AK2_subset_blockmean_0.05.txt';
%
load AlaskaBorder;
state=[]; state(1).polygon(:,1)=Alaska.lon;state(1).polygon(:,2)=Alaska.lat; 

%plot topography and locations of the cross-sections
clear topo;
topo=load(topofile);
dxlon=0.02; dylat=0.02;
clear tlon tlat topogrid;
tlon=min(topo(:,1)):dxlon:max(topo(:,1));
tlat=min(topo(:,2)):dxlon:max(topo(:,2));
% topogrid=griddata(topo(:,1),topo(:,2),topo(:,3),tlon,tlat','cubic');
% save('AlaskaTopoGrid_0.02.mat','topogrid');
load AlaskaTopoGrid_0.02.mat;
%%
close all;
ite_nm = ['ite_0.05deg_05'];
testbaselist={'layer20km','layer30km','layer40km'};%,'layer20km','layer30km'};
% testbaselist={'forearc5'};
resulttype='input'; %input or output
if strcmp(resulttype,'input')
    caxislimit=[-10 10];
    cint=5;
else
    caxislimit=[-7.5 7.5];
    cint=2.5;
end
%
maparea.lon=[-158,-139]; %  
maparea.lat=[57, 65];
%the following locations work for plotting segmentsDVG, which was designed
%to following the orientation of profile B-B'.
% ptLon1 = [-156.6,   -153.5,   -149.0, -142.3]; ptLon1=ptLon1+360;
% ptLon2 = [-150.0,   -147.5,   -146.5, -145.6];ptLon2=ptLon2+360;
% ptLat1 = [59.9,       64.2,     64.65,   63.0];
% ptLat2 = [58.7,       60.5,     60.30,   60.0];

%the following profiles work for velocity profiles as shown in the paper.
ptLon1 = [-156.6,   -152.7,   -149.15,   -142.3]; ptLon1=ptLon1+360;
ptLon2 = [-150.0,   -147.5,   -146.5,   -145.6];ptLon2=ptLon2+360;
ptLat1 = [60.1,     63.8,     64.85,  63.0];
ptLat2 = [58.6,     60.5,     60.30,  60.0];
%subplotpar=[3,1];
plongridnum = length(ptLon1);
platgridnum = length(ptLat1);

proflabel={'A','A'' ';'B','B'' ';'C','C'' ';'D','D'' ';'E','E'' ';'F','F'' ';'G','G'' '};
fout_xsection_locations='vertical_profile_locations4gmt.txt';
fout_xsection_labels='vertical_profile_labels4gmt.txt';

fido_locations=fopen(fout_xsection_locations,'w');
fido_labels=fopen(fout_xsection_labels,'w');
%
for np=1:plongridnum
% start point (slat slon)  and end point (elat elon) for plotting vertical image profile
    slat = ptLat1(np); 
    slon = ptLon1(np); 
    elat = ptLat2(np); 
    elon = ptLon2(np);

    proflength = deg2km(distance(slat, slon, elat, elon));

    % number of pts on this profile
    profptnum = 4*round((proflength/5));
    [flat,flon]=track2(slat,slon-360,elat,elon-360,wgs84Ellipsoid,'degrees',profptnum); %get great cirlce path
    flon=360+flon;
    % save points along the profile.
    fprintf(fido_locations,'> %s\n',char(proflabel{np}));
    for kk=1:length(flat)
        fprintf(fido_locations,'%g   %g\n',flon(kk)-360,flat(kk));
    end
    fprintf(fido_labels,'%g   %g   %d   %d   %d   %s   %s\n',...
        slon-.5,slat-.5,14,0,1,'CB',proflabel{np,1});
    fprintf(fido_labels,'%g   %g   %d   %d   %d   %s   %s\n',...
        elon+.5,elat,14,0,1,'CB',proflabel{np,2});

    % continue to plot the profiles.
    plot(flon-360,flat,'-','linewidth',3,'color','k');
    text(slon-360,slat+.1,proflabel{np,1},'fontsize',16,'color','k');
    text(elon-360,elat+.1,proflabel{np,2},'fontsize',16,'color','k');
end
fclose(fido_locations);
fclose(fido_labels);  

%****************************
% get topo data
vert_exaggeration=1.0;
mindepth=0;
mindepthaxis=-55;
maxdepth=150; %20 for P, 30 for S
figheightscale=1.1; %7 for P, 5.4 FOR S
figwidthscale=.87;
xlabelunit='deg'; %deg (for degrees) or km (for km distance)
xlabeldegtype='auto'; %auto (automatically decides which of lat or lon is used;
                      %lon: use lontitude
                      %lat: use latitude.
nmarkpt=10; %degree to km coefficient.
plotcontour=0;
if plotcontour
    contouritv=4.2:0.2:4.4;
    clinecolor=[0.5 0.5 0.5];
    clinewidth=0.5;
    clabelcolor=[.5 .5 .5];
    clabelspace=500;
end

plotmoho=0;

plottopo=0;
topocolor=[.5 .5 .5];
%get topo data
topofile='etopo1_bedrock_AK2_subset_blockmean_0.05.txt';
topo=load(topofile);
toposhift=25; % scale and shift the topo data by (mindepth - toposhift)
topomax=.8; % scale facter for the maximum topo values. 

project_quakes=0;
earthquakefile='AK_quakes_1970_gt4.5_NEIC_matlab.txt';
quakeedgecolor='k';
quakelinewidth=1.5;
quakesizescale=2;

depthidx=42; %range: 36:40

close all;
for t=1:length(testbaselist)
    testnamebase=testbaselist{t};

    % clear vck;
    % vck=reshape(load([testnamebase '.txt']),[mx my mz]);
    if strcmp(resulttype,'output')
        datafile=['./result.1th.' testnamebase '/try.damp2.smot2.st0.ev0.lo0.dat'];
        ncmp=2;
        smoth=[1 1 1];
    elseif strcmp(resulttype,'input')
        datafile=[testnamebase '.txt'];
        ncmp=1; %only Vs perturbations.
        smoth=[1 1 1];
    end
    clear dataraw datavs data0;
    data0=load(datafile);
    dataraw=reshape(data0((end-mx*my*mz*ncmp+1):end),[mx my mz ncmp]);
    if ncmp==2
        datavs=100*squeeze(dataraw(:,:,:,2));
    elseif ncmp==1
        datavs=100*squeeze(dataraw(:,:,:,1));
    end
    % size(datavs)
    mv=smooth3(datavs,'box',smoth);
    % plot mapview of the checkboard: project all
    mapviewall=zeros(size(z,1),size(z,2));
    for i=1:size(z,1)
        for j=1:size(z,2)
    %         mapviewall(i,j)=max(mv(i,j,:));
            mapviewall(i,j)=mv(i,j,depthidx);
        end
    end

    %plot map view
    figure('Position',[400 400 350 250]);
    % h0=image(tlon,tlat,topogrid); hold on;
    h0=image(y(1,:,1)-360,x(:,1,1),mapviewall); hold on;
    colormap('jetwr');
    h0.CDataMapping='scaled';
    % demcmap([-6000,2000]);
    hcbar=colorbar('eastoutside');
    set(hcbar,'YDir','normal','TickDirection','out');
    hcbar.Label.String='%';
    hcbar.FontSize=12;
    set(gca,'TickDir','out','YDir','normal','CLim',caxislimit);

    axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
    daspect([1 cosd(mean(maparea.lat)) 1]);

    for sb=1:length(state)
        plot(state(sb).polygon(:,1)*0.997-0.2, state(sb).polygon(:,2)*0.993+0.25,'color',[.0 .0 .0],'LineWidth',1.5);
    end

    for np=1:plongridnum
    % start point (slat slon)  and end point (elat elon) for plotting vertical image profile
        slat = ptLat1(np); 
        slon = ptLon1(np); 
        elat = ptLat2(np); 
        elon = ptLon2(np);

        proflength = deg2km(distance(slat, slon, elat, elon));

        % number of pts on this profile
        profptnum = 4*round((proflength/5));
        [flat,flon]=track2(slat,slon-360,elat,elon-360,wgs84Ellipsoid,'degrees',profptnum); %get great cirlce path
        flon=360+flon;

        % continue to plot the profiles.
        plot(flon-360,flat,'-','linewidth',3,'color','k');
        text(slon-360,slat+.1,proflabel{np,1},'fontsize',16,'color','k');
        text(elon-360,elat+.1,proflabel{np,2},'fontsize',16,'color','k');
    end

    title([testnamebase,': ',resulttype,' at ',num2str(z(1,1,depthidx)),' km']);
    set(gca,'FontSize',14);
    hold off;
    saveas(gca,strcat('modelckb_',testnamebase,'_',ite_nm,'_mapview_',resulttype,num2str(round(z(1,1,depthidx))),'.eps'),'epsc');
  
    
    %plot vertical profiles.
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
        markptdeltadist = proflength/(nmarkpt-1);
        markptdist = 0:markptdeltadist:proflength;

        % number of pts on this profile
        profptnum = 4*round((proflength/5));
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
        depth=abs(z(1,1,:));
        [distmesh,depthmesh] = meshgrid(dist, depth);
        idepth=mindepth:1:maxdepth;
        [intpdistmesh,intpdepthmesh] = meshgrid(intpdist, idepth);

        intpvVert = nan(size(mv,3),profptnum);

        for k = 1:size(mv,3)
            v=squeeze(mv(:,:,k));
            intpvVert(k, 1:profptnum) = interp2(y(1,:,end), x(:,1,end), ...
                v, flon, flat, 'spline'); 
        end
        clear intpvVert2;
        intpvVert2 = interp2(distmesh,depthmesh, intpvVert, intpdistmesh,intpdepthmesh, 'spline');

        if plottopo
            clear tx ty tz td;
            tx=flon-360;ty=flat;td=intpdist;
            disp('  Getting topographic line ...');
            tz=interp2(tlon,tlat,topogrid,flon-360,flat); %this is much faster.
            tz=tz/8000; %normalization with 8000, maximum amplitude of the topography/bathymetry in the region.
            tz=-topomax*abs(mindepthaxis-mindepth)*tz+(mindepth-toposhift);
            tz(tz<0)=tz(tz<0)*1.2; %force to scale the land topography, otherwise the trench is too deep to make the
            %land topography visible.
        end

        figure('Position',[400 400 figwidthscale*(maxdepth-mindepth)*max(intpdist)/(vert_exaggeration*(maxdepth-mindepth)) ...
            figheightscale*(maxdepth-mindepth)]), hold on, box on, axis on;

        %
        if plotmoho
            clear mx my mz md;
            disp('  Getting xyz1 projection ...');
            mx=flon-360;my=flat;md=intpdist;
            mz=interp2(xyz1_x,xyz1_y,xyzgrid1,flon-360,flat);
        end
        if project_quakes
            clear qdist qlon qlat qdep qmag;
            disp('  Projecting earthquakes ...');
           [qdist,~,qlon,qlat,qdep,qmag]=project_points2profile(earthquakefile,...
               slon-360,slat,elon-360,elat,25,'yes','depth');
        end

        if strcmp(xlabelunit,'km')
            h2=image(intpdist, idepth, intpvVert2);
            plot([intpdist(1) intpdist(end)],[mindepth,mindepth],'k','linewidth',.5);
            h2.CDataMapping='scaled';

            daspect([1 1/vert_exaggeration 1]);
            if plotmoho
               plot(md(:,1),mz,'w-','linewidth',2); 
            end

            if plotcontour
                [C, h]=contour(intpdist,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                clabel(C,h,'FontSize',12,'Color',clabelcolor,'labelspacing', clabelspace);
                end
            if project_quakes
               scatter(qdist,qdep,quakesizescale*qmag.^2,'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth);
            end
            xlabeltag='Distance (km)';
            if plottopo
               %plot(td,tz,'k-','linewidth',1);
               area(td,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
               axis([min(qdist) max(qdist) mindepthaxis maxdepth]);
            else
                axis([min(qdist) max(qdist) mindepth maxdepth]);
            end
        elseif strcmp(xlabelunit,'deg')
            if strcmp(xlabeldegtype,'auto')
                if (az>=45 && az<=135) || (az >=225 && az<=315 )%range(flon) >= range(flat)
                    h2=image(flon-360, idepth, intpvVert2);
                    plot([flon(1)-360 flon(end)-360],[mindepth,mindepth],'k','linewidth',.5);
                    h2.CDataMapping='scaled';
                    %axis([min(flon-360) max(flon-360) -12 maxdepth]);
                    daspect([1 (max(intpdist)/range(flon))/vert_exaggeration 1]);
                    if plotmoho
                       plot(mx,mz,'m-','linewidth',2); 
                       %plot(mx2,mz2,'m--','linewidth',2);
                    end
                    if plotcontour
                        [C, h]=contour(flon-360,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                        clabel(C,h,'FontSize',12,'Color',clabelcolor,'labelspacing', clabelspace);
                    end
                    if project_quakes
                       scatter(qlon,qdep,quakesizescale*qmag.^2,...
                           'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth); 
                    end
                    xlabeltag='Longitude (degree)';
                    if plottopo
                       %plot(tx,tz,'k-','linewidth',1);
                       area(tx,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
                        axis([min(flon)-360 max(flon)-360 mindepthaxis maxdepth]);
                    else
                        axis([min(flon)-360 max(flon)-360 mindepth maxdepth]);
                    end
                    if flon(end) < flon(1)
                       set(gca,'XDir','reverse'); 
                    end
                else
                    h2=image(flat, idepth, intpvVert2);
                    plot([flat(1) flat(end)],[mindepth,mindepth],'k','linewidth',.5);
                    h2.CDataMapping='scaled';
                    %axis([min(flat) max(flat) 0 maxdepth]);
                    daspect([1 (max(intpdist)/range(flat))/vert_exaggeration 1]);
                    if plotmoho
                       plot(my,mz,'m-','linewidth',2); 
                       %plot(my2,mz2,'m--','linewidth',2);
                    end
                    if plotcontour
                        [C, h]=contour(flat,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                        clabel(C,h,'FontSize',12,'Color',clabelcolor,'labelspacing', clabelspace);
                    end
                    if project_quakes
                       scatter(qlat,qdep,...
                           quakesizescale*qmag.^2,'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth); 
                    end

                    xlabeltag='Latitude (degree)';
                    %clabel(C,h,'LabelSpacing',.2);
                    if plottopo
    %                    plot(ty,tz,'k-','linewidth',1);
                        area(ty,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
                        axis([min(flat) max(flat) mindepthaxis maxdepth]);
                    else
                        axis([min(flat) max(flat) mindepth maxdepth]);
                    end
                    if flat(end) < flat(1)
                       set(gca,'XDir','reverse'); 
                    end
                end
            elseif strcmp(xlabeldegtype,'lon')
                h2=image(flon-360, idepth, intpvVert2);
                plot([flon(1)-360 flon(end)-360],[mindepth,mindepth],'k','linewidth',.5);
                h2.CDataMapping='scaled';
                %axis([min(flon-360) max(flon-360) 0 maxdepth]);
                daspect([1 (max(intpdist)/range(flon))/vert_exaggeration 1]);
                if plotmoho
                   plot(mx,mz,'w-','linewidth',2); 
                end
                if plotcontour
                    [C, h]=contour(flon-360,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                    clabel(C,h,'FontSize',12,'Color',clabelcolor,'labelspacing', clabelspace);
                end
                if project_quakes
                   scatter(qlon,qdep,...
                       quakesizescale*qmag.^2,'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth);
                end
                xlabeltag='Longitude (degree)';
                if plottopo
                   %plot(tx,tz,'k-','linewidth',1);
                   area(tx,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
                    axis([min(flon)-360 max(flon)-360 mindepthaxis maxdepth]);
                else
                    axis([min(flon)-360 max(flon)-360 mindepth maxdepth]);
                end
                if flon(end) < flon(1)
                   set(gca,'XDir','reverse'); 
                end
            elseif strcmp(xlabeldegtype,'lat')
                h2=image(flat, idepth, intpvVert2);
                plot([flat(1) flat(end)],[mindepth,mindepth],'k','linewidth',.5);
                h2.CDataMapping='scaled';
                %axis([min(flat) max(flat) 0 maxdepth]);
                daspect([1 (max(intpdist)/range(flat))/vert_exaggeration 1]);
                if plotmoho
                   plot(my,mz,'w-','linewidth',2); 
                end
                if plotcontour
                    [C, h]=contour(flat,idepth,intpvVert2,contouritv,'LineColor',clinecolor,'linewidth',clinewidth);
                    clabel(C,h,'FontSize',12,'Color',clabelcolor,'labelspacing', clabelspace);
                end
                if project_quakes
                   scatter(qlat,qdep,...
                       quakesizescale*qmag.^2,'markeredgecolor',quakeedgecolor,'linewidth',quakelinewidth); 
                end
                xlabeltag='Latitude (degree)';
                if plottopo
                   %plot(ty,tz,'k-','linewidth',1);
                   area(ty,tz,mindepth-toposhift,'FaceColor',topocolor,'EdgeColor','none');
                    axis([min(flat) max(flat) mindepthaxis maxdepth]);
                else
                    axis([min(flat) max(flat) mindepth maxdepth]);
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

        title(strcat(proflabel{np,1}, ' - ',proflabel{np,2},' (',num2str(slon-360),', ',num2str(slat),') to (',num2str(elon-360),...
            ', ',num2str(elat),')'),'FontSize',14);

        set(gca,'YTick',max([0 mindepth]):50:maxdepth,'YTickLabel',max([0 mindepth]):50:maxdepth);
        set(gca,'TickDir','out');
        set(gca,'FontSize',14);
        set(gca,'YDir','reverse');

        xlabel(xlabeltag,'FontSize',14);

        ylabel('Depth (km)','FontSize',14);
    %     shading flat;
        colormap('jetwr'); 
        %colormap(hot);

        set(gca,'CLim',caxislimit);
        hcbar=colorbar;
        hcbar.TickDirection='out';
        hcbar.Ticks=caxislimit(1):cint:caxislimit(2);

        hold off

        saveas(gca,strcat('modelckb_',testnamebase,'_',ite_nm,'_',num2str(slon-360),'_',num2str(slat),'_to_',num2str(elon-360),...
            '_',num2str(elat),'_',resulttype,'.eps'),'epsc');


    end  
end