%
clear all; close all;
% /opt/matlab/2007b/bin/matlab

set_mfiles_path
set_netcdf

ite_nm = ['ite_0.015deg_01'];
fnm_conf=['./SeisFD3D.conf_' ite_nm];

iplot = 0; % 0, absolute velocity; 1, velocity perturbation

rd=[(0:31)/31,ones(1,32)];
gn=[(0:31)/31,(31:-1:0)/31];
bl=[ones(1,32),(31:-1:0)/31];
rwb=[rd',gn',bl'];
rwb=flipud(rwb);

id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];

%read updated model (same conf and coord)
dir_media=['./updated_input_' ite_nm];
dir_coord=['./updated_input_' ite_nm];

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
mvs{n}=smooth3(mvs{n},'box',[5 5 1]);

% load Topography, ../misc/CascadiaTopo.dat
% origgrid=0.05;
% [lon,lat,elev]=textread('../misc/CascadiaTopo.dat','%f %f %f');
% Topo=-elev/1000;

% set starting and ending points for plotting vertical profile images
ptLon1 = [229];
ptLon2 = [244];
ptLat1 = [46];
ptLat2 = [46];
plongridnum = length(ptLon1);
platgridnum = length(ptLat1);
figlabel={'(a)','(b)'};

for np=1:plongridnum,
% start point (slat slon)  and end point (elat elon) for plotting vertical image profile
    slat = ptLat1(np); 
	slon = ptLon1(np); 
    elat = ptLat2(np); 
	elon = ptLon2(np);
				    
	proflength = deg2km(distance(slat, slon, elat, elon));
    nmarkpt = round(proflength/111.13);
	markptdeltadist = proflength/(nmarkpt-1);
	markptdist = zeros(1, nmarkpt);
    markptdist = 0:markptdeltadist:proflength;
	markptlon = zeros(1, nmarkpt);
	markptlat = zeros(1, nmarkpt);
												    
	% number of pts on this profile
	profptnum = 10*round((proflength/5));
	% final pts (flat flon) that are shown on the profile
	flat = slat + (0:(profptnum-1))*(elat - slat)/(profptnum - 1);
	flon = slon + (0:(profptnum-1))*(elon - slon)/(profptnum - 1);
																		       
	dist = zeros(1, profptnum);
    for i = 1:(profptnum-1)
		ptdist = deg2km(distance(flat(i), flon(i), flat(i+1), flon(i+1)));
		dist(i+1) = ptdist + dist(i);
	end
	markptlon = interp1(dist, flon, markptdist, 'cubic');
	markptlat = interp1(dist, flat, markptdist, 'cubic');
																												       
	intpdist = 0:(proflength/(profptnum-1)):proflength;
	proflength = dist(profptnum);
	depth=6371-abs(ZSIM{1}(npml,npml,:)/1000);
	[distmesh,depthmesh] = meshgrid(dist, depth);
    idepth=5:1:100;
	[intpdistmesh,intpdepthmesh] = meshgrid(intpdist, idepth);
																																	       
    intpvVert = [];

	for k = 1:size(mvs{1},3)
        v=squeeze(mvs{1}(:,:,k));
		intpvVert(k, 1:profptnum) = interp2(YSIM{1}(1,:,end), XSIM{1}(:,1,end), v, flon, flat, 'spline')/1000; 
    end
    
    intpvVert2 = interp2(distmesh,depthmesh, intpvVert, intpdistmesh,intpdepthmesh, 'spline');

    figure('Position',[400 400 600 200]), hold on, box on, axis on
    
    imagesc(intpdist, idepth, intpvVert2);
    axis([0 ceil(max(intpdist)) -3 60]) 

    % plot plate interface
    fslab = fopen('../misc/CascadeSlabContours.dat','r');  % ?
    XYZ = fscanf(fslab,'%f %f %f',[3 inf]);
    XYZ = XYZ';
    xx = XYZ(:,1); % the west hemisphere
    yy = XYZ(:,2);
    zzc = XYZ(:,3);
    [XC YC] = meshgrid(235:0.25:max(xx),ptLat1(np) );
    ZC=griddata(xx,yy,zzc,XC,YC);
    slabpts = deg2km(distance(slat, slon, elat, 235:0.25:max(xx)));	
    
    plot(slabpts, ZC,'b-','LineWidth',2);
    text(360,52,'subducting slab','FontSize',16,'Rotation',-30,'Color','w')

    % plot seafloor
    [XC YC] = meshgrid(slon:0.25:max(xx),slat );
    Eelev=griddata(lon,lat,Topo,XC,YC);
    topopts = deg2km(distance(slat, slon, elat, slon:0.25:max(xx)));	

    plot(topopts, Eelev,'k-','LineWidth',2);
    
    % plot the oceanic Moho, assuming the crustal thickness ~ 5 km;
    oceanM = 6; Tslab = 35;
    dd=find(Eelev>1.85);
    
    ss = find(ZC<=25);

    text(70, 13, 'Oceanic Moho','FontSize',13,'Color','w');
    text(510, 34, 'Continental Moho','FontSize',13,'Color','k');
    
    plot([topopts(dd) slabpts(1:ss) slabpts(ss+1:end)],[Eelev(dd)+oceanM ZC(1:ss)+oceanM/cosd(20) ZC(ss+1:end)+oceanM/cosd(40) ],'m--','LineWidth',2);
        
    [tlon,tlat]=textread('../misc/JdFTrenchPoints.txt','%f %f');
    tlon = tlon+360;
    idx=find( abs(tlat-slat)==min(abs(tlat-slat)) );
    tdist = deg2km(distance(slat, slon, ptLat1(np),tlon(idx)));
    plot([tdist tdist],[0 10],'b','LineWidth',2);
    text(tdist+5, 7, 'Trench','Color','k','FontSize',14);
    
    % Axial seamount, Cobb-Eickelberg chain
    if np==3
    clon=230; clat=46.06;
    cdist = deg2km(distance(slat, slon, clat, clon));
    plot(cdist,3.0,'k^','MarkerSize',12,'MarkerFaceColor','k');
    end
    
    title(['(b) W-E profile along Latitude ' num2str(ptLat1(np))],'FontSize',16);
    
    set(gca,'XTick',[0:100:max(intpdist)],'XTickLabel',[0:100:max(intpdist)]);
    set(gca,'YTick',[0:20:100],'YTickLabel',[0:20:100]);
	set(gca,'TickDir','out');
    set(gca,'FontSize',13);
    set(gca,'YDir','reverse');
    %xlabel('Longitude (deg.)','FontSize',17);
    xlabel('Distance (km)','FontSize',17);
    ylabel('Depth (km)','FontSize',17);
    shading interp;
    colormap('jetwr'); 
    %colormap(rwb);
    if iplot==0
    %caxis([3.2 5.0]);
    caxis([3.0 4.5]);
    elseif iplot==1
    caxis([8 8])
    end
%     h=colorbar;
%     set(h,'fontsize',13);
   
 
%%%% save figure
set(gcf,'PaperPositionMode','auto');   
switch iplot
    case 0
        figname = ['VelModel_EW_profile_OBS_' num2str(ptLat1(np)) '.eps']
    case 1
        figname = ['VelModel_Percent_EW_profile_OBS_' num2str(ptLat1(np)) '.eps']
end
eval(['print -depsc ' figname])
pause
close all
 
end

