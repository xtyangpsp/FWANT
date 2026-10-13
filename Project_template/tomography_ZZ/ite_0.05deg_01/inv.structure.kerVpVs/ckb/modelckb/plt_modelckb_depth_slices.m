function plt_modelckb_depth_slices(testnamebase,depthindex)
%plot simple along-grid vertical profile and slices.
%clear all;
state=[];
load us_states;
colorlimit=[-6 6];
%nz=depthindex;
mx=160;
my=106;
mz=25;
fnm_blk='../../masterfiles/block.160x106x25.1x1x1.1x1x1.nc';

X=nc_varget(fnm_blk,'x');
Y=nc_varget(fnm_blk,'y');
Z=nc_varget(fnm_blk,'z');

x=90-reshape(X,[mx,my,mz])/pi*180;
y=reshape(Y,[mx,my,mz])/pi*180;
R=reshape(Z,[mx,my,mz]);
z=(6371-reshape(Z,[mx,my,mz])/1e3);

% minlon=min(y(1,:,1))-360;
% maxlon=max(y(1,:,1))-360;
% minlat=min(x(:,1,1));
% maxlat=max(x(:,1,1));
minlon=-76; maxlon=-68;
% minlon=minlon+360; maxlon=maxlon+360; 
minlat=40;maxlat=46;
X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
% convert to Cartesian coord
[yc,xc,zc]=sph2cart(Y,pi/2-X,R);
%%
%testnamebase='crust_uppermantle_sandwich';
%testnamebase='am_sne_low'
clear vck;
vck=reshape(load([testnamebase '.txt']),[mx my mz]);
datafile=['./result.1th.' testnamebase '/try.damp2.smot2.st0.ev0.lo0.dat'];
ncmp=2;
clear dataraw datavs;
dataraw=reshape(load(datafile),[mx my mz ncmp]);
datavs=squeeze(dataraw(:,:,:,2));
datavs=smooth3(datavs,'box',[11 11 1]);
%%
%E-W profile for check model input.
%nx=64;
nfig=0;
for i=1:length(depthindex)
    nz=depthindex(i);
    if rem(i-1,2)==0
        figure('Position',[400 400 600 500]);
        nfig=nfig+1;
    end
    clear vplotck;
    vplotck=squeeze(vck(:,:,nz)*100);

    subplot(2,2,2*(i-2*(nfig - 1))-1);
    pcolor(squeeze(y(1,:,1))-360,squeeze(x(:,1,1)),vplotck);shading flat;
%     h1.CDataMapping='scaled';
    %h.AlphaData=amaskcut;
    colormap('jetwr'); 
    hcbar1=colorbar('eastoutside');
    set(gca,'CLim',colorlimit);
    set(hcbar1,'TickDirection','out','Ticks',-10:5:10,'Fontsize',12);
    hcbar1.Label.String='dV_S (%)';
    %xlabel('Longitude (degree)'); ylabel('Depth (km)');
    %colorbar;
    hold on;
    for sb=1:length(state)
        plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[0 0 0],'LineWidth',1);
    end
    hold off;
    xlim([minlon maxlon]);
    ylim([minlat maxlat]);
    if 2*(i-2*(nfig - 1))-1 == 1
        figlabel='(a) ';
    elseif 2*(i-2*(nfig - 1))-1 == 3
        figlabel='(c) ';
    end
    title([figlabel 'Input: depth= ' num2str(round(z(1,1,nz))) ' km']);
    set(gca,'YDir','normal','Fontsize',12);
    set(gca,'TickDir','out');
    grid off
    set(gca,'layer','top')
    daspect([1 cosd(mean([minlat maxlat])) 1]);

    %E-W profiles for inversion result.
    clear vplot_inv;
    vplot_inv=squeeze(datavs(:,:,nz))*100;

    subplot(2,2,2*(i-2*(nfig - 1)));
    pcolor(squeeze(y(1,:,1))-360,squeeze(x(:,1,1)),vplot_inv); shading flat;
%     h2.CDataMapping='scaled';
    %h.AlphaData=amaskcut;
    colormap('jetwr'); 
    hcbar2=colorbar('eastoutside');
    set(gca,'CLim',colorlimit);
    set(hcbar2,'TickDirection','out','Ticks',-10:5:10,'Fontsize',12);
    hcbar2.Label.String='dV_S (%)';
    %xlabel('Longitude (degree)'); ylabel('Depth (km)');
    %colorbar;
    hold on;
    for sb=1:length(state)
        plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[0 0 0],'LineWidth',1);
    end
    hold off;
    xlim([minlon maxlon]);
    ylim([minlat maxlat]);
    %ylim([5 100]);
    if 2*(i-2*(nfig - 1)) == 2
        figlabel='(b) ';
    elseif 2*(i-2*(nfig - 1)) == 4
        figlabel='(d) ';
    end
    title([figlabel 'Recovered: depth= ' num2str(round(z(1,1,nz))) ' km']);
    set(gca,'YDir','normal','Fontsize',12);
    set(gca,'TickDir','out');
    grid off
    set(gca,'layer','top')
    daspect([1 cosd(mean([minlat maxlat])) 1]);
end
end