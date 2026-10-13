function plt_modelckb_simple_profiles(testnamebase,nx,ny)
%plot simple along-grid vertical profile and slices.
% clear all;
maparea.lon=[-76,-68];
maparea.lat=[40, 46];

colorlimit=[-6 6];
depthlimit=[15 100];
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

X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
% convert to Cartesian coord
% [yc,xc,zc]=sph2cart(Y,pi/2-X,R);
%%
% testnamebase='crust_uppermantle_sandwich';
clear vck;
vck=reshape(load([testnamebase '.txt']),[mx my mz]);
datafile=['./result.1th.' testnamebase '/try.damp2.smot2.st0.ev0.lo0.dat'];
% datafile=[testnamebase '.txt_joint_test'];
ncmp=2;
clear dataraw datavs;
dataraw=reshape(load(datafile),[mx my mz ncmp]);
datavs=squeeze(dataraw(:,:,:,2));
% size(datavs)
datavs=smooth3(datavs,'box',[7 7 1]);
%%
%E-W profile for check model input.
% nx=64;
figure('Position',[400 400 800 350]);
clear vplotck;
vplotck=squeeze(vck(nx,:,:)*100);

subplot(2,2,1);
pcolor(squeeze(y(1,:,1))'-360,squeeze(z(1,1,:))',vplotck');
shading flat;
% h1.CDataMapping='scaled';
%h.AlphaData=amaskcut;
colormap('jetwr'); 
hcbar1=colorbar('eastoutside');
set(gca,'CLim',colorlimit);
set(hcbar1,'TickDirection','out','Ticks',-10:2:10,'Fontsize',12);
hcbar1.Label.String='dV_S (%)';
xlabel('Longitude (degree)'); ylabel('Depth (km)');
%colorbar;
ylim(depthlimit);
xlim(maparea.lon);
title(['(a) Input: E-W at latitude= ' num2str(x(nx,1,1),3)]);
set(gca,'YDir','reverse','Fontsize',12)
set(gca,'TickDir','out');
grid on;
set(gca,'layer','top')
%E-W profiles for inversion result.
clear vplot_inv;
vplot_inv=squeeze(datavs(nx,:,:))*100;

subplot(2,2,2);
pcolor(squeeze(y(1,:,1))'-360,squeeze(z(1,1,:))',vplot_inv');shading flat;
% h2.CDataMapping='scaled';
%h.AlphaData=amaskcut;
colormap('jetwr'); 
hcbar2=colorbar('eastoutside');
set(gca,'CLim',colorlimit);
set(hcbar2,'TickDirection','out','Ticks',-10:2:10,'Fontsize',12);
hcbar2.Label.String='dV_S (%)';
xlabel('Longitude (degree)'); ylabel('Depth (km)');
%colorbar;
ylim(depthlimit);
title(['(b) Recovered: E-W at latitude= ' num2str(x(nx,1,1),3)]);
set(gca,'YDir','reverse','Fontsize',12);
set(gca,'TickDir','out');
xlim(maparea.lon);
grid on;
set(gca,'layer','top')
%%
% ny=20;
% figure('Position',[400 400 700 150]);
clear vplotck;
vplotck=squeeze(vck(:,ny,:)*100);

subplot(2,2,3);
pcolor(squeeze(x(:,1,1))',squeeze(z(1,1,:))',vplotck');shading flat;
% h1.CDataMapping='scaled';
%h.AlphaData=amaskcut;
colormap('jetwr'); 
hcbar1=colorbar('eastoutside');
set(gca,'CLim',colorlimit);
set(hcbar1,'TickDirection','out','Ticks',-10:2:10,'Fontsize',12);
hcbar1.Label.String='dV_S (%)';
xlabel('Latitude (degree)'); ylabel('Depth (km)');
%colorbar;
ylim(depthlimit);
title(['(c) Input: N-S at latitude= ' num2str(y(1,ny,1)-360,3)]);
set(gca,'YDir','reverse','Fontsize',12);
set(gca,'TickDir','out');
xlim(maparea.lat);
grid on;
set(gca,'layer','top')

%N-S profile
clear vplot_inv;
vplot_inv=squeeze(datavs(:,ny,:))*100;

subplot(2,2,4);
pcolor(squeeze(x(:,1,1))',squeeze(z(1,1,:))',vplot_inv');shading flat;
% h2.CDataMapping='scaled';
%h.AlphaData=amaskcut;
colormap('jetwr'); 
hcbar2=colorbar('eastoutside');
set(gca,'CLim',colorlimit);
set(hcbar2,'TickDirection','out','Ticks',-10:2:10,'Fontsize',12);
hcbar2.Label.String='dV_S (%)';
xlabel('Latitude (degree)'); ylabel('Depth (km)');
%colorbar;
ylim(depthlimit);
title(['(d) Recovered: N-S at longitude= ' num2str(y(1,ny,1)-360,3)]);
set(gca,'YDir','reverse','Fontsize',12);
set(gca,'TickDir','out');
xlim(maparea.lat);
grid on;
set(gca,'layer','top')
end