%generate anomaly bodies
clear all;
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
% convert to Cartesian coord
% [yc,xc,zc]=sph2cart(Y,pi/2-X,R);

%%
outmodlefile='layer40kmB.txt';
vout=inv_make_anomaly(outmodlefile,mx,my,mz,...
    [1 mx 1 my mz-9 mz -0.1],...
    [1 mx 1 my mz-13 mz-10 0],...
    [1 mx 1 my mz-17 mz-14 -0.1],...
    [1 mx 1 my mz-20 mz-18 0],...
    [1 mx 1 my mz-23 mz-21 -0.1]); % 5 layer

%
%
load AlaskaBorder;
state=[]; state(1).polygon(:,1)=Alaska.lon;state(1).polygon(:,2)=Alaska.lat; 
maparea.lon=[-158,-139]; % temporary 
maparea.lat=[57, 65];
figure('Position',[400 400 700 650]);
vout=load(outmodlefile);
vck=reshape(vout,[mx my mz]);
nx=25; ny=75;nz=[41,40];
clear vplot;
vplot=squeeze(vck(nx,:,:)*100);
subplot(3,2,[1 2])
pcolor(squeeze(y(1,:,1))'-360,squeeze(z(1,1,:))',vplot');
colormap('jetwr');
set(gca,'CLim',[-10 10])
set(gca,'YDir','reverse')
colorbar;
ylim([0 150]);
grid on;
xlim(maparea.lon)
title(['latitude: ' num2str(x(nx,1,1))]);

clear vplot;
vplot=squeeze(vck(:,ny,:)*100);
subplot(3,2,[3 4])
pcolor(squeeze(x(:,1,1))',squeeze(z(1,1,:))',vplot');
colormap('jetwr');
set(gca,'CLim',[-10 10])
set(gca,'YDir','reverse')
colorbar;
ylim([0 150]);
grid on;
xlim(maparea.lat)
title(['longitude: ' num2str(y(1,ny,1)-360)]);
%
clear vplot;
vplot=squeeze(vck(:,:,nz(1))*100);
subplot(3,2,5); 
pcolor(squeeze(y(1,:,1))'-360,squeeze(x(:,1,1))',vplot);
hold on;
for sb=1:length(state)
    plot(state(sb).polygon(:,1)*0.997-0.2, state(sb).polygon(:,2)*0.993+0.25,'color',[.3 .3 .3],'LineWidth',1);
end
hold off;
colormap('jetwr');
set(gca,'CLim',[-10 10])
set(gca,'YDir','normal')
colorbar;
ylim([56 66.8]);
grid on;
xlim([-162 -135])
axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
daspect([1 cosd(mean(maparea.lat)) 1]);
title(['depth: ' num2str(z(1,1,nz(1)))]);
shading flat;

clear vplot;
% nz=36;
vplot=squeeze(vck(:,:,nz(2))*100);
subplot(3,2,6)
pcolor(squeeze(y(1,:,1))'-360,squeeze(x(:,1,1))',vplot);
hold on;
for sb=1:length(state)
    plot(state(sb).polygon(:,1)*0.997-0.2, state(sb).polygon(:,2)*0.993+0.25,'color',[.3 .3 .3],'LineWidth',1);
end
hold off;
colormap('jetwr');
set(gca,'CLim',[-10 10])
set(gca,'YDir','normal')
colorbar;
ylim([56 66.8]);
grid on;
xlim([-162 -135])
axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
daspect([1 cosd(mean(maparea.lat)) 1]);
title(['depth: ' num2str(z(1,1,nz(2)))]);
shading flat;