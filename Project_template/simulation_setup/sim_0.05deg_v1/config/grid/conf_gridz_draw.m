clear all

fnm_grid='gridz_0.025.dat';

fid=fopen(fnm_grid,'r');
r=fscanf(fid,'%f');
fclose(fid);

%nk=250;
nk=length(r);
r=flipdim(r,1);

r0=6371;
r=r0-r;

figure
nabs=30;
subplot(1,2,1)
plot(r); hold on;
%plot([nk-30:nk],r(nk-30:nk),'r');
%legend('grid depth','ABS grid')
title('Radius Grid Location');
set(gca,'ydir','reverse')
xlabel('Grid Index')
ylabel('Depth (km)')
grid on
ylim([-10,1020])
xlim([0 260])

subplot(1,2,2)
plot(diff(r),r(1:end-1)); 
title('Radius Grid Spacing');
set(gca,'ydir','reverse')
ylabel('Depth (km)')
xlabel('Grid Spacing (km)')
grid on
ylim([-10,1010])
xlim([0 12])

%r(111)=267.5+5;
%r(112:250)=[1:250-112+1]*5+r(111);

%r=r0-fliplr(r);
%fid=fopen('gridz_upto652.dat','w');
%fprintf(fid,'%f %f %f %f\n',r)
%fclose(fid)
