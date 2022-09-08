%
close all; clear all

%nk=100; %total number of vertical grid

r0=6371;
define_latlon;
% set the vertical grid spacing near the surface (within half of surface wave wavelength) less than 
% 1/3 of the horizontal grid spacing
dx=dlat*110;

rbottom=605;
% top n1-1 layers: layer thickness increases by dh1 km per layer
dz1=round(dx/3*1000)/1000;
dh1=0.0;n1=15;
r(1)=0;r(2)=dz1;
for k=3:n1
	r(k)=r(k-1)+r(k-1)-r(k-2)+dh1;
end


rlamda=dx*8; 
n2=round((rlamda-r(k))/dz1/2);

dh2=round(2*dz1/n2*10000)/10000;
for k=n1+1:n1+n2
    r(k)=r(k-1)+r(k-1)-r(k-2)+dh2;
end
dz2=r(k)-r(k-1);

n3=round((rbottom-r(k))/dz2/1.5);
dh3=round(dz2/n3*100000)/100000;

for k=n1+n2+1:n1+n2+n3
	    r(k)=r(k-1)+r(k-1)-r(k-2)+dh3;
end
nk=length(r);

iflag=0;
for k=1:nk
    if r(k) > rbottom
    if iflag == 0
    iflag=1;
    kcut=k;
    end
    end
end
r(kcut+1:end)=[];
nk=length(r)

if 1
figure
plot(r,'*');
set(gca,'ydir','reverse');
end
r=r0-r;
r=flipdim(r,2);
fid=fopen(['gridz_' num2str(dx/110) '.dat'],'w');
fprintf(fid,'%f %f %f %f\n',r);
fclose(fid);

diff(r)
