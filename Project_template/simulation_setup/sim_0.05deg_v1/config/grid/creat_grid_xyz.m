% modified from the file in the directory /net/fs01/data/tibet/ite_00_v163/config/grid
% Y.S. 03/13/2010

clear all; close all;

define_latlon;
r0=6371;

max_distance=(pi/180)*distance(minlat,minlon, maxlat,  maxlon)*r0;
max_traveltime=max_distance/3;
mincolat = 90-maxlat; maxcolat = 90-minlat;

%%%%%%%%%%%%% x grid: latitude direction
nx=4*round((maxcolat-mincolat)/dlat/4);
x(1:nx)=mincolat+[0:nx-1]*dlat;
fid=fopen(['gridx_' num2str(dlat) '.dat'],'w');
fprintf(fid,'%f %f %f %f\n',x);
fclose(fid);

%%%%%%%%%%%%% y grid: y = longitude
dlon=dlat./sind((minlat+maxlat)/2);
ny=4*round((maxlon-minlon)/dlon/4);
y(1:ny)=minlon+[0:ny-1]*dlon;
fid=fopen(['gridy_' num2str(dlat) '.dat'],'w');
fprintf(fid,'%f %f %f %f\n',y);
fclose(fid);

%%%%%%%%%%%% z grid: depth
dx=dlat*110;
rbottom=1000;
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
nz=length(r);

iflag=0;
for k=1:nz
    if r(k) > rbottom
        if iflag == 0
            iflag=1;
            kcut=k;
        end
    end
end
r(kcut+1:end)=[];
nz=length(r);

if 1
    figure
    plot(r,'*');
    set(gca,'ydir','reverse');
end
r=r0-r;
r=flip(r,2);
fid=fopen(['gridz_' num2str(dx/110) '.dat'],'w');
fprintf(fid,'%f %f %f %f\n',r);
fclose(fid);

nx,ny,nz

%save grid metadata.
save('gridxyz_metadata.mat','dlat','minlat','minlon', 'maxlat', 'maxlon',...
    'max_distance','max_traveltime','nx','ny','nz','x','y','r','rbottom');