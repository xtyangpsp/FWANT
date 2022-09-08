% modified from the file in the directory /net/fs01/data/tibet/ite_00_v163/config/grid
% Y.S. 03/13/2010

clear all; close all;

define_latlon;

max_distance=(pi/180)*distance(minlat,minlon, maxlat,  maxlon)*6371
max_traveltime=max_distance/3
mincolat = 90-maxlat; maxcolat = 90-minlat;

dx=dlat;
nx=4*round((maxcolat-mincolat)/dx/4)
x(1:nx)=mincolat+[0:nx-1]*dx;
fid=fopen(['gridx_' num2str(dx) '.dat'],'w');
fprintf(fid,'%f %f %f %f\n',x);
fclose(fid);

% y = longitude
dy=dx./sind((minlat+maxlat)/2);
ny=4*round((maxlon-minlon)/dy/4)
y(1:ny)=minlon+[0:ny-1]*dy;
fid=fopen(['gridy_' num2str(dx) '.dat'],'w');
fprintf(fid,'%f %f %f %f\n',y);
fclose(fid);


