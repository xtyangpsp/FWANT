clear all

fid=fopen('gridz.SeisGrid.dat','r');
grdz=fscanf(fid,'%f');
fclose(fid);

grd=grdz(19:end);
fid=fopen('gridz_new.dat','w');
fprintf(fid,'%f %f %f %f\n',grd)
fclose(fid)

