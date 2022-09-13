%check if the assembled kernel for the inversion blocks is OK
%modified from the same name file in /net/fs01/tibet/ite_02/inv.structure

clear all; close all

MFILE_ROOT='../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

fnm_blk='block.51x59x20.1x1x1.1x1x1.nc';
mx=51; my=59; mz=20;

x=nc_varget(fnm_blk,'x');
y=nc_varget(fnm_blk,'y');
z=nc_varget(fnm_blk,'z');

%----- in spherical coordinate (X=lat; Y=lon; D=depth)
X=90-reshape(x,[mx,my,mz])/pi*180;
Y=reshape(y,[mx,my,mz])/pi*180;
Z=(6371-reshape(z,[mx,my,mz])/1e3);

% sensitivities to Vp
%fnm_G='./G_spool/TA.E07A_TA.B04A_BHZ_2_T1T2.P2.Kap.unformat'
% sensitivities to Vs
%fnm_G='./G_spool/US.DUG_CN.BMSB_BHZ_2_T1T2.P2.Kbp.unformat'

! ls -f ./G_spool/*_T1T2.P2.Kap.unformat > Kap_list
! ls -f ./G_spool/*_T1T2.P2.Kbp.unformat > Kbp_list

afiles = textread('Kap_list','%s');
bfiles = textread('Kbp_list','%s');

nfiles = length(afiles)
for cc = 1:nfiles
fnm_Ga = char(afiles(cc))
fnm_Gb = char(bfiles(cc))

fid1 = fopen(fnm_Ga,'r');
fid2 = fopen(fnm_Gb,'r');

tmp = fread(fid);

%{
% n and v0
t=fread(fid,1,'float32');
num_blk=fread(fid,1,'int');
kv0=fread(fid,1,'float32');
t=fread(fid,1,'float32');

t=fread(fid,1,'float32');
npix=fread(fid,1,'int');
t=fread(fid,1,'float32');

% indx and kernel
t=fread(fid,1,'float32');
pos=ftell(fid);
indx=fread(fid,npix,'int',4);
stat=fseek(fid,pos+4,'bof');
ker=fread(fid,npix,'float32',4);
%}
fclose(fid1);
fclose(fid2);

V=zeros(mx,my,mz);
V1=zeros(1,mx*my*mz);
for n=1:npix
V1(indx(n))=ker(n);
end

V=reshape(V1,[mx,my,mz]);

end

