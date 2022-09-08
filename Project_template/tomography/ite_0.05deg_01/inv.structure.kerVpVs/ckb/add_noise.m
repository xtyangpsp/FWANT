clear all; close all;

[Kap, Kbp, dt, weight, stno,evno,xlocn,xlocv,ylocn,ylocv,zlocn,zlocv,RL]=textread('inv_Gd_list_8x8x8.ori','%s %s %f %f %d %d %d %d %d %d %d %d %s');

nt=length(dt)
% standard deviation of travel times in second
%ntstd=0.5;
ntstd=textread('dataerror_list','%f');
% noise=(ntstd-0.8).*randn(nt,1);
% ttwn=dt+noise;
ttwn = dt-ntstd;

outfile = ['inv_Gd_list_8x8x8.nosynerr'];
fidout = fopen(outfile,'w');
for ii=1:nt,
fprintf(fidout,'%s %s %f %f %d %d %d %d %d %d %d %d %s\n',Kap{ii}, Kbp{ii}, ttwn(ii), weight(ii), stno(ii),evno(ii),xlocn(ii),xlocv(ii),ylocn(ii),ylocv(ii),zlocn(ii),zlocv(ii),RL{ii});
end
fclose(fidout);

