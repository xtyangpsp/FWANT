clear all; close all;

[Kap, Kbp, dt, weight, stno,evno,xlocn,xlocv,ylocn,ylocv,zlocn,zlocv,RL]=textread('inv_Gd_list_LVV.noise','%s %s %f %f %d %d %d %d %d %d %d %d %s');

nt=length(dt)
cc=0;
for ii=1:nt
if dt(ii)>=25
cc=cc+1;
disp(Kap(ii));
end
end


