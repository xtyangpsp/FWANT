close all
% PROJHOME = '/depot/xtyang/data/projects/xtyang/craton';
PROJHOME = '/Users/xtyang/Work/Research/Projects/Craton';
stainfo = [PROJHOME '/STinfo/craton_station_withdata.txt'];
ite_tag = 'ite_0.05deg_01';
rawdelaydatafile=[ite_tag,'_measure_result_craton.dat'];
%%
get_delay_latlon(stainfo,{rawdelaydatafile})
%%
delayfile=[rawdelaydatafile '_latlon'];
pband=[7.5 15; 10 20; 15 30; 20 40; 30 60; 40 75; 60 100; 75 125];
fband=flip(flip(1./pband),2);
maxdt=8;
cutoffmethod='no';
plt_data_dist_multiple({delayfile},fband,[2,4],maxdt,0.+maxdt,cutoffmethod)
%%
ray_min=20;
cmax=100;
plt_raypath(delayfile,stainfo,fband,ite_tag,ray_min,cmax)


%% create filter file for kernel calculation
% syndt=0.8; % dt of the synthetics. (saved dt, not the simulation dt).
% kern_createfilterfile(fband,syndt)

