% script to measure phase delays between the observed and synthetic waveforms
% modified from measure_phase_delay.csh
% Further modified by Xiaotao Yang @ UMASS
% Modifications:
% 1. changed to power relation ship time window length as a funtion of
% period;
% 2. SNR-based xcorrelation shift. if SNRs for the positive and negative
% data differ too much, use the larger SNR result. Otherwise, use weighted
% by SNR values to get the final shift (maximum correlation coefficient).
% 3. Removed cycle-skipping related selection criterion. This can be
% applied later when inspecting the overall measurements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Set up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise/src'))
PROJHOME = '/depot/xtyang/data/projects/xtyang/craton';
%PROJHOME = '/Users/xtyang/Work/Research/Projects/Craton';
ite = '/ite_0.05deg_01';
wkdir = [PROJHOME ite];
egfdir = [PROJHOME '/data/data_craton/PAIRS_TWOSIDES_gaussian_a2t6'];
syndir = [wkdir '/syn.seismograms'];
outdir = [wkdir '/measure/measure_stnpair'];
plotdir = strcat(wkdir,'/measure/plots_test');
stainfo = [PROJHOME '/STinfo/craton_station_withdata.txt'];
%%%%%%%%%
% parpool('local',15);
% parpool('local',50);
%%%%%%%%%%%%%%%%%END OF BEHAVIOR CONTROL PARAMETERS
%period and frequency band
pband=[7.5 15; 10 20; 15 30; 20 40; 30 60; 40 75; 60 100; 75 125];
%%%%% define parameters %%%%
max_dV = 0.15;
max_dT = 15; % maximum delay in second
snr_cutoff = 6.0; % snr limit
xcoeff_cutoff = 0.7; % cross correlation coefficient limit

min_substack = 4; % minimum number of substacks.
% define time windows
tminimum=min(min(pband)); % shortest period in second, should be the taper window length
tmaximum=1200; % length of synthetic green's function
tmaximumegf=1200; % length of egfs
waveletshift=6; %this is the time factor in simulation and convolution.

%%% interpolate and decimate both syn and egf to the same sampling rate
dt_resample=0.2;

%%%%%%%%%
cmin=2.;cmax=5.5; % km/s group velocity
synerr=0.05*110/2/mean([cmin,cmax]); % about half of the grid spacing divided by group velocity (~0.05*110 km /2/ 4.5 km/s)
verylargenumber = 1.e9; %to check numerically unstable values. 

%%% define the filter and taper for the egf and syn
taperfraction=tminimum/tmaximumegf; % taper tminimum at the ends
syntaperfraction=tminimum/tmaximum; % taper tminimum at the ends
N=2;
%%%%%%%%%%%%%%%%%BEHAVIOR CONTROL PARAMETERS
fig_flag = 0;  % 1 == plot figure; else no figure
savefig=0;
saveresult=1; %1, save measurements to output file.
