% script to compare seismograms from models with different depths to see if
% the limited depth of similation box affects the long-period arrivals.
% Background: Long-period waves (200-400 s) from easthemi are systematically 
% earlier than predicted.  
% Test: A model with the maximum depth of 1200 km and another with maxdep = 1800.
% All other parameters the same.  If the wave in the thicker model arrives earlier, then
% the thinner model is not deep enough for the long-period waves
% Results: Indeed, at 200-400s period, Uz1800 arrives earlier than Uz1200 by ~8s (BJT2SUR).
% at 300-150 s, the two waves still have noticeable shift.
% at 200-100 s, there is no noticeable difference for Rayleigh waves.
% Conclusion: the 1200 km deep model is OK for Rayleigh waves at period less than 200 s. (03-31-2011)
 
%load BJT2SUR_d1200
load HYB2SUR_d1200
Ux1200=Ux;Uy1200=Uy;Uz1200=Uz;
%load BJT2SUR_d1800
load HYB2SUR_d1800
Ux1800=Ux;Uy1800=Uy;Uz1800=Uz;

%T=[1:1001]*4;
T=[1:801]*4;
dt=4;
N=2; 
Feff=1/dt/2;

freq=[1/200 1/100];freq=freq/Feff;

[b,a]=butter(N,freq);

Uz1200fb=filtfilt(b,a,Uz1200);
Ux1200fb=filtfilt(b,a,Ux1200);
Uy1200fb=filtfilt(b,a,Uy1200);
Uz1800fb=filtfilt(b,a,Uz1800);
Ux1800fb=filtfilt(b,a,Ux1800);
Uy1800fb=filtfilt(b,a,Uy1800);

figure(1)
clf
subplot(3,1,1)
plot(T,Uz1200fb,'k');hold on
plot(T,Uz1800fb,'r');hold on
subplot(3,1,2)
plot(T,Ux1200fb,'k');hold on
plot(T,Ux1800fb,'r');hold on
subplot(3,1,3)
plot(T,Uy1200fb,'k');hold on
plot(T,Uy1800fb,'r');hold on

% time window to measure delays
t1=2000;t2=3200;
it1=round(t1/dt); it2=round(t2/dt);
ref=Uz1800fb(it1:it2);
sig=Uz1200fb(it1:it2);
maxlag=15;
ttc=[-maxlag:maxlag]*dt;
c=xcorr(sig,ref,maxlag);
[dum ixcmax]=max(c);
phase=ttc(ixcmax)

