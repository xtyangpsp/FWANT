function gravity_corrected=correct_gravity(modelfile,correctmodeldepth,gravity,rou_ref)
%correct gravity anomaly by giving a model file including the following format:
%layer_num vp vs density thickness ? ? ? ?
% 1 5.5 2.7 2.65 1.03366
% 2 3.3 1.9 2.4 0.0888223
% 3 5 2.7 2.65 0.623239
% 4 5 2.2 2.65 0.0693071
% 5 6.5 3.5 2.75 2.34026
% 6 6 3.3 2.6 0.208585
% 7 5.5 3.2 2.65 0.311623
%
%USAGE:
%modelfile: filename of the model file
%correctmodeldepth: depth profile of used for correction
%gravity: original gravity needs to be corrected
%rou_ref: reference density used in calculating the Bouguer gravity anomaly
%
%** corrected gravity will be returned.
%
%Xiaotao Yang @ Indiana University 10/14/2015

%correct gravity for sediments
%https://en.wikipedia.org/wiki/Bouguer_anomaly
G=6.67e-11;
%rou_ref=2.67e+3;

vmodel=load(modelfile);
rou=1000*vmodel(:,4);
h=vmodel(:,5); h_total=sum(h);
rou_effective=sum(rou.*h)/h_total;

g_correction_gradient=2*pi*G*(rou_effective - rou_ref)*1000;
g_correction_gradient=g_correction_gradient*100; %convert to Gal
g_correction_gradient=g_correction_gradient*1000; %convert to mGal

gravity_corrected=gravity - g_correction_gradient*abs(correctmodeldepth);

return

end

