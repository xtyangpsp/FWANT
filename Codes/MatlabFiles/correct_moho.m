function moho_corrected=correct_moho(modelfile,correctmodeldepth,mohodepth,vp_ref,vs_ref)
%correct moho depth by giving a model file including the following format:
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
%mohodepth: original moho depth needs to be corrected (***NOTE: NEGATIVE DEPTH VALUE)
%vp_ref, vs_ref: reference velocities, usually should be the model used in
%migration/imaging
%corrected moho depth profile will be returned.
%
%Xiaotao Yang @ Indiana University 10/11/2015

vmodel=load(modelfile);
vp=vmodel(:,2); vs=vmodel(:,3); h=vmodel(:,5); h_total=sum(h);

vp_effective=h_total/(sum(h./vp)); vs_effective=h_total/(sum(h./vs));

h_ref=abs(correctmodeldepth);
Tps_sed=h_ref/vs_effective - h_ref/vp_effective;
Tps_ref=h_ref/vs_ref - h_ref/vp_ref;

Tps_diff=Tps_sed - Tps_ref;

moho_correction=8.8*Tps_diff;

moho_corrected=mohodepth+moho_correction;

return

end

