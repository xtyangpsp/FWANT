% draw_kernel
%
% $Date: 2008-03-03 09:16:59 +0800 (Mon, 03 Mar 2008) $
% $Revision: 2 $
% $LastChangedBy: zhangw $

clear all

%MFILEROOT='/net/fs01/data/tibet/code/mfiles';
MFILEROOT='/home/yang/Proj/easthemi/mfiles';
path([MFILEROOT '/fun-spool'],path);

OUTPUT_ROOT=['./'];

% ----------------------- parameter -----------------------
flag_surf=0;
flag_pcolor=1;
flag_overlap = 0;

flag_subplot = 0;
flag_print = 0;
flag_jetwr = 1;

OUTPUT_ROOT='./';
fnm_conf=[OUTPUT_ROOT 'SeisFD3D.conf'];
%pnm_metric='/net/fs01/data/tibet/sim.input/';
pnm_metric='/home/yang/Proj/easthemi/ite_01/sim.input/';

%pnm_out='./IC.KMI/2002.180.06.54.42/BHR/1/T1T2.P1/';
%pnm_out='IC.KMI/2002.180.06.54.42/BHZ/3/T1T2.P1/';
%pnm_out='IC.KMI/2002.180.06.54.42/BHR/3/T1T2.P1/';
pnm_out='CD.KMI/CD.BJI/BHZ/3/T1T2.P2/';

id=1;
% [x,y,z] order
%subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
%subs=[ 101 1  1  ]; subc=[ 1   -1 -1]; subt=[ 1   1  1 ];
subs=[ 1 1 21]; subc=[ -1 -1 1]; subt=[ 1   1  1 ];
%subs=[ 1 140 1]; subc=[ -1 1 -1]; subt=[ 1 1 1 ];

var_list=[]; scl_caxis_list=[];
%var_list{end+1}='phase_Vp'; scl_caxis_list{end+1}=[-1e-15,1e-15];
var_list{end+1}='phase_Vs'; scl_caxis_list{end+1}=[-1e-15,1e-15];
%var_list{end+1}='Kbpz'; scl_caxis_list{end+1}=[-1e-11,1e-11];
%var_list{end+1}='Kbqz'; scl_caxis_list{end+1}=[-2e-11,1e-21];

% -------------------- load data --------------------------

[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

%-- get coord data
[X,Y,Z]=gather_coord(snapinfo,'coorddir',pnm_metric);
X=90-X/pi*180; Y=Y/pi*180; Z=6371e3-Z; Z=Z/1e3;
nx=size(X,1);ny=size(X,2);nz=size(X,3);
%X=permute(X,[2,1,3]); Y=permute(Y,[2,1,3]); Z=permute(Z,[2,1,3]);
%[x,y,z]=sph2cart(X,Y,Z);
%[x,y,z]=sph2cart(Y,X-pi/2,Z);

% ----------------------- plot kernel -----------------------------------
nvar=length(var_list);

for n=1:nvar

[V,varnm]=gather_dist(snapinfo,id,var_list{n},'outdir',pnm_out);
%v=permute(V,[2,1,3]);

if flag_overlap==1
   hold on
elseif flag_subplot==0 | n==1
   hid=figure;set(hid,'renderer','zbuffer');
end
if flag_subplot==1
   subplot(2,2,n)
elseif flag_subplot==2
   subplot(4,1,n)
end

if flag_surf==1
   sid=surf(squeeze(Y),squeeze(X),squeeze(Z),squeeze(V));
   %View = [-37.5 70];
else
   if nx==1
      hid=pcolor(squeeze(Y),squeeze(Z),squeeze(V));
      set(gca,'ydir','reverse');
      %hid=pcolor(squeeze(Y/pi*180),squeeze(Z),squeeze(v));
   elseif ny==1
      hid=pcolor(squeeze(X),squeeze(Z),squeeze(V));
      set(gca,'ydir','reverse');
      %hid=pcolor(squeeze(X/pi*180-90),squeeze(Z),squeeze(v));
   else
      hid=pcolor(squeeze(Y),squeeze(X),squeeze(V));
      %hid=pcolor(squeeze(Y/pi*180),squeeze(X/pi*180-90),squeeze(v));
   end
   %view(-90,90)
end % flag_surf

%axis image
shading flat;
if flag_jetwr==1
   colormap('jetwr');
end
if exist('scl_caxis_list') & ~ isempty(scl_caxis_list)
    caxis(scl_caxis_list{n});
end

colorbar('vert')

title(var_list{n})

if flag_subplot==0 & flag_print==1
   print(gcf,'-dpng',[var_list{n} '_ak135_ricker05.png']);
end

end % nvar

% -------------------- save figures ------------------------
if flag_subplot>0 & flag_print==1
   print -dpng kall_ak135_ricker05.png
end
