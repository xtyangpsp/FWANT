% draw_kernel
%
% $Date: 2008-03-03 09:16:59 +0800 (Mon, 03 Mar 2008) $
% $Revision: 2 $
% $LastChangedBy: zhangw $

clear all

MFILEROOT='../../../mfiles';
path([MFILEROOT '/fun-spool'],path);

% ----------------------- parameter -----------------------
load ../../../misc/topogrid.mat
topo=qt;

% topo lat long defined as the simulation grid
%RUN_ROOT='/net/fs01/data/yang/easthemi/ite_01/sim.station/skel/fx';
RUN_ROOT='../../sim.station/skel/fx';
fnm_conf   =[RUN_ROOT '/' 'SeisFD3D.conf'];
dir_coord  =[RUN_ROOT '/' 'input'];
dir_metric =[RUN_ROOT '/' 'input'];
dir_media  =[RUN_ROOT '/' 'input'];
dir_source =[RUN_ROOT '/' 'input.src'];
dir_station=[RUN_ROOT '/' 'input'];
dir_out    =[RUN_ROOT '/' 'output'];

id = 2; % snapshot id (see SeisFD3d.conf) 
subs=[1,1,1];subc=[-1,-1,1];subt=[1,1,1];
[snapinfosurf]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

    [CLATSURF,LONSURF,RSURF]=gather_coord(snapinfosurf,'coorddir',dir_coord);
    nxsurf=size(CLATSURF,1);nysurf=size(CLATSURF,2);nzsurf=size(CLATSURF,3);

    LATSURF=pi/2-CLATSURF;
    clatsurf=CLATSURF*180/pi;
    lonsurf=LONSURF*180/pi;
    latsurf=90-clatsurf;

Zin=topo;Zlat=LATSURF;Zlon=LONSURF;
Zclat=CLATSURF;
% define area to be plotted
%flatmin=-45*pi/180;
flatmin=Zlat(end,1);flatmax=Zlat(1,1);
%flatmin=Zclat(1,1); flatmax=Zclat(end,1);
%flonmin=-30*pi/180;
flonmin=Zlon(1,1);flonmax=Zlon(1,end);
%nindlat=find(Zlat<flatmin);
%nindlon=find(Zlon<flonmin);
%Zin(nindlat)=NaN; Zin(nindlon)=NaN;
phi_min=flatmin; phi_max=flatmax;
theta_min=flonmin; theta_max=flonmax;
Zin=double(Zin);theta_min=double(theta_min);theta_max=double(theta_max);
phi_min=double(phi_min);phi_max=double(phi_max);

% plate boundaries
!/bin/cp /home/yang/Proj/Shared/plate_boundaries/Bird_plate_boundaries.xy plate.dat
load plate.dat
plon=plate(:,1)*pi/180; plat=plate(:,2)*pi/180;
pr(1:length(plon))=6371e3; pr=pr';
nindplon=find(plon>flonmax | plon < flonmin);
plon(nindplon)=[]; plat(nindplon)=[]; pr(nindplon)=[];
nindplat=find(plat>flatmax | plat < flatmin);
plon(nindplat)=[]; plat(nindplat)=[];pr(nindplat)=[];
[py,px,pz]=sph2cart(plon,plat,pr);
py=py/1e3; px=px/1e3; pz=pz/1e3;

mx=47; my=90; mz=17;
fnm_blk='../block.47x90x17.1x1x1.1x1x1.nc'

smot_list=[16];
damp_list=[64];
zval_list=[ 1 ];
zindx_list=[10];
res_list={'./result.1th'};
cmp_list={'Vp','Vs'};

num_z=numel(zindx_list);
num_smot=numel(smot_list);
num_damp=numel(damp_list);
num_res=numel(res_list);
num_cmp=numel(cmp_list);

X=nc_varget(fnm_blk,'x');
Y=nc_varget(fnm_blk,'y');
Z=nc_varget(fnm_blk,'z');

R=reshape(Z,[mx,my,mz]);
X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
% convert to Cartesian coord
[y,x,z]=sph2cart(Y,pi/2-X,R);
x=x/1e3; y=y/1e3;; z=z/1e3; %in km

dep=6371-R(1,1,zindx_list(1))/1.e3;

% ----------------------- plot kernel -----------------------------------

    pnm_result=res_list{1};
    %pnm_fig=[pnm_result '.fig'];
    KSMOTNM=num2str(smot_list(1));
    KDAMPNM=num2str(damp_list(1));

    fnm_try=[pnm_result '/try.damp' KDAMPNM '.smot' KSMOTNM '.st0.ev0.lo0.dat']

    M=load(fnm_try);
    W=reshape(M,mx,my,mz,2);
%    W=permute(W,[2,1,3,4]);

    W=W*100; %percent

nk=1;
for ncmp=2:2
    KCMPNM=cmp_list{ncmp};
    Vel=squeeze(W(:,:,:,ncmp));

        k=zindx_list(nk);
        d=zval_list(nk);
        V=squeeze(Vel(:,:,k));
        xs=squeeze(x(:,:,k));
        ys=squeeze(y(:,:,k));
        zs=squeeze(z(:,:,k));
        KDNM=num2str(zval_list(nk));
        KINM=num2str(zindx_list(nk));
        Vmax=max(max(abs(V)));

        V=double(V);

   hid=figure;set(hid,'renderer','zbuffer');
   surf(squeeze(permute(ys,[2 1 3])), ...
            squeeze(permute(xs,[2 1 3])), ...
            squeeze(permute(zs,[2 1 3])), ...
            squeeze(permute(V,[2 1 3])));hold on
shading interp;
   colormap('jetwr');
view(180,0)
set(gca,'box','off');
   camlight(0,-45,'local');
   lighting phong
   material dull
caxis([-3 3]);
colorbar('location','manual','position',[0.85 0.2 0.02 0.2])
hold on
% depth of horizontal surface in km
%dep=(6371000-R(1,1))/1.e3;
%cR=6371*1.e3;
cR=6371;
sphere3d1c10m(Zin,theta_min,theta_max,phi_min,phi_max,cR,1,'contour','linear',1);hold on
%sphere3d1c2000m(Zin,theta_min,theta_max,phi_min,phi_max,cR,1,'contour','nearest',0);hold on

plot3(py,px,pz,'+','MarkerSize',2);hold on

view(-10,5);
axis off; axis tight;

title([cmp_list{ncmp}, ' at dpeth = ', num2str(dep,4)]); 

end
%end
%end

% -------------------- save figures ------------------------
