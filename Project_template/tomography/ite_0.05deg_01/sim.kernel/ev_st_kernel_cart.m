clear all
close all

area = [ '/net/fs01/data/yang/n.cascadia' ];
ite_num = [ 'ite_04' ];
MFILE_ROOT='/home/yang/matlab/mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);
path([MFILE_ROOT '/fileexchange'],path);

global X Y Z V

OUTPUT_ROOT=['.'];
rcd=[];

rcd{1}={1 'P2' 125.95 227.95 './TA.E07A/TA.B04A/BHZ/2/T1T2.P2' [1 1 1 1] 1};
%rcd{1}={1 'P2' 125.95 227.95 './CC.JRO/UW.TTW/BHZ/2/T1T2.P2' [1 1 1 1] 1};

KSTNM='TA.E07A'
KEVNM='TA.B04A'
STLA= 46.5585; STLO= -119.8548; STEL=561.0; STEL=STEL/1e3; 
STLA2= 48.0575 ; STLO2= -123.504; STEL2= 293.8; STEL2=STEL2/1e3;%m to km
%KSTNM='CC.JRO'
%KEVNM='UW.TTW'
%STLA= 46.2751; STLO= -122.217798; STEL=1280.0; STEL=STEL/1e3;
%STLA2= 47.69445 ; STLO2= -121.69011; STEL2= 542.0; STEL2=STEL2/1e3;%m to km

EVLA=STLA2;EVLO=STLO2;EVDP=STEL2; 
%--------------- unchanged -----------------
nrcd=numel(rcd);
nrcd=1;
%pnm_metric='/net/fs01/data/yang/n.cascadia/ite_01/sim.kernel/input/';
pnm_metric=[ area '/' ite_num '/sim.kernel/input/'];
fnm_conf=['./' 'SeisFD3D.conf'];
id=1;%snap_id
subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
scl_daspect=[1 1 110];
%scl_ylim=[45.5 49.3]; scl_xlim=[-124.85 -118.75];
scl_ylim=[45.8 49.0]; scl_xlim=[-124.35 -119.25];
scl_zlim=[-50 2];
ypln=[45.5 49.3];
xpln=[-124.85 -118.75];
%ypln=[45.0 50]; xpln=[-125 -118];
zpln=[1,50];
zval=[1];
[ANGRT,MDLA,MDLO,xd,yd,zd]=fun_plane4slice(xpln,ypln,zpln,STLA,STLO,EVLA,EVLO,zval);
%[px,py,pz]=fun_border4slice(xd,yd,zd,scl_xlim,scl_ylim,MDLA,MDLO,'hide');
[px,py,pz]=fun_border4slice(xd,yd,zd,scl_xlim,scl_ylim,scl_zlim,MDLA,MDLO,'hide');
xd0=xd; yd0=yd; zd0=zd; px0=px; py0=py; pz0=pz;
var_list{1}='phase_Vp'; var_list{2}='phase_Vs';
scl_caxis_list{1}=[-1e-13,1e-13]; scl_caxis_list{2}=[-1e-13,1e-13];
nvar=length(var_list);
% config
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
% see FD3Dtopo/*/run/n.cascadia/define_xy_coords.m
minlon=-124.85; maxlon=-118.75; minlat=45.5;  maxlat=49.3;
lat0=(minlat+maxlat)/2; lon0=(minlon+maxlon)/2; 
x0=0;y0=0; 
alpha=90;
%-- get coord data --
[X,Y,Z]=gather_coord(snapinfo,'coorddir',pnm_metric);
[Y,X]=cart2geo(X,Y,x0,y0,lat0,lon0,alpha); 
Z=-Z/1e3;scl_zlim=[-scl_zlim(2) -scl_zlim(1)];
nx=size(X,1)
ny=size(X,2)
nz=size(X,3)
Y=permute(Y,[2 1 3]);
X=permute(X,[2 1 3]);
Z=permute(Z,[2 1 3]);

for m=1:nrcd
%    close all;
    pnm_out=[OUTPUT_ROOT '/' rcd{m}{5} '/']
    FIG_ROOT=['./ker_figs/' rcd{m}{5} '/']
    fctzlim=rcd{m}{7};
    if fctzlim == 1
       zpln=[1 50]*fctzlim;
       zval=[1.0]*fctzlim;
       [ANGRT,MDLA,MDLO,xd,yd,zd]=fun_plane4slice(xpln,ypln,zpln,STLA,STLO,EVLA,EVLO,zval);
       [px,py,pz]=fun_border4slice(xd,yd,zd,scl_xlim,scl_ylim,MDLA,MDLO,'hide');
%       [px,py,pz]=fun_border4slice(xd,yd,zd,scl_xlim,scl_ylim,scl_zlim,MDLA,MDLO,'hide');
    else
       xd=xd0; yd=yd0; zd=zd0; px=px0; py=py0; pz=pz0;
    end
for n=1:nvar
    fctclim=rcd{m}{6}(n);
    scl_caxis=scl_caxis_list{n}*fctclim;
    [V,varnm]=gather_dist(snapinfo,id,var_list{n},'outdir',pnm_out);
    V=permute(V,[2 1 3]);
    V(find(X<scl_xlim(1) | X>scl_xlim(2) | Y<scl_ylim(1) | Y>scl_ylim(2)))=NaN;
%    [h,sid]=fun_slice_border(xd,yd,zd,px,py,pz, ...
%     'caxis',scl_caxis, ...
%     'xlim',scl_xlim,'ylim',scl_ylim, 'daspect', scl_daspect.*[1,1,fctzlim], ...
%     'jetwr', 'bgcolor', [0.95,0.95,0.95], ...
%     'station',[STLO,STLA],'event',[EVLO,EVLA,EVDP]);
    [h,sid]=fun_slice_border(xd,yd,zd,px,py,pz, ...
     'caxis',scl_caxis, ...
     'xlim',scl_xlim,'ylim',scl_ylim, 'zlim',scl_zlim,'daspect', scl_daspect.*[1,1,fctzlim], ...
     'jetwr', 'bgcolor', [0.95,0.95,0.95], ...
     'station',[STLO,STLA,STEL],'station2',[STLO2,STLA2,STEL2]);
    cid=colorbar
    set(cid,'position',[0.85, 0.65, 0.02, 0.15])
    view([20 30]);
end
end

%
