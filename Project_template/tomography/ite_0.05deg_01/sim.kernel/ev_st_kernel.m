clear all

MFILEROOT='../../mfiles';
path([MFILEROOT '/fun-spool'],path);

global X Y Z V

OUTPUT_ROOT=['.'];
rcd=[];
%rcd{1}={1 'P1' 72.966 309.36 './IC.KMI/2002.220.11.42.07/BHZ/1/T1T2.P1' [1 1 1 1] 1};
rcd{1}={1 'P2' 72.966 309.36 './CD.LSA/CD.BJI/BHZ/3/T3T4.P2' [1 1 1 1] 1};
%rcd{1}={1 'P1' 72.966 309.36 './YA.MC07/2004.087.18.47.29/BHZ/2/T3T4.P1' [1 1 1 1] 1};
%rcd{1}={1 'P1' 72.966 309.36 './YA.MC07/2004.087.18.47.29/BHZ/2/T5T6.P1' [1 1 1 1] 1};
KSTNM='CD.LSA';
KEVNM='CD.BJI';
STLA= 29.70; STLO= 091.15;
EVLA= 40.0403; EVLO= 116.175; EVDP= 0.0;
%[KSTNM KEVNM STLA STLO EVLA EVLO]=system(sprintf('ev_st_loc.sh %s',rcd{1}(5)));
%[stat loc_str]=system(sprintf('./ev_st_loc.sh %s',char(rcd{1}(5))));
%[KSTNM remain]=strtok(loc_str);
%[KEVNM remain]=strtok(remain);
%[STLA remain]=strtok(remain);
%[STLO remain]=strtok(remain);
%[EVLA remain]=strtok(remain);
%[EVLO remain]=strtok(remain);
%[EVDP remain]=strtok(remain);
%STLA=str2num(STLA);
%STLO=str2num(STLO);
%EVLA=str2num(EVLA);
%EVLO=str2num(EVLO);
%EVDP=str2num(EVDP);
STLA= 90-STLA; 
EVLA= 90-EVLA;
EVDP=6371-EVDP;
%--------------- unchanged -----------------
nrcd=numel(rcd);
nrcd=1;
%pnm_metric='../../model.iran.volume/';
%fnm_conf=['../../event.CMT.BOB/' KEVNM '/' 'SeisFD3D.conf'];
pnm_metric='../sim.input/';
fnm_conf=['./SeisFD3D.conf'];
id=1;
%subs=[1,1,2];subc=[-1,-1,-1];subt=[3,3,2];
subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
scl_daspect=[1 1 110];
%scl_ylim=[19,39]; scl_xlim=[81,108];
scl_ylim=[19,45]; scl_xlim=[81,120];
ypln=[19,45];
xpln=[81,120];
%scl_ylim=[23,36]; scl_xlim=[86,106];
%ypln=[23-5,36+5];
%xpln=[86-5,106+5];
zpln=[0,500];
zval=[0 700 200];
%scl_daspect=[1 1 110];
%scl_xlim=[39,65.8]; scl_ylim=[22,43.8];
%xpln=[39-10,66+10];
%ypln=[22-10,44+10];
%zpln=[0,400];
%zval=[0 150 80];
[ANGRT,MDLA,MDLO,xd,yd,zd]=fun_plane4slice(xpln,ypln,zpln,STLA,STLO,EVLA,EVLO,zval);
[px,py,pz]=fun_border4slice(xd,yd,zd,scl_xlim,scl_ylim,MDLA,MDLO,'hide');
xd0=xd; yd0=yd; zd0=zd; px0=px; py0=py; pz0=pz;
var_list{1}='phase_Vp'; var_list{2}='phase_Vs';
%var_list{3}='amplitude_Vp'; var_list{4}='amplitude_Vs';
scl_caxis_list{1}=[-0.5e-14,0.5e-14]; scl_caxis_list{2}=[-0.5e-14,0.5e-14];
%scl_caxis_list{3}=[-1e-14,1e-14]; scl_caxis_list{4}=[-1e-14,1e-14];
nvar=length(var_list);
% config
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
%-- get coord data --
[Y,X,Z]=gather_coord(snapinfo,'coorddir',pnm_metric);
Y=90-Y/pi*180; X=X/pi*180; Z=6371e3-Z; Z=Z/1e3;
%Y=Y/pi*180-90; X=X/pi*180; Z=6371e3-Z; Z=Z/1e3;
for m=1:nrcd
%    close all;
    pnm_out=['./' rcd{m}{5} '/']
%    FIG_ROOT=['./ker_figs/' rcd{m}{5} '/']
    fctzlim=rcd{m}{7};  
  if fctzlim~=1
       zpln=[0,500*fctzlim];
       zval=[0 150 80]*fctzlim;
       %zval=[0];
       [ANGRT,MDLA,MDLO,xd,yd,zd]=fun_plane4slice(xpln,ypln,zpln,STLA,STLO,EVLA,EVLO,zval);
       [px,py,pz]=fun_border4slice(xd,yd,zd,scl_xlim,scl_ylim,MDLA,MDLO,'hide');
    else
       xd=xd0; yd=yd0; zd=zd0; px=px0; py=py0; pz=pz0;
    end
for n=1:nvar
    fctclim=rcd{m}{6}(n);
    scl_caxis=scl_caxis_list{n}*fctclim;
    [V,varnm]=gather_dist(snapinfo,id,var_list{n},'outdir',pnm_out);
    V=V/10^7;
    V(find(X<scl_xlim(1) | X>scl_xlim(2) | Y<scl_ylim(1) | Y>scl_ylim(2)))=NaN;
    [h,sid]=fun_slice_border(xd,yd,zd,px,py,pz, ...
     'caxis',scl_caxis, ...
     'xlim',scl_xlim,'ylim',scl_ylim, 'daspect', scl_daspect.*[1,1,fctzlim], ...
     'jetwr', 'bgcolor', [0.95,0.95,0.95], ...
     'station',[STLO,STLA],'event',[EVLO,EVLA,EVDP]);
     title([rcd{1}(5),'  ',var_list{n}])
     xlabel('Longitude (E)')
     ylabel('Latitude (N)')
     zlabel('Depth (km)')
     cid=colorbar('vert','location','SouthOutSide');
     set(cid,'xtick',scl_caxis);
     view([35 25])
    %fun_print(sid,FIG_ROOT,[KSTNM '_' KEVNM '_' varnm '_' 'border'],'png');
    %[h,sid]=fun_slice_light(xd,yd,zd, ...
    % 'caxis',scl_caxis, ...
    % 'xlim',scl_xlim,'ylim',scl_ylim, 'daspect', scl_daspect.*[1,1,fctzlim], ...
    % 'jetwr', 'bgcolor', [0.95,0.95,0.95], ...
    % 'station',[STLO,STLA],'event',[EVLO,EVLA,EVDP]);
    %fun_print(sid,FIG_ROOT,[KSTNM '_' KEVNM '_' varnm '_' 'light'],'png');
    %fun_colorbar_print(sid,FIG_ROOT,[KSTNM '_' KEVNM '_' varnm],'png','png',scl_caxis);

      scrsz=get(0,'ScreenSize');
      h=figure('position',[1 1 scrsz(3)/3 2*scrsz(4)/3]);
      set(h,'renderer','zbuffer');
      set(h,'menubar','none');
      set(h,'toolbar','figure');
    for ip=1:3
      ip
      vid=interp3(X,Y,Z,V,xd{ip},yd{ip},zd{ip});
      if (ip==1)
        subplot(4,1,1)
        pcolor(xd{ip},zd{ip},vid);shading interp;colormap('jetwr')
        daspect([1,110/abs(cos(ANGRT*pi/180)),1])
        axis ij
        axis tight
        xmin=min(STLO,EVLO)-2.5*abs(cos(ANGRT*pi/180));
        xmax=max(STLO,EVLO)+2.5*abs(cos(ANGRT*pi/180));
        zmin=min(min(zd{ip}));
        zmax=max(max(zd{ip}));
        axis([xmin xmax zmin zmax])
        xlabel('Longitude (E)')
        ylabel('Depth (km)')
      elseif (ip==2)
        subplot(4,1,2)
        pcolor(xd{ip},zd{ip},vid);shading interp;colormap('jetwr')
        daspect([1,110/abs(sin(ANGRT*pi/180)),1])
        axis ij
        axis tight
        xmin=min(STLO,EVLO)-2.5*abs(sin(ANGRT*pi/180));
        xmax=max(STLO,EVLO)+2.5*abs(sin(ANGRT*pi/180));
        zmin=min(min(zd{ip}));
        zmax=max(max(zd{ip}));
        axis([xmin xmax zmin zmax])
        xlabel('Longitude (E)')
        ylabel('Depth (km)')
      elseif (ip==3)
        subplot(4,1,[3 4])
        pcolor(X(:,:,1),Y(:,:,1),V(:,:,end));shading interp;colormap('jetwr')
        %daspect([1,110/abs(sin(ANGRT*pi/180)),1])
        axis tight
        axis equal
        xmin=min(STLO,EVLO)-2.5;
        xmax=max(STLO,EVLO)+2.5;
        zmin=min(STLA,EVLA)-2.5;
        zmax=max(STLA,EVLA)+2.5;
        axis([xmin xmax zmin zmax])
        xlabel('Longitude (E)')
        ylabel('Latitude (N)')
      end
      caxis([-1e-14 1e-14]);%colorbar;
    end
    
end
end
