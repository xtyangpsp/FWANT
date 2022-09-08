% compare media in the input directory and the reference model 
% modified from draw_medis_surf_all.m, Y.Shen, 02-2011

clear all

set_mfiles_path
[fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
  dir_station,dir_out]=get_simul_path

% ----------------------- parameter -----------------------
flag_overlap = 0;

flag_km=0;
flag_print = 0;
flag_jetwr = 1;
flag_clb=0; 
flag_title=1; 
flag_emlast=1;
flag_light=1;

% --------------------- output parameter -------------------
%varnm='rho';    %scl_caxis=[0,5000];
%varnm='mu';     %scl_caxis=[0,1e10];
%varnm='lambda'; %scl_caxis=[0,1e10];
varnm='Vp';     
scl_caxis=[5.5,9.0];
%varnm='Vs';     
%scl_caxis=[0,5];
%scl_caxis=[0,4];

scl_daspect=[1 1 10];
%scl_daspect=[5 5 1];

% model dimensions
!cat SeisFD3D.conf | grep total_grids_in_x | awk '{print $5}' > mx
!cat SeisFD3D.conf | grep total_grids_in_y | awk '{print $5}' > my
!cat SeisFD3D.conf | grep total_grids_in_z | awk '{print $5}' > mz
load mx; load my; load mz;
!/bin/rm mx my mz

id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];
% a north-south profile in the middle of longitude grids
id{end+1} = 0; subs{end+1}=[1,round(my/2)+1,1];subc{end+1}=[-1,1,-1];subt{end+1}=[1,1,1];
               %indxem{end+1}=[1,100,1,1,1,200];
               indxem{end+1}=[];
               indxkp{end+1}=[];
%		indxkp{end+1}=[1,mx,1,1,round(mz/2),mz]; %upper half
nsnap=numel(id);

% -------------------- load data --------------------------
for n=1:nsnap
    [snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
    [CLAT{n},LON{n},R{n}]=gather_coord(snapinfo{n},'coorddir',dir_coord);
    nx{n}=size(CLAT{n},1);ny{n}=size(CLAT{n},2);nz{n}=size(CLAT{n},3);

    clat{n}=CLAT{n}*180/pi;
    lat{n}=90-clat{n};

    lon{n}=LON{n}*180/pi;

%    [y{n},x{n},z{n}]=sph2cart(LON{n},pi/2-CLAT{n},R{n});
    
    y{n}=lon{n};
    x{n}=lat{n};

    z{n}=-6371+R{n}/1e3;

    str_unit='m';
    if flag_km
       x{n}=x{n}/1e3;y{n}=y{n}/1e3;z{n}=z{n}/1e3; str_unit='km';
    end
end

% ---------------------- plot model ----------------------
hid=figure;
set(hid,'BackingStore','on');
set(hid,'renderer','zbuffer');
set(hid,'menubar','none');
set(hid,'toolbar','figure');
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf,'PaperUnits','points');
%set(gcf,'PaperPosition',[0 0 1024 768]);
%set(0, 'DefaultFigurePaperType', 'A4');

for n=1:nsnap

    switch varnm
    case 'Vp'
       rho=gather_media(snapinfo{n},'rho','mediadir',dir_media);
       mu=gather_media(snapinfo{n},'mu','mediadir',dir_media);
       lambda=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
       v{n}=( (lambda+2*mu)./rho ).^0.5;
       v{n}=v{n}/1e3;
    case 'Vs'
       rho=gather_media(snapinfo{n},'rho','mediadir',dir_media);
       mu=gather_media(snapinfo{n},'mu','mediadir',dir_media);
       v{n}=( mu./rho ).^0.5;
       v{n}=v{n}/1e3;
    case 'rho'
       v{n}=gather_media(snapinfo{n},varnm,'mediadir',dir_media);
       v{n}=v{n}/1e3;
    otherwise
       v{n}=gather_media(snapinfo{n},varnm,'mediadir',dir_media);
    end

    if ~ isempty(indxem{n})
       i1=indxem{n}(1);i2=indxem{n}(2);
       j1=indxem{n}(3);j2=indxem{n}(4);
       k1=indxem{n}(5);k2=indxem{n}(6);
       v{n}(i1:i2,j1:j2,k1:k2)=NaN;
    end
    if ~ isempty(indxkp{n})
       i1=indxkp{n}(1);i2=indxkp{n}(2);
       j1=indxkp{n}(3);j2=indxkp{n}(4);
       k1=indxkp{n}(5);k2=indxkp{n}(6);
       vtmp=v{n}(i1:i2,j1:j2,k1:k2);
       v{n}(:,:,:)=NaN;
       v{n}(i1:i2,j1:j2,k1:k2)=vtmp;
    end

    % matlab 2009b
    v{n}=double(v{n});x{n}=double(x{n});y{n}=double(y{n});z{n}=double(z{n});

    if flag_emlast
       sid{n}=surf(squeeze(y{n}), ...
                   squeeze(x{n}), ...
                   squeeze(z{n}), ...
                   squeeze(v{n}));
    else
       sid{n}=surf(flipdim(squeeze(y{n}),3), ...
                   flipdim(squeeze(x{n}),3), ...
                   flipdim(squeeze(z{n}),3), ...
                   flipdim(squeeze(v{n}),3));
    end

    %if n>1
    %set(sid{n},'DiffuseStrength',1.0,'SpecularStrength',0.2, ...
    %    'SpecularExponent',50, ...
    %    'SpecularColorReflectance',0.1)
    %end
    hold on
end

% -- axis daspect --
%axis image
if exist('scl_daspect'); daspect(scl_daspect); end
axis tight

% -- colormap and colorbar
%c_spec=colormap('jetwr');
%c_spec=flipud(c_spec);
%colormap(c_spec);
if flag_jetwr; colormap(jetwr); end
if exist('scl_caxis','var'); caxis(scl_caxis); end
if flag_clb, cid=colorbar; end

if flag_light
   view(-30,35)
   set(gca,'box','off');
   camlight(0,10,'local');
   lighting phong
end

% -- shading --
%shading interp;
shading flat;

xlabel('Longitude,deg.');ylabel('Latitude,deg.');zlabel('Depth,km');
if flag_title, title(varnm); end


figure(2)
clf
load /home/yang/FD3Dspher/FD3Dspher.090305.s2fr.PML/run/easthemi/config/ref_models/ak135.tvel
plot(ak135(:,1),ak135(:,2),'k-');hold on
plot(-z{1}(1,:),v{1}(1,:),'b-+');hold on
axis([0 1000 0 12]);

% -------------------- save figures ------------------------
if flag_print==1
   print(gcf,'-dpng',[varnm '_ak135.png']);
end

%figure;set(gcf,'renderer','zbuffer');
%magesc(squeeze(v{1}))
%colorbar

