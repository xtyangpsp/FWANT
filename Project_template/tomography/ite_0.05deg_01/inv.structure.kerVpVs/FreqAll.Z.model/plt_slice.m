
clear all;
close all;
addpath(genpath('/depot/xtyang/data/codes/mexcdf'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise'))
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))
MFILE_ROOT='../../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

% -- parameter --
mx=8;my=6;mz=8;
fnm_blk=['../block.' num2str(mx) 'x' num2str(my) 'x' num2str(mz) '.1x1x1.1x1x1.nc'];

smot_list=[1.6];
damp_list=[1.6];
zindx_list=[8 7 6 5 4 3 2 1];
res_list={'./result.1th'};

cmp_list={'Vp','Vs'};

savefig_flag=1;

% ----------------------------------
num_z=numel(zindx_list);
num_smot=numel(smot_list);
num_damp=numel(damp_list);
num_res=numel(res_list);

X=nc_varget(fnm_blk,'x');
Y=nc_varget(fnm_blk,'y');
Z=nc_varget(fnm_blk,'z');

x=90-reshape(X,[mx,my,mz])/pi*180;
y=reshape(Y,[mx,my,mz])/pi*180;
R=reshape(Z,[mx,my,mz]);
z=(6371-reshape(Z,[mx,my,mz])/1e3);

X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
% convert to Cartesian coord
[yc,xc,zc]=sph2cart(Y,pi/2-X,R);
xc=xc/1e3; yc=yc/1e3; zc=zc/1e3; %in km

minlon=min(y(1,:,1))-360; maxlon=max(y(1,:,1))-360; minlat=min(x(:,1,1));maxlat=max(x(:,1,1));
%% stations

% ----------------------------------
%for nres=1:num_res
for nres=1
    pnm_result=res_list{nres};
    pnm_fig=[pnm_result '.fig'];
    if ~ isdir(pnm_fig)
       mkdir(pnm_fig);
    end
for nsmot=1:num_smot
    KSMOTNM=num2str(smot_list(nsmot));
for ndamp=1:num_damp
    KDAMPNM=num2str(damp_list(ndamp));

    fnm_try=[pnm_result '/try.damp' KDAMPNM '.smot' KSMOTNM '.st0.ev0.lo0.dat'];

    M=load(fnm_try);
    W=reshape(M,mx,my,mz,2);
%    W=permute(W,[2,1,3,4]);

    W=W*100; %percent

% z slice
for ncmp=1:length(cmp_list)
    KCMPNM=cmp_list{ncmp};
    if strcmp(KCMPNM, 'Vp')
        Vel=squeeze(W(:,:,:,1));
    elseif strcmp(KCMPNM, 'Vs')
        Vel=squeeze(W(:,:,:,2));
    else
        error('***ERROR: Wrong component tag, only: Vp and/or Vs are valid!');
    end

    for nk=1:num_z
%    for nk=8
        k=zindx_list(nk);
%         d=zval_list(nk);
        V=squeeze(Vel(:,:,k));
		% remove the mean, 04/26/2011, Y.Shen
		V=V-mean(mean(V));
		%
	XS=squeeze(x(:,:,k));
	YS=squeeze(y(:,:,k));
	zs=squeeze(zc(:,:,k));
%         KDNM=num2str(zval_list(nk));
        KINM=num2str(zindx_list(nk));
        Vmax=max(max(abs(V)));

    if strcmp(KCMPNM, 'Vp')
    %fnm_fig=['figure_Vp_' int2str(nk) '.eps']
    fnm_fig=['figure_Vp_s' num2str(smot_list(nsmot)) 'd' num2str(damp_list(ndamp)) '_z' num2str(z(round(mx/2),round(my/2),k),2)];
    fnm_outfile=['Vp_' int2str(nk) '.txt'];
    Vmax=10;
    end;
    if strcmp(KCMPNM, 'Vs')
    fnm_fig=['figure_Vs_s' num2str(smot_list(nsmot)) 'd' num2str(damp_list(ndamp)) '_z' num2str(z(round(mx/2),round(my/2),k),2)];
    fnm_outfile=['Vs_' int2str(nk) '.txt'];
    Vmax=10;
    end;

hid=figure('Position',[400 -200 700 700]); hold on, box on, axis on
set(hid,'renderer','zbuffer');
set(gca,'TickDir','out');

pcolor(YS-360,XS,V);

shading flat
colormap('jetwr');

daspect([1,cosd((minlat+maxlat)/2),1]);
xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
axis([minlon maxlon minlat maxlat]);
colorbar

title(['smooth: ' num2str(smot_list(nsmot)) ', damp: ' num2str(damp_list(ndamp)) ','...
    cmp_list{ncmp} ', z: ',num2str(z(round(mx/2),round(my/2),k),2),' km'])
caxis([-Vmax,Vmax]);

if savefig_flag
    print('-dpng',['./' pnm_fig '/' fnm_fig '.png']);
    close;
end

    end
end

end
end
end
