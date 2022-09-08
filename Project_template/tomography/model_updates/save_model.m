%
clear all; close all;

set_mfiles_path
set_netcdf

ite_nm = ['ite_0.025deg_02'];

%read previous model
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./updated_input_' ite_nm];
dir_media=['./updated_input_' ite_nm];
disp(['Read current model... ' ite_nm]);

id=[];subs=[];subc=[];subt=[];indxem=[];indxkp=[];
id{end+1} = 0; subs{end+1}=[1,1,1];subc{end+1}=[-1,-1,-1];subt{end+1}=[1,1,1];
               indxem{end+1}=[];
               indxkp{end+1}=[];
n=1;
[snapinfo{n}]=locate_snap(fnm_conf,id{n},'start',subs{n},'count',subc{n},'stride',subt{n});
[XSIM{n},YSIM{n},ZSIM{n}]=gather_coord(snapinfo{n},'coorddir',dir_coord);
% convert from radian to degrees
XSIM{n}=90-XSIM{n}*180/pi; %latitude
YSIM{n}=YSIM{n}*180/pi;

%define the area of plot (exclude pmls)
npml=12; %number of pml layers

mrh{n}=gather_media(snapinfo{n},'rho','mediadir',dir_media);
mmu{n}=gather_media(snapinfo{n},'mu','mediadir',dir_media);
mla{n}=gather_media(snapinfo{n},'lambda','mediadir',dir_media);
mvp{n}=((mla{n}+2*mmu{n})./mrh{n}).^0.5;
mvs{n}=(mmu{n}./mrh{n}).^0.5;

NX=size(ZSIM{1},1);
NY=size(ZSIM{1},2);
NZ=size(ZSIM{1},3); 

dep=squeeze(6371-abs(ZSIM{1}(npml,npml,:)/1000)); 
 
Y = double(XSIM{n}(:,:,:));
X = double(YSIM{n}(:,:,:));
Z= double(6371-abs(ZSIM{1}(:,:,:)/1000));
V = double(mvs{n}(:,:,:))/1000;
save NewEngland_VpVs_Gao_deg0.025_02.mat X Y Z V

disp('save velocity model'); 

fid=fopen(['NewEngland_VpVs_Gao_' num2str(ite_nm) '.dat'],'w');
fprintf(fid,'%s %s %s %s %s %s\n','lat','lon','depth','Vp','Vs','Density');
for xx=1:NX, xx
for yy=1:NY, yy

for dd=length(dep):-3:1
%for dd=length(dep):-3:50
     fprintf(fid,'%3.3f %3.3f %3.2f %1.2f %1.2f %1.2f \n',XSIM{n}(xx,yy,dd), YSIM{n}(xx,yy,dd),dep(dd),mvp{n}(xx,yy,dd)/1000,mvs{n}(xx,yy,dd)/1000,mrh{n}(xx,yy,dd)/1000);
end

end
end
fclose(fid);    

