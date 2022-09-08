%seis3d_conf_helper:
%this is a simple script as a helper in setting up frequently modified
%parameters in seis3d configuration files.
clc;
%inline function to compute snap parameters
snappar=@(t,s,b,d) (t+1-s-b)./d; 
%t: total number of grid points; s: start index; b: number of boundary
%grids at the end; d: interval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get study area range: minlon, minlat, maxlon, maxlat
define_latlon;
gridmetadata=load('./config/grid/gridxyz_metadata.mat');
gridxsize=gridmetadata.dlat; %dlat in define_latlon.
nx=gridmetadata.nx;
ny=gridmetadata.ny;
nz=gridmetadata.nz;

%saving option for full strain tensor.
snap_T_start=[13 13 40];
snap_T_interval=[4 4 2];
snap_T_time_interval=5;

%saving option for surface velocity only.
snap_V_start=[13 13 nz];
snap_V_interval=[1 1 1];
snap_V_time_interval=5;

receiver_z=9000E3;
inline_receiver_interval_x=[0.15 0.0 0.0];
inline_receiver_interval_y=[0.0 0.15 0.0];

inlineposition='center';

receiver_geolocation=[maxlat-0.5 maxlon-0.5;minlat minlon];%;maxlat-0.5 minlon;minlat maxlon-0.5]; %we don't use colatitude here. 
%it will be converted to colatitude for the configuration file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cat grid files into one single file.
if 1
    gridfile='gridall.dat';

    gridx_nm=strcat('./config/grid/gridx_',num2str(gridxsize),'.dat');
    gridy_nm=strcat('./config/grid/gridy_',num2str(gridxsize),'.dat');
    gridz_nm=strcat('./config/grid/gridz_',num2str(gridxsize),'.dat');

    unix(['echo "# x grid" > ' gridfile]);
    unix(['echo "<x grid>" >>' gridfile]);
    unix(['cat ' gridx_nm '>>' gridfile]); 

    unix(['echo "# y grid" >> ' gridfile]);
    unix(['echo "<y grid>" >>' gridfile]);
    unix(['cat ' gridy_nm '>>' gridfile]);

    unix(['echo "# z grid" >> ' gridfile]);
    unix(['echo "<z grid>" >>' gridfile]);
    unix(['cat ' gridz_nm '>>' gridfile]);

    disp(['all grid info saved in a single file: ' gridfile]);
end
% generate snapshot lines:
%wron usually, check manually
snap_T_par=snappar([nx ny nz],snap_T_start,[12 12 0],snap_T_interval);
snap_V_par=[nx-24 ny-24 1];

disp(['snap_001 = ' num2str(snap_T_start) ' ' num2str(floor(snap_T_par)) ' ' num2str(snap_T_interval) ...
    ' ' num2str(snap_T_time_interval) ' 10000 T']);
disp(['snap_002 = ' num2str(snap_V_start) ' ' num2str(snap_V_par) ' ' num2str(snap_V_interval) ...
    ' ' num2str(snap_V_time_interval) ' 10000 V']);

% generate inline receivers;
%line_001 =  52  288 9000E3 | 0.05 0.0 0.0 | 198

if(strcmp(inlineposition,'center'))
    inline_start_x=[90-minlat mean([minlon maxlon]) receiver_z];
    inline_end_x=[90-maxlat mean([minlon maxlon]) receiver_z];
    npointsx=floor((inline_end_x-inline_start_x)./inline_receiver_interval_x)+1;
    disp(['line_001 = ' num2str(inline_start_x) ...
        ' | ' num2str(inline_receiver_interval_x) ' | ' num2str(abs(npointsx(1)))]);
    
    inline_start_y=[90-mean([minlat maxlat]) minlon receiver_z];
    inline_end_y=[90-mean([minlat maxlat]) maxlon receiver_z];
    npointsy=floor((inline_end_y-inline_start_y)./inline_receiver_interval_y)+1;
    disp(['line_002 = ' num2str(inline_start_y) ...
        ' | ' num2str(inline_receiver_interval_y) ' | ' num2str(abs(npointsy(2)))]);
end

%generate separate receivers
%recv_001 = 42.1 283.5 9000E3
[nrecv,temp]=size(receiver_geolocation);

for i=1:nrecv
    fprintf('recv_%-3.3d = %7.3f %8.3f %g\n',i,90-(receiver_geolocation(i,1)),...
        receiver_geolocation(i,2),receiver_z);
end