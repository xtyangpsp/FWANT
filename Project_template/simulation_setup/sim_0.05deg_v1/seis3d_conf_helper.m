%seis3d_conf_helper:
%this is a simple script as a helper in setting up frequently modified
%parameters in seis3d configuration files.
clc;
%inline function to compute snap parameters
snappar=@(t,s,b,d) floor((t+d-s-b)./d); 
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

%boundary layers
sim_boundary_layers=[12 12 0];
%saving option for full strain tensor.
snap_T_start=[sim_boundary_layers(1)+3 sim_boundary_layers(2)+3 22];
snap_T_interval=[10 10 5];
snap_T_time_interval=5;

allow_missing_grids=1; %THIS, IF TRUE (1), WILL TOLERATE THE PARAMETERS WHEN SOME GRIDS ARE NOT SAVED.
%saving option for surface velocity only.

snap_V_start=[sim_boundary_layers(1)+1 sim_boundary_layers(2)+1 nz];
snap_V_interval=[1 1 1];
snap_V_time_interval=snap_T_time_interval;

receiver_z=9000E3;
inline_receiver_interval_x=[dlat*2 0.0 0.0];
inline_receiver_interval_y=[0.0 dlat*2 0.0];

inlineposition='center';

receiver_geolocation=[maxlat-dlat*10 maxlon-dlat*10;minlat minlon];%;maxlat-0.5 minlon;minlat maxlon-0.5]; %we don't use colatitude here. 
%it will be converted to colatitude for the configuration file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cat grid files into one single file.
if 1
    gridfile='gridall.dat';

    gridx_nm=strcat('./config/grid/gridx_',num2str(gridxsize,'%.5f'),'.dat');
    gridy_nm=strcat('./config/grid/gridy_',num2str(gridxsize,'%.5f'),'.dat');
    gridz_nm=strcat('./config/grid/gridz_',num2str(gridxsize,'%.5f'),'.dat');

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
%wrong usually, check manually
snap_T_par=snappar([nx ny nz],snap_T_start,sim_boundary_layers,snap_T_interval);
snap_T_par
%check the snap_T_par
totalx = snap_T_start(1) + (snap_T_par(1) - 1)*snap_T_interval(1) + sim_boundary_layers(1);
totaly = snap_T_start(2) + (snap_T_par(2) - 1)*snap_T_interval(2) + sim_boundary_layers(2);
totalz = snap_T_start(3) + (snap_T_par(3) - 1)*snap_T_interval(3) + sim_boundary_layers(3);

totalx,totaly,totalz
if allow_missing_grids && (totalx < nx || totaly<ny || totalz < nz)
    warning('Some grids are missing. Allowing it, as specified by user.')
else
    if totalx < nx
        error('Not all grids (excluding boundary layers) are saved in X direction.')
    end
    if totaly < ny
        error('Not all grids (excluding boundary layers) are saved in Y direction.')
    end
    if totalz < nz
        error('Not all grids (excluding boundary layers) are saved in Z direction.')
    end
end
snap_V_par=[nx-2*sim_boundary_layers(1) ny-2*sim_boundary_layers(2) 1];

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