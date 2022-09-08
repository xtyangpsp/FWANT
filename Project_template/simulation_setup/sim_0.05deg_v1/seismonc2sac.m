clear all

set_mfiles_path
%MFILEROOT='/net/fs01/data/tibet/code/mfiles';
%path([MFILEROOT '/fun-spool'],path);
%path([MFILEROOT '/saclab'],path);

%fnm_nc='2001.300.05.35.42.synthetic.interp2.linear.nc';
fnm_nc='2003.230.09.03.02.synthetic.interp2.linear.nc';

pnm_sac = './output_sac';
%event_name='2001.300.05.35.42';
event_name='2003.230.09.03.02';
if ~ isdir(pnm_sac)
   mkdir(pnm_sac);
end

ndim=nc_getdiminfo(fnm_nc,'number_of_station'); NSTAT=ndim.Length;

for n=1:NSTAT
    t=nc_varget(fnm_nc,'time');
    Vx=nc_varget(fnm_nc,'Vx',[n-1,0],[1,-1]);
    Vy=nc_varget(fnm_nc,'Vy',[n-1,0],[1,-1]);
    Vz=nc_varget(fnm_nc,'Vz',[n-1,0],[1,-1]);
    lat=nc_varget(fnm_nc,'latitude',[n-1],[1]);
    lon=nc_varget(fnm_nc,'longitude',[n-1],[1]);
    snm=nc_varget(fnm_nc,'station_name',[n-1,0],[1,-1]);

    Vx=-Vx/10^7;
    Vy=Vy/10^7;
    Vz=Vz/10^7;


    Sx=bsac(t,Vx); Sy=bsac(t,Vy); Sz=bsac(t,Vz);
    ch(Sx,'STLA',lat,'STLO',lon);
    ch(Sy,'STLA',lat,'STLO',lon);
    ch(Sz,'STLA',lat,'STLO',lon);

    wsac([pnm_sac '/' event_name '.' strtrim(snm) '.FHN.SAC'],Sx);
    wsac([pnm_sac '/' event_name '.' strtrim(snm) '.FHE.SAC'],Sy);
    wsac([pnm_sac '/' event_name '.' strtrim(snm) '.FHZ.SAC'],Sz);
end
