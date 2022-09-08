%script to read outputs from sim.stations and write seismograms at
% receivers in sac
addpath(genpath('/depot/xtyang/data/codes/mexcdf'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise'))
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))

PROJHOME = ['/depot/xtyang/data/projects/xtyang/craton'];
SCRATCHROOT = ['/scratch/bell/xtyang/FWANT_craton'];
ite = ['/ite_0.05deg_01'];
ite_num = [PROJHOME ite];
sim_dir = [SCRATCHROOT '/' ite '/sim.station'];
syn_dir = [ite_num '/syn.seismograms'];
cd (syn_dir);

MFILE_ROOT=[PROJHOME '/mfiles'];
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

% ----------------- station info -----------------------

stnfile = [ PROJHOME '/STinfo/craton_station_withdata.txt'];

fid = fopen(stnfile);
recv = textscan(fid,'%s%s%f%f%f'); %ntwk stnm lon lat elevation
fclose(fid);
nrecv=length(recv{1});
ntwk=recv{1}; stnm=recv{2}; 
slon=recv{3}; slat=recv{4};
xrecv=nan(nrecv); yrecv=nan(nrecv);
for n=1:nrecv
    xrecv(n)=(90-slat(n))/180*pi;
    yrecv(n)=slon(n)/180*pi;
end

fcomp = ['fz']; 
recvcomp = ['Vz']; % particle velocity
recvcompdisp = ['Uz']; % displacement
sim_input = [ ite_num '/sim.input/' ];
pnm_metric = sim_input;
% ----------------------- parameter -----------------------
% snap_002 is the velocity output on the free surface
id=2;subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
n1=1;dn=1;
n2=1500; % stepssaved every dn steps. n2 should be <= saved time steps.
layers=n1:dn:n2;
nlayers=length(layers);
parfor nstsrc=1:nrecv %
    sim_stn = [ sim_dir '/' char(ntwk(nstsrc)) '.' char(stnm(nstsrc)) ];
    sim_fcmp = [ sim_stn '/' fcomp];
    event_name=[char(ntwk(nstsrc)) '.' char(stnm(nstsrc))];
    evtid = exist(event_name,'dir');
    if ( evtid == 0 ) 
        mkdir(event_name);
    end

    OUTPUT_ROOT=[ syn_dir '/' event_name ];
    %check if already have files
    if length(dir([OUTPUT_ROOT '/*.SAC']))<2*nrecv
        disp([ '--> extracting waveforms for virtual source  ' event_name ]);
    else
        disp([event_name ' already finished. Skip!'])
        continue
    end     
    
    INPUT_ROOT=[ sim_fcmp ];
    fnm_conf=[INPUT_ROOT '/SeisFD3D.conf'];
    pnm_out =[INPUT_ROOT '/output/'];

    % -------------------- load coord --------------------------
    [snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

    %-- get coord data
    [x,y,z]=gather_coord(snapinfo,'coorddir',pnm_metric);
    nx=size(x,1);ny=size(x,2);nz=size(x,3);

    % --------------------- load data --------------------------
    Vz=[];t=[];v=[];Wz=nan(nlayers,nrecv); T=nan(nlayers,1);
    for j=1:nlayers
        ly=layers(j);
        [Vz,t]=gather_snap(snapinfo,id,ly,recvcomp,'outdir',pnm_out);
        V=interp2(squeeze(y),squeeze(x),squeeze(Vz),yrecv,xrecv,'linear');
        Wz(ly,:)=V(1:nrecv)';
        T(j)=t;
    end  % end of nlayer

    tt=[0; T];dt=tt(2)-tt(1);
    for nstrec=1:nrecv
        vzrecv=[0 Wz(:,nstrec)'];
        Szrecv=bsac(tt,vzrecv);
        [disttemp,~]=distance(slat(nstsrc),slon(nstsrc),slat(nstrec),slon(nstrec));
        Szrecv=ch(Szrecv,'STLA',slat(nstrec),'STLO',slon(nstrec),'EVLA',slat(nstsrc),'EVLO',slon(nstsrc),...
            'DIST',disttemp);
        recvname = [ char(ntwk(nstrec)) '.' char(stnm(nstrec)) ];
        wsac([event_name '/' event_name '.to.' recvname '.' fcomp '.' recvcomp '.SAC'],Szrecv);
        % for displacement
        uzrecv=cumtrapz(vzrecv)*dt;
        Szrecv=bsac(tt,uzrecv);
        Szrecv=ch(Szrecv,'STLA',slat(nstrec),'STLO',slon(nstrec),'EVLA',slat(nstsrc),'EVLO',slon(nstsrc));
        wsac([event_name '/' event_name '.to.' recvname '.' fcomp '.' recvcompdisp '.SAC'],Szrecv);
    end  % end of receiver

end  % end of source





