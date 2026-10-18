%! /bin/bash
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Love-wave (TT) version of kern_synthetic4kernel.m
%%%
%%% Copies the synthetics needed by the kernel calculation into the
%%% kernel directories (ascii: time, velocity, displacement).
%%%
%%% For each station pair (source = statref, receiver = statmes) three
%%% files are written, in the order the kernel code reads them:
%%%     <statmes>/<statref>/BHN/synthetic.vel.disp.dat   <- fT.Vx  (x = colatitude, SOUTH)
%%%     <statmes>/<statref>/BHE/synthetic.vel.disp.dat   <- fT.Vy  (y = EAST)
%%%     <statmes>/<statref>/BHZ/synthetic.vel.disp.dat   <- fT.Vz  (z = up)
%%% The names BHN/BHE/BHZ are only the slot names used by
%%% kern_measure2kernelconf; the content is the x,y,z response to a force
%%% transverse to the path at the source (written by snap2sac_sphere_fxfy.m).
%%% The kernel code rotates these to T at the receiver itself, with the
%%% rot_mat given in the *_conf files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%| set paths to use matlab routines |%%%%%%%%%%%%%%%%%%%%
MFILEROOT='../../mfiles';
path([MFILEROOT '/fun-spool'],path);
path([MFILEROOT '/saclab'],path);

%|   set directories (CHANGE for the TT project) |%%%%%%%%%%%%%%%%%%%%%
projdir = '/depot/xtyang/data/projects/xtyang/RioGrande/FWANT_ZZ/';
ite_num = 'ite_0.025deg_01';
synthdir=strcat(projdir,ite_num, '/syn.seismograms');
kerndir=strcat(projdir,ite_num, '/sim.kernel');

%|   which SAC files fill the three synthetic slots |%%%%%%%%%%%%%%%%%
fcomp   = 'fT';                       % source-side transverse force
sacvel  = {'Vx','Vy','Vz'};           % x (south), y (east), z (up)
kerncmp = {'BHN','BHE','BHZ'};        % slot directories, same order as in the conf file
% The Z slot is only read by the ORIGINAL SI_ker_sta.F90 (it reads Z whenever it
% uses rotated picks); the weight of Z in T is zero.  SI_ker_sta_fxyz.F90 never
% reads it.  With allow_missing_Z = true a missing fT.Vz (snap2sac_sphere_fxfy.m
% run with save_Z = false) is replaced by zeros instead of failing.
allow_missing_Z = true;

unix(['cat ',kerndir,'/*_conf | grep "/" | awk ',' ''{if(NF==1) print $1}''',...
    '| awk -F/ ','''{print $(NF-3), $(NF-2)}''', '| sort | uniq > .temp']);

fid=fopen('.temp','r');
texttemp=textscan(fid,'%s %s');
fclose(fid);
[starec, stasrc]=texttemp{1:2};
%correspond to shift of source time function, see SeisSource.conf in sim.station
tshift=3;
stept=0.2; % time interval of the saved synthetics (s)
nt=4000;   % number of samples read by the kernel code (SGT_TIMEDIM_SIZE)
t=(1:nt)*stept;
parfor i=1:length(starec)
    stat1=char(stasrc(i)); % source (virtual source, statref)
    stat2=char(starec(i)); % receiver (statmes)
    disp([stat1,'-',stat2])
    for c=1:3
        w=[];T=[];V=[];v=[];s=[];
        sfile=strcat(synthdir, '/', stat1, '/', stat1, '.to.', stat2, '.', fcomp, '.', sacvel{c}, '.SAC');
        if allow_missing_Z && c==3 && ~exist(sfile,'file')
            v=zeros(size(t));
        else
            w=readsac(sfile);
            T=w.B:w.DELTA:(w.E)+w.DELTA; T=T+tshift;
            V=w.DATA1;
            nn=min(numel(T),numel(V)); T=T(1:nn); V=V(1:nn);

            T=[0;T';T(end)+10]; % pad the ends before interpolation
            V=[0;V(:);0];
            v=interp1(T,V,t);
            v(isnan(v))=0;      % beyond the padded end: no signal
        end
        u=cumtrapz(v)*stept;
        s=[t' v' u'];
        soutfile=strcat(kerndir, '/', stat2, '/', stat1, '/', kerncmp{c}, '/synthetic.vel.disp.dat');
        save_variable2file(s,soutfile,'-ascii');
    end
end

unix(strcat('rm -f .temp'));
