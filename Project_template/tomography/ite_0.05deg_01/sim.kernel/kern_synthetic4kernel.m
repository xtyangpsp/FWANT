%! /bin/bash
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to copy synthetics with time windows around Rayleigh phase
%%% to kernel calculation directories
%%% files are written out in ascii format
%%% Modified from kern_synthetic4kernel.sh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%| set paths to use matlab routines |%%%%%%%%%%%%%%%%%%%%
MFILEROOT='../../mfiles';
path([MFILEROOT '/fun-spool'],path);
path([MFILEROOT '/saclab'],path);

%|   set directories              |%%%%%%%%%%%%%%%%%%%%%
ite_num = 'ite_0.05deg_01';
synthdir=strcat('/depot/xtyang/data/projects/xtyang/craton/',ite_num, '/syn.seismograms');
%kerndir=strcat('/scratch/bell/xtyang/FWANT_craton/',ite_num, '/sim.kernel');
kerndir=strcat('/depot/xtyang/data/projects/xtyang/craton/',ite_num, '/sim.kernel');

unix(['cat ',kerndir,'/*_conf | grep "/" | awk ',' ''{if(NF==1) print $1}''',...
    '| awk -F/ ','''{print $(NF-3), $(NF-2)}''', '| sort | uniq > .temp']);

fid=fopen('.temp','r');
texttemp=textscan(fid,'%s %s');
fclose(fid);
[starec, stasrc]=texttemp{1:2};
%correspond to shift of source time function, see SeisSource.conf in sim.station
tshift=6;
stept=0.8; % 0.2 s per time step, saved every 4 steps
nt=1500; %5650 steps total, save every 5 steps
t=(1:nt)*stept;
parfor i=1:length(starec)
    w=[];T=[];V=[];v=[];s=[];
	stat1=char(stasrc(i)); % 
	stat2=char(starec(i)); % 
	%comp=char('BHZ');
        disp([stat1,'-',stat2])	
	sfile=strcat(synthdir, '/', stat1, '/', stat1, '.to.', stat2, '.fz.Vz.SAC');
	w=readsac(sfile); %T=w(:,1); V=w(:,2); T=T+tshift; 
	T=w.B:w.DELTA:(w.E)+w.DELTA; T=T+tshift;
	V=w.DATA1;
 
	T=[0;T';1210]; %synthetics is 1200; pad the ends before interpolation
	V=[0;V;0];
    v=interp1(T,V,t); 
	u=cumtrapz(v)*stept;
	s=[t' v' u'];
	soutfile=strcat(kerndir, '/', stat2, '/', stat1, '/BHZ/synthetic.vel.disp.dat');
	save_variable2file(s,soutfile,'-ascii');
	%save(soutfile, 's', '-ascii');
%     clear w T V v s
end

unix(strcat('rm -f .temp'));
