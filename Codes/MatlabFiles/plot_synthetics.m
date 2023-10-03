addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))
R0=6371;

% define filter parameters. Currently not used.
T_min = 20; T_max = 50; %period min and max.
fband=[1/T_max, 1/T_min]; % frequency band

% km/s, minimum and maximum group velocities
cmin=2.5; cmax=4.5;
signal_window_extend=50;
npoints=1501;

%%%%%%%%%%%%%% END OF - Setting up global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sta = 'CN.SADO';
unix(['ls -f ' sta '/*.fz.Uz.SAC > filelist3']);
synfile1=textread('filelist3' ,'%s');
nsynfile=length(synfile1);

synfbdata=zeros(npoints,nsynfile);
evla=nan(nsynfile,1);
evlo=nan(nsynfile,1);
stla=nan(nsynfile,1);
stlo=nan(nsynfile,1);
distall=nan(nsynfile,1);
for i=1:nsynfile
    filename1=char(synfile1(i));
    %clear dd1 tt
    dd3=readsac(filename1);
    dtsyn=dd3.DELTA;
    ttsyn=dd3.B:dtsyn:dd3.E;
    nptsyn=dd3.NPTS;
    Feff=1/dtsyn/2; % in increment of 2

    evla(i)=dd3.EVLA;evlo(i)=dd3.EVLO;
    stla(i)=dd3.STLA;stlo(i)=dd3.STLO;
   
    [distall(i),~]=distance(evla(i),evlo(i),stla(i),stlo(i));
    distall(i)=round(R0*distall(i)*pi/180.0);
    synfbdata(:,i)=dd3.DATA1;
end
%%
figure('Position',[300, -200, 800, 600]); hold on;

scale=100;

for i=1:nsynfile
    tmin=distall(i)/cmax-T_max;tmax=signal_window_extend + distall(i)/cmin+T_max;
    itmin=round(tmin/dtsyn); itmax=round(tmax/dtsyn);
    if itmin < 1; itmin = 1; end
    if itmax > nptsyn; itmax = nptsyn; end
    plot(ttsyn(itmin:itmax),scale*synfbdata(itmin:itmax,i)+distall(i),'k');
end
hold off;
xlim([0 400]);
ylim([150 800]);
box on;
axis on;
xlabel('lag time (s)')
ylabel('inter-station distance (km)')
set(gca,'tickdir','out')
title(['synthetics: ',sta],'FontSize',14);

saveas(gca,['synthetics_',sta,'.png'])