%plot great circle paths between station pairs having useful data
clear all
close all;

%addpath('~/FWT/ANT//Proj/cascadia/matlab/');
MFILE_ROOT='../../mfiles';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);


tempdelaydatafile='measure_result_NEANT_withlatlon.dat';
minlon=283.8; maxlon=293.2; minlat=38.2;maxlat=47.8;

fband=[0.0067 0.01333;0.01 0.02;0.01333 0.0286;0.02 0.04;0.0286 0.0667; 0.05 .1; 0.0667 0.1333; .1 .2];
pband=flip(1./fband,2);
[nfb, nc]=size(fband);
      
fidtemp=fopen(tempdelaydatafile,'r');
tempdelaydata=textscan(fidtemp,'%s %s %f %f %f %f %f %f %s');
fclose(fidtemp);
[src,rcv,slat,slon,rlat,rlon,delay,err,fb] = tempdelaydata{1:9};
clear tempdelaydata;
nd = length(slat);

%data_var=0;
[dist,az]=distance(slat,slon,rlat,rlon,[6378.14 0.0818]);
gcplat=nan(nd,11);
gcplon=nan(nd,11);
for i=1:nd
    %clear gcpdist;
    gcpdist=0:dist(i)/10:dist(i);

%for k=1:11
    [gcplat(i,:),gcplon(i,:)]=reckon(slat(i),slon(i),gcpdist,az(i),[6378.14 0.0818]);
%end
    gcplon(i,:)=360+gcplon(i,:);
end

%%
% for i=1:nd
%     %i
%     %data_var=data_var+delay(i)^2;
%     [dist(i),az(i)]=distance(slat(i),slon(i),rlat(i),rlon(i),[6378.14 0.0818]);
%     gcpdist=0:dist(i)/10:dist(i);
%     for k=1:11
%         [gcplat(i,k),gcplon(i,k)]=reckon(slat(i),slon(i),gcpdist(k),az(i),[6378.14 0.0818]);
%     end
% end
%data_var=data_var/(nd-1);
%% get values for each frequence band
disp('--> getting values for each frequence band ...');
idf=nan(length(delay),nfb); %idf: index of frequency band, initiated as full length;
% dataf.delay=idf;
% dataf.dist=idf;
% dataf.az=idf;
% %dataf.rms=nan(nfb,1);
dataf.num=nan(nfb,1);
% dataf.src=cell(length(delay),nfb);
% dataf.rcv=cell(length(delay),nfb);
% dataf.rlon=idf;
% dataf.rlat=idf;
% dataf.slon=idf;
% dataf.slat=idf;
% dataf.gcplon=nan(length(delay),nfb,11);
% dataf.gcplat=nan(length(delay),nfb,11);

for i=1:nfb
    disp(strcat(num2str(int16(pband(i,1))),'-',num2str(int16(pband(i,2))),' s'));
    ftag=strcat('f',num2str(i));
    clear idftemp;
    idftemp=strmatch(ftag,fb);
    dataf.num(i)=length(idftemp);
    idf(1:dataf.num(i),i)=idftemp;
%     dataf.rlon(1:dataf.num(i),i)=rlon(idftemp);
%     dataf.rlat(1:dataf.num(i),i)=rlat(idftemp);
%     dataf.slon(1:dataf.num(i),i)=slon(idftemp);
%     dataf.slat(1:dataf.num(i),i)=slat(idftemp);
%     dataf.delay(1:dataf.num(i),i)=delay(idftemp);
%     dataf.dist(1:dataf.num(i),i)=dist(idftemp);
%     dataf.az(1:dataf.num(i),i)=az(idftemp);
%     %dataf.rms(i)=rms(dataf.delay(1:dataf.num(i),i));
%     dataf.gcplon(1:dataf.num(i),i,:)=gcplon(1:dataf.num(i),:);
%     dataf.gcplat(1:dataf.num(i),i,:)=gcplat(1:dataf.num(i),:);
%     
%     for j=1:dataf.num(i)
%         dataf.src{j,i}=src(idftemp(j));
%         dataf.rcv{j,i}=rcv(idftemp(j));
%     end
    
end

%%
ff = {'f1','f2','f3','f4','f5','f6','f7','f8'};
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) '};
disp('--> plotting ...');
figure('Position',[400 450 1200 650]);
%rcvloc=unique([rlon rlat],'rows');
for i=1:nfb
    %plot
    
    subplot(2,4,i) 
    hold on, box on
    clear rcvloc srcloc;
    rcvloc=nan(dataf.num(i),2);
    srcloc=nan(dataf.num(i),2);
    rcvloc(:,1)=rlon(idf(1:dataf.num(i),i));
    rcvloc(:,2)=rlat(idf(1:dataf.num(i),i));
    srcloc(:,1)=slon(idf(1:dataf.num(i),i));
    srcloc(:,2)=slat(idf(1:dataf.num(i),i));
    clear stations;
    stations=unique([rcvloc; srcloc],'rows');
%     rcvloc=unique(rcvloc,'rows');
%     srcloc=unique(srcloc,'rows');
    for j=1:dataf.num(i)
%         dlon=squeeze(dataf.gcplon(j,i,:));
%         dlat=squeeze(dataf.gcplat(j,i,:));
        dlon=gcplon(idf(j,i),:);
        dlat=gcplat(idf(j,i),:);
        plot(dlon,dlat,'k-'); 
%         disp(['source: ',dataf.src{j}]);
%         srcloc(j,:)
%         disp(['receiver: ',dataf.rcv{j}]);
%         rcvloc(j,:)
%         plot(rcvloc(j,1),rcvloc(j,2),'r^','MarkerSize',5);
%         plot(srcloc(j,1),srcloc(j,2),'bo','MarkerSize',5);
        %pause
    end
    plot(stations(:,1),stations(:,2),'r^','MarkerSize',5);
%     plot(rcvloc(:,1),rcvloc(:,2),'r^','MarkerSize',5);
    %plot(srcloc(:,1),srcloc(:,2),'bo','MarkerSize',5);
    
    xlabel('Longitude','FontSize',13)  
    ylabel('Latitude','FontSize',13)  
    axis([minlon  maxlon minlat maxlat])
    set(gca,'XTick',minlon:3:maxlon, 'YTick',minlat:2:maxlat);
    set(gca,'TickDir','out');
    title(strcat(figlabel{i},num2str(int16(pband(i,1))),'-',num2str(int16(pband(i,2))),' s:  ',...
        num2str(dataf.num(i))),'FontSize',14);
    hold off;
    %clear dlon dlat;
    drawnow;
    %clear rcvloc srcloc;
end

