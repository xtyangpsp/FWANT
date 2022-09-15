function plt_greatcirclepath(tempdelaydatafile,fband,subplotpar)
%plot great circle paths between station pairs having useful data

pband=flip(1./fband,2);
[nfb, ~]=size(fband);
      
fidtemp=fopen(tempdelaydatafile,'r');
tempdelaydata=textscan(fidtemp,'%s %s %f %f %f %f %f %f %s');
fclose(fidtemp);
[~,~,slat,slon,rlat,rlon,delay,~,fb] = tempdelaydata{1:9};
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

dataf.num=nan(nfb,1);

for i=1:nfb
    disp(strcat(num2str(int16(pband(i,1))),'-',num2str(int16(pband(i,2))),' s'));
    ftag=strcat('f',num2str(i));
    clear idftemp;
    idftemp=strmatch(ftag,fb);
    dataf.num(i)=length(idftemp);
    idf(1:dataf.num(i),i)=idftemp;
    
end

%%
% ff = {'f1','f2','f3','f4','f5','f6','f7','f8'};
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) '};
disp('--> plotting ...');
figure('Position',[400 450 1200 650]);
for i=1:nfb
    %plot
    
    %if ~isempty(subplotpar)
    subplot(subplotpar(1),subplotpar(2),i) 
    %end
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
    for j=1:dataf.num(i)
        dlon=gcplon(idf(j,i),:);
        dlat=gcplat(idf(j,i),:);
        plot(dlon,dlat,'k-'); 
    end
    plot(stations(:,1),stations(:,2),'r^','MarkerSize',5);
    xlabel('Longitude','FontSize',13)  
    ylabel('Latitude','FontSize',13)  
%     axis([minlon  maxlon minlat maxlat])
%     set(gca,'XTick',minlon:3:maxlon, 'YTick',minlat:2:maxlat);
    set(gca,'TickDir','out');
    title(strcat(figlabel{i},num2str(int16(pband(i,1))),'-',num2str(int16(pband(i,2))),' s:  ',...
        num2str(dataf.num(i))),'FontSize',14);
    hold off;
    %clear dlon dlat;
    drawnow;
    %clear rcvloc srcloc;
end

end
