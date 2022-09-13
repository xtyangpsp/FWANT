% close all; %clear all;

tempdelaydatafilelist={'ite_0.05deg_03_measure_result_Alaska_flat9.0.dat_latlon',...
    'ite_0.05deg_04_measure_result_Alaska_flat9.0.dat_latlon'};
numinput=numel(tempdelaydatafilelist); 
%%%% define parameters %%%%
max_dT = 9; % maximum delay in second (see mesure_phase_delay.csh)

fband=[0.005 0.01;0.0067 0.01333;0.01 0.02;0.01333 0.0286;0.02 0.04;0.0286 0.0667; 0.04 .1; 0.0667 .1333];
%fband=[0.0067 0.01333;0.01 0.02;0.01333 0.0286;0.02 0.04;0.0286 0.0667];
pband=flip(1./fband,2);
[nfb, nc]=size(fband);
      
%% 
figure('Position',[400 400 1200 500]);
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) '};
colorlist={'k','r'};
for k=1:numinput
    tempdelaydatafile=tempdelaydatafilelist{k};
    fidtemp=fopen(tempdelaydatafile,'r');
    tempdelaydata=textscan(fidtemp,'%s %s %f %f %f %f %f %f %s %*f %*f');
    fclose(fidtemp);
    [src,rcv,slat,slon,rlat,rlon,delay,err,fb] = tempdelaydata{1:9};
    clear tempdelaydata;
    nd = length(slat);

    %% get values for each frequence band
    idf=nan(length(delay),nfb); %idf: index of frequency band, initiated as full length;
    dataf.delay=idf;
    %dataf.dist=idf;
    % dataf.rms=nan(nfb,1);
    dataf.num=nan(nfb,1);
    % dataf.src=cell(length(delay),nfb);
    % dataf.rcv=cell(length(delay),nfb);

    for i=1:nfb
        disp(strcat(num2str(int16(pband(i,1))),'-',num2str(int16(pband(i,2))),' s'));
        ftag=strcat('f',num2str(i));
        idftemp=strmatch(ftag,fb);
        dataf.num(i)=length(idftemp);
        idf(1:dataf.num(i),i)=idftemp;

        dataf.delay(1:dataf.num(i),i)=delay(idftemp);
        %dataf.dist(1:dataf.num(i),i)=dist(idftemp);
    %     dataf.rms(i)=rms(dataf.delay(1:dataf.num(i),i));

    %     for j=1:dataf.num(i)
    %         dataf.src{j,i}=src(idftemp(j));
    %         dataf.rcv{j,i}=rcv(idftemp(j));
    %     end
    end

    %%
    disp('Plotting ...');
    for i=1:nfb
        %plot

        subplot(2,4,i), hold on, box on

        [n,xout]=hist(dataf.delay(:,i),round((max(dataf.delay(:,i))-min(dataf.delay(:,i)))));
        bar(xout, n, 'EdgeColor',colorlist{k},'FaceColor','none','LineWidth',1);
        xlim([-max_dT max_dT]);
        c = length(xout);
        w = xout(c)-xout(1);
        t = linspace(xout(1)-w/c,xout(end)+w/c,c+1);
        dt = diff(t);
        Fvals = cumsum([0,n.*dt]);
        F = spline(t, [0, Fvals, 0]);
        DF = fnder(F);  % computes its first derivative
        fnplt(DF, colorlist{k}, 2);
        
        maxfnval(k,i)=max(fnval(DF,xout(1):.1:xout(end)));
        if k==numinput
            plot([0 0],[0 max(maxfnval(:,i))+10],'k--','LineWidth',.5)
            ylim([0 max(maxfnval(:,i))+10]);
            %text(max_dT-1,100,['RMS = ' num2str(dataf.rms(i),3)]);
            xlabel('Delay, s','FontSize',12);
            set(gca,'XTick',-max_dT:1:max_dT,'XTickLabel',-max_dT:1:max_dT);
            title(strcat(figlabel{i}, num2str(int16(pband(i,1))),'-',num2str(int16(pband(i,2))),' s:  ',...
                num2str(dataf.num(i))),'FontSize',14);
        end
        %pause;
    end
    drawnow;
end 
hold off;
legend('ite-03','ite-03','ite-04','ite-04');
% set(gcf,'PaperPositionMode','auto');   
% print -depsc Histogram_dt.eps