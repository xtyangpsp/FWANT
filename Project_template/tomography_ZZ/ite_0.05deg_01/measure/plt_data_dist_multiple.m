function plt_data_dist_multiple(delayfilelist,fband,subplotpar,maxdt,maxy,cutoffmethod)
%USAGE: plt_data_dist_yang_multiple(stationinfofile,rawdelaydatafilelist,fband,subplotpar)
%modified from plt_data_dist.m from Haiying Gao
%Main modifications:
%   1. use loop in extracting data and plotting each frequency band
%   2. avoid hard-coded frequency IDs.
%   3. read frequency band authomatically from specified array/vector.
%   4. improved cutoff methods. added two more methods. and changed from
%   directly removing data to only writing removing list into a text file
%   for later use. this makes the cutoff procedure more migratable,
%   which means copying the shell script (seperate file: exclude_measurements.sh) and the list
%   from a machine to another machine and running the script will actually 
%   apply the cutoffs (data selection procedure).
%   5. loop on multiple data sources.
% cutoffmethod='draw'; 
            %(1) flat (single threshold for all period bands), uses the 
            % value specified for 'delaycutoff'. 
            %(2) norm: cut by confidence level of the associated normal fit.
            %NEED to give the number of sigmas (standard deviation) below.
            %(3) nfalc: norm fit after linear correction. This method
            %calculate the mean and sigma after linear correction. The raw
            %data will not be corrected. Cutoff will be applied with the
            %linear fit corrections.
            %(4) draw: manually draw two bounding lines (lower and upper
            %limits) with cursor.
% Xiaotao Yang @ UMASS
% 12/13/2016
% 1/24/2017
%
%clear all; close all;
%%
%addpath('~/FWT/ANT/Proj/cascadia/matlab/');
%%%%%%%%%%%%%%%%%%%GLOBAL PARAMETERS%%%%%%%%%%%%%%%
%stationinfofile='../00masterdatafiles/Alaska_stations_withdata20170831.txt';
%rawdelaydatafilelist={'ite_0.05deg_01_measure_result_Alaska_simpleaverage.dat','ite_0.05deg_06_measure_result_Alaska.dat'};
numinput=numel(delayfilelist);
%Alaska fband
% pband=[100, 200; 75,150;50,100;35,75;25,50;15,35;10,25;7.5,15];fband=flip(1./pband,2);
% fband=[0.005 0.01;0.0067 0.01333;0.01 0.02;0.01333 0.0286;0.02 0.04;0.0286 0.0667; 0.04 .1; 0.0667 .1333];

pband=flip(1./fband,2);
[nfb, ~]=size(fband);
cmax=5.5;
max_dT = maxy; % maximum delay in second (see mesure_phase_delay.csh)
            %this is actually used as yaxis limit.
delaycutoff = maxdt;%7.5; %this is used to plot the bounding lines in each figure.

if strcmp(cutoffmethod,'no')
    cleanbydelaycutoff = 0; %if 1, the program removes measurement files following the delaycutoff.
            %I think it is better to apply cutoff on the assembled file
            %instead of the original measurement files.
    plotcutoff=0;
else
    cleanbydelaycutoff = 1; 
    plotcutoff=1;
end
% cutoffmethod='draw'; 
            %(1) flat (single threshold for all period bands), uses the 
            % value specified for 'delaycutoff'. 
            %(2) norm: cut by confidence level of the associated normal fit.
            %NEED to give the number of sigmas (standard deviation) below.
            %(3) nfalc: norm fit after linear correction. This method
            %calculate the mean and sigma after linear correction. The raw
            %data will not be corrected. Cutoff will be applied with the
            %linear fit corrections.
            %(4) draw: manually draw two bounding lines (lower and upper
            %limits) with cursor.
nsigma=1.96;  %confidence interval: [mu-nsigma*sigma, mu+nsigma*sigma];   
%%%%%%%%%%%%%%%%%%%END OF GLOBAL PARAMETERS%%%%%%%%%%%%%%%
%%
figure('Position',[400 400 1300 650]);
% colorsymbollist={'k.','r.'};
% cmap=colormap(jet(numinput));
cmap=[0.4 0.4 0.4;1 0 0];
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) ','(i) ','(j) ','(k) ','(l) '};
for k=1:numinput
    tempdelayfile=delayfilelist{k};

    fidtemp=fopen(tempdelayfile,'r');
    tempdelaydata=textscan(fidtemp,'%s %s %f %f %f %f %f %*f %s %*f %*f');
    fclose(fidtemp);
    [src,rcv,slat,slon,rlat,rlon,delay,fb] = tempdelaydata{1:8};
    clear tempdelaydata;
    %std(delay);

%     nd = length(slat);
    [dist,~]=distance(slat,slon,rlat,rlon,[6378.14 0.0818]);
%     data_var=0;
%     for i=1:nd
%         data_var=data_var+delay(i)^2;
%         %dist(i)=geo2dist(slat(i),slon(i),rlat(i),rlon(i));
%     end
%     data_var=data_var/(nd-1);

    % get values for each frequence band
    idf=nan(length(delay),nfb); %idf: index of frequency band, initiated as full length;
    dataf.delay=idf;
    dataf.dist=idf;
    dataf.rms=nan(nfb,1);
    dataf.num=nan(nfb,1);
    dataf.src=cell(length(delay),nfb);
    dataf.rcv=cell(length(delay),nfb);

    for i=1:nfb
        disp(strcat(num2str(int16(pband(i,1))),'-',num2str(int16(pband(i,2))),' s'));
        ftag=strcat('f',num2str(i));
        idftemp=find(strcmp(ftag,fb)==1);
        
        dataf.num(i)=length(idftemp);
        idf(1:dataf.num(i),i)=idftemp;

        dataf.delay(1:dataf.num(i),i)=delay(idftemp);
        dataf.dist(1:dataf.num(i),i)=dist(idftemp);
        dataf.rms(i)=rms(dataf.delay(1:dataf.num(i),i));

        for j=1:dataf.num(i)
            dataf.src{j,i}=src(idftemp(j));
            dataf.rcv{j,i}=rcv(idftemp(j));
        end
    end

    % get interval used for cutoff later
    if cleanbydelaycutoff
        if strcmp(cutoffmethod,'norm')
            cinv=nan(2,nfb);
            for i=1:nfb
                [mu0,sigma0]=normfit(dataf.delay(1:dataf.num(i),i));
                cinv(:,i)=[mu0-nsigma*sigma0;mu0+nsigma*sigma0];
            end
        elseif strcmp(cutoffmethod,'nfalc') 
            cinv=nan(length(delay),3,nfb);%3 values: 1, lower limit, 2, fit (prediction or correction), 
                                            % 3, upper limit
            delaytemp=nan(length(delay),nfb);
            for i=1:nfb
                p=polyfit(dataf.dist(1:dataf.num(i),i),dataf.delay(1:dataf.num(i),i),1);
                cinv(1:dataf.num(i),2,i)=p(1).*dataf.dist(1:dataf.num(i),i) + p(2); 
                            %correction or predicted value based on linear fit
                delaytemp(1:dataf.num(i),i)=dataf.delay(1:dataf.num(i),i) - cinv(1:dataf.num(i),2,i);
                            %correct delay data for normal fit next.
                [mu0,sigma0]=normfit(delaytemp(1:dataf.num(i),i));
                cinv(1:dataf.num(i),1,i)=cinv(1:dataf.num(i),2,i)+(mu0-nsigma*sigma0);
                cinv(1:dataf.num(i),3,i)=cinv(1:dataf.num(i),2,i)+(mu0+nsigma*sigma0);
            end
        elseif strcmp(cutoffmethod,'draw')
            disp('--> for <draw> method, proceed to each plot.')
        elseif strcmp(cutoffmethod,'dv')
            disp('--> for dv max cutoff, proceed to each plot')
        end
    end
    %
    disp('Plotting ...');
    if strcmp(cutoffmethod,'draw')
        cinv=nan(length(delay),2,nfb);
        pl_all=cell(nfb,1); %save the slope and intercept of the lines.
        pu_all=cell(nfb,2);
    end
    for i=1:nfb
        %plot
        subplot(subplotpar(1),subplotpar(2),i)
        hold on, box on
        mindist=4.5*max(pband(i,:));
        maxdist=max(dataf.dist(1:dataf.num(i),i));
        plot(dataf.dist(:,i),dataf.delay(:,i),'k.','color',cmap(k,:),'MarkerSize',8);
        xlim([mindist maxdist+25]);
        ylim([-max_dT max_dT]);
    %     plot([dist(idf) dist(idf)],[delay(idf)+err(idf) delay(idf)-err(idf)],'k-'); 
    %    errorbar(dist(idf),delay(idf),err(idf),'k.','MarkerSize',7);
        if plotcutoff
            if strcmp(cutoffmethod,'flat') 
                plot([mindist maxdist],[-delaycutoff -delaycutoff],'b-');
                plot([mindist maxdist],[delaycutoff delaycutoff],'b-');
            elseif strcmp(cutoffmethod,'norm') 
                plot([mindist maxdist],[cinv(1,i) cinv(1,i)],'b-');
                plot([mindist maxdist],[cinv(2,i) cinv(2,i)],'b-');
            elseif strcmp(cutoffmethod,'nfalc') 
                plot(dataf.dist(1:dataf.num(i),i),cinv(1:dataf.num(i),2,i),'r-');
                plot(dataf.dist(1:dataf.num(i),i),cinv(1:dataf.num(i),1,i),'b-');
                plot(dataf.dist(1:dataf.num(i),i),cinv(1:dataf.num(i),3,i),'b-');
            elseif strcmp(cutoffmethod,'draw')
                disp('choose two points for lower bound line:')
                [x,y]=ginput(1);
                plot(x,y,'bo')
                xl(1)=x;yl(1)=y;
                [x,y]=ginput(1);
                plot(x,y,'bo')
                xl(2)=x;yl(2)=y;
                pl=polyfit(xl,yl,1);
                yl_pred=polyval(pl,dataf.dist(1:dataf.num(i),i));
                pl_all{i}=pl;
                hold on;
                plot(dataf.dist(1:dataf.num(i),i),yl_pred,'b');
            
                disp('choose two points for upper bound line:')
                [x,y]=ginput(1);
                plot(x,y,'bo')
                xu(1)=x;yu(1)=y;
                [x,y]=ginput(1);
                plot(x,y,'bo')
                xu(2)=x;yu(2)=y;
                plot(xu,yu,'bo')
                pu=polyfit(xu,yu,1);
                pu_all{i}=pu;
                yu_pred=polyval(pu,dataf.dist(1:dataf.num(i),i));
            
                plot(dataf.dist(1:dataf.num(i),i),yu_pred,'b');
                cinv(1:dataf.num(i),1,i)=yl_pred;
                cinv(1:dataf.num(i),2,i)=yu_pred;
            elseif strcmp(cutoffmethod,'dv')
                dv_max=input(['--> Type dv max (0-1)  for f' num2str(i) ': ']);

                pl=[-0.95*dv_max/cmax, 0];
                yl_pred=polyval(pl,dataf.dist(1:dataf.num(i),i));
                pl_all{i}=pl;
                hold on;
                plot(dataf.dist(1:dataf.num(i),i),yl_pred,'b');
            
                pu=[0.95*dv_max/cmax, 0];
                pu_all{i}=pu;
                yu_pred=polyval(pu,dataf.dist(1:dataf.num(i),i));
            
                plot(dataf.dist(1:dataf.num(i),i),yu_pred,'b');
                cinv(1:dataf.num(i),1,i)=yl_pred;
                cinv(1:dataf.num(i),2,i)=yu_pred;
            end
        end
        
        if k==numinput
            
            plot([mindist maxdist],[0 0],'k-');
            xlabel('Distance (km)','FontSize',12);
            ylabel('Delay (s)','FontSize',12);
            title([figlabel{i},'  ',num2str(pband(i,1),4),'-',num2str(pband(i,2),4),' s'],'FontSize',14);
            
            if length(delayfilelist)>1
                legend('ref.','ite-1')
            end
            axis on;
            grid on;
            set(gca,'YTick',-maxy:4:maxy,'fontsize',14)
        end
        %hold off;
        drawnow;
    end
    set(gcf,'PaperPositionMode','auto');
    eval(['print -dpng -painters ' tempdelayfile '_plot.png']);
    % cutoff procedures
    if cleanbydelaycutoff==1
        if strcmp(cutoffmethod,'flat') 
            cleanupscript=strcat(tempdelayfile, '_droplist_flat',num2str(delaycutoff),'.txt'); 
            disp(['Generating clean-up list: [ ' cleanupscript ' ] ...']);
            %open clean-up file
            fidclean=fopen(cleanupscript,'w');
            % cutoff by single threshold
            cc = 0;
            for i=1:nfb
              for j=1:dataf.num(i)
                 if abs(dataf.delay(j,i)) > delaycutoff
                    cc=cc+1;
                    fprintf(fidclean,'%s_%s  f%s\n',cell2mat(dataf.src{j,i}),...
                                cell2mat(dataf.rcv{j,i}),num2str(i));
                 end
              end
            end
            fclose(fidclean);
            
            %write cleanup procedures and parameters
            fidlog=fopen([tempdelayfile '_cleanup_log.txt'],'w');
            fprintf(fidlog,'original file: %s\n',tempdelayfile);
            fprintf(fidlog,'method: %s\n',cutoffmethod);
            fprintf(fidlog,'parameter: %f\n',delaycutoff);
            fclose(fidlog);
        elseif strcmp(cutoffmethod,'norm')
            cleanupscript=strcat(tempdelayfile,'_droplist_norm',num2str(nsigma),'.txt'); 
            disp(['Generating clean-up list: [ ' cleanupscript ' ] ...']);
            %open clean-up file
            fidclean=fopen(cleanupscript,'w');
            % cutoff by confidence interval based on statistical 
            % analysis of the delay for each frequency group.
            cc=0;
            for i=1:nfb
              for j=1:dataf.num(i)
                 if dataf.delay(j,i) > cinv(2,i) || dataf.delay(j,i) < cinv(1,i)
                    cc=cc+1;
                    fprintf(fidclean,'%s_%s  f%s\n',cell2mat(dataf.src{j,i}),...
                                cell2mat(dataf.rcv{j,i}),num2str(i));
                 end
              end
            end
            fclose(fidclean);

            %write cleanup procedures and parameters
            fidlog=fopen([tempdelayfile '_cleanup_log.txt'],'w');
            fprintf(fidlog,'original file: %s\n',tempdelayfile);
            fprintf(fidlog,'method: %s\n',cutoffmethod);
            fprintf(fidlog,'parameter: %f\n',nsigma);
            fclose(fidlog);

        elseif strcmp(cutoffmethod,'nfalc') 
            cleanupscript=strcat(tempdelayfile,'_droplist_nfalc',num2str(nsigma),'.txt'); 
            disp(['Generating clean-up list: [ ' cleanupscript ' ] ...']);
            %open clean-up file
            fidclean=fopen(cleanupscript,'w');
            % cutoff by confidence interval based on statistical 
            % analysis of the delay for each frequency group.
            cc=0;
            for i=1:nfb
              for j=1:dataf.num(i)
                 if dataf.delay(j,i) > cinv(j,3,i) || dataf.delay(j,i) < cinv(j,1,i)
                    cc=cc+1;
                    fprintf(fidclean,'%s_%s  f%s\n',cell2mat(dataf.src{j,i}),...
                                cell2mat(dataf.rcv{j,i}),num2str(i));
                 end
              end
            end
            fclose(fidclean);

            %write cleanup procedures and parameters
            fidlog=fopen([tempdelayfile '_cleanup_log.txt'],'w');
            fprintf(fidlog,'original file: %s\n',tempdelayfile);
            fprintf(fidlog,'method: %s\n',cutoffmethod);
            fprintf(fidlog,'parameter: %f\n',nsigma);
            fclose(fidlog);
         elseif strcmp(cutoffmethod,'draw') || strcmp(cutoffmethod,'dv')
             if strcmp(cutoffmethod,'draw')
                cleanupscript=strcat(tempdelayfile,'_droplist_draw.txt');
             elseif strcmp(cutoffmethod,'dv')
                 cleanupscript=strcat(tempdelayfile,'_droplist_dv.txt');
             end
            disp(['Generating clean-up list: [ ' cleanupscript ' ] ...']);
            %open clean-up file
            fidclean=fopen(cleanupscript,'w');
            % cutoff by confidence interval based on statistical 
            % analysis of the delay for each frequency group.
            cc=0;
            if strcmp(cutoffmethod,'dv')
                fprintf(fidclean,'cmax: %g',cmax);
            end
            for i=1:nfb
              for j=1:dataf.num(i)
                 if dataf.delay(j,i) > cinv(j,2,i) || dataf.delay(j,i) < cinv(j,1,i)
                    cc=cc+1;
                    fprintf(fidclean,'%s_%s  f%s\n',cell2mat(dataf.src{j,i}),...
                                cell2mat(dataf.rcv{j,i}),num2str(i));
                 end
              end
            end
            fclose(fidclean);

            %write cleanup procedures and parameters
            fidlog=fopen([tempdelayfile '_cleanup_log.txt'],'w');
            fprintf(fidlog,'original file: %s\n',tempdelayfile);
            fprintf(fidlog,'method: %s\n',cutoffmethod);
            for i=1:nfb
                fprintf(fidlog,'parameter f%d: lower line [%f %f], upper line [%f %f]\n',...
                    i,pl_all{i}(1),pl_all{i}(2),pu_all{i}(1),pu_all{i}(2));
            end
            fclose(fidlog);
        else
            error(['*** Unknown value of cutoffmethod: ',cutoffmethod]);
        end

        disp(['---> Need to remove [ ',num2str(cc),' ] data points']);

        

    end
end
hold off;

end