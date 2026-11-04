%plot tradeoff curve
close all
load tradeoff.dat
smot_list=[2 4 6 8 10 12 16]; %[0.5 1 2 3 4 5 6 8 10];
damp_list=smot_list; %[2 4 8 12 16];
sel_damp=10; sel_smooth=8;
nsm=length(smot_list);  ndamp=length(damp_list);
isel_smooth=find(smot_list == sel_smooth); %corresponds to 
jsel_damp=find(damp_list == sel_damp);

nk=nsm*ndamp;
ii=0;
for j=1:ndamp
    for i=1:nsm
        ii=ii+1;
        kai(i,j)=tradeoff(ii,3);
        damp(i,j)=tradeoff(ii,1);
        smooth(i,j)=tradeoff(ii,2);
        var(i,j)=tradeoff(ii,4);
        model(i,j)=tradeoff(ii,5);
    end
end

figure('Position',[40 20 500 700])
subplot(3,1,1), hold on, box on, axis on
cs=contour(damp,smooth,kai,'linewidth',1);
clabel(cs);
plot(sel_damp,sel_smooth,'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('Damping','FontSize',12);
ylabel('Smoothness','FontSize',12);
title('a','FontSize',14)
% axis([1 17 1 17]);
colormap("turbo")
cbar=colorbar('location','eastoutside');
cbar.Label.String='\chi^2/N';
set(gca,'FontSize',12,'TitleHorizontalAlignment','left');
hold off;

% model norm versus variance reduction
subplot(3,1,2), hold on, box on, axis on
for i=1:nsm
    if i == isel_smooth
        plot(var(i,:),model(i,:),'b-','LineWidth',2);
    else
        plot(var(i,:),model(i,:),'k-','LineWidth',1); 
    end
end
plot(var(isel_smooth,jsel_damp),model(isel_smooth,jsel_damp),'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('Variance reduction (%)','FontSize',12);
ylabel('Model norm','FontSize',12);
title('b','FontSize',14)
% axis([67 76 0.01 0.06]);
set(gca,'FontSize',12,'TitleHorizontalAlignment','left');
hold off;

% model norm versus kai square
subplot(3,1,3), hold on, box on, axis on
for i=1:nsm
    if i == isel_smooth
        plot(kai(i,:),model(i,:),'b-','LineWidth',2); 
    else
        plot(kai(i,:),model(i,:),'k-','LineWidth',1);
    end
end
plot(kai(isel_smooth,jsel_damp),model(isel_smooth,jsel_damp),'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('\chi^2/N','FontSize',12);
ylabel('Model norm','FontSize',12);
title('c','FontSize',14)
%axis([0 1 0 0.06]);
set(gca,'FontSize',12,'TitleHorizontalAlignment','left');
hold off;

set(gcf,'renderer','Painters','PaperPositionMode','auto'); 
print -dpdf -r300 -bestfit -vector tradeoff.pdf

%% L-curve only, with labeled smooth and damping parameters
figure('Position',[40 20 500 350]);

% model norm versus variance reduction
hold on, box on, axis on
for i=1:nsm
    if i == isel_smooth
        plot(var(i,:),model(i,:),'bo-','LineWidth',2,'MarkerSize',10);
    else
        plot(var(i,:),model(i,:),'ko-','LineWidth',1,'MarkerSize',10); 
    end
    text(var(i,1)-0.04,0.002+model(i,1),['s',num2str(smot_list(i))],'fontsize',12);
end

for j=1:ndamp
    if j == jsel_damp
        plot(var(:,j),model(:,j),'bo:','LineWidth',2,'MarkerSize',10);
    else
        plot(var(:,j),model(:,j),'ko:','LineWidth',1,'MarkerSize',10); 
    end
    text(var(1,j)+0.02,model(1,j)-0.001,['d',num2str(damp_list(j))],'fontsize',12);
end

h=plot(var(isel_smooth,jsel_damp),model(isel_smooth,jsel_damp),'bo','MarkerSize',15,'MarkerFaceColor','b',...
    'DisplayName','Selected');
xlabel('Variance reduction (%)','FontSize',12);
ylabel('Model norm','FontSize',12);
% ylim([0.01,0.05])
% title('b','FontSize',14)
% axis([67 76 0.01 0.06]);
set(gca,'FontSize',12,'TickDir','out','TitleHorizontalAlignment','left');
hold off;
legend(h);
print -dpdf -r300 -bestfit -vector tradeoff_Lcurve.pdf
