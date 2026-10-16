%plot tradeoff curve
load tradeoff.dat
smot_list=[2 4 8 12 16];
damp_list=[2 4 8 12 16];
sel_damp=8; sel_smooth=8;
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

figure('Position',[40 20 700 600])
subplot(3,1,1), hold on, box on, axis on
cs=contour(damp,smooth,kai,'k');
clabel(cs);
plot(sel_damp,sel_smooth,'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('damping parameter','FontSize',12);
ylabel('smoothness parameter','FontSize',12);
axis([0 20 0 20]);
set(gca,'FontSize',12);
hold off;

% model norm versus variance reduction
subplot(3,1,2), hold on, box on, axis on
for i=1:nsm
    if i == isel_smooth
        plot(var(i,:),model(i,:),'b-','LineWidth',2);
    else
        plot(var(i,:),model(i,:),'k-'); 
    end
end
plot(var(isel_smooth,jsel_damp),model(isel_smooth,jsel_damp),'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('variance reduction','FontSize',12);
ylabel('model norm','FontSize',12);
axis([67 76 0.01 0.06]);
set(gca,'FontSize',12);
hold off;

% model norm versus kai square
subplot(3,1,3), hold on, box on, axis on
for i=1:nsm
    if i == isel_smooth
        plot(kai(i,:),model(i,:),'b-','LineWidth',2); 
    else
        plot(kai(i,:),model(i,:),'k-');
    end
end
plot(kai(isel_smooth,jsel_damp),model(isel_smooth,jsel_damp),'bo','MarkerSize',8,'MarkerFaceColor','b');
xlabel('x^2/N','FontSize',12);
ylabel('model norm','FontSize',12);
%axis([0 1 0 0.06]);
set(gca,'FontSize',12);
hold off;

print -depsc tradeoff.eps
