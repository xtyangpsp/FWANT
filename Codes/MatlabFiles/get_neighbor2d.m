function goodpoints=get_neighbor2d(x,y,mv2d,x0,y0,maxdistxy)
%get value from neighbors in 2d grid
% x: x array (longitude).
% y: y array (latitude).
% mv2d: 2d matrix showing the values to extract. Size of mv2d should be
% equal to [length(x),length(y)].
% x0,y0: the point location to extract
% maxdistxy: maximum distance away from x0,y0.
%
%By Xiaotao Yang
%Email: stcyang@gmail.com

[r,c]=size(mv2d);
if r~= length(x) || c~= length(y)
    error('The size of mv2d does not match that of x and y. Required to satisfy size(mv2d)=[length(x),length(y)].')
end
% R0=6371; %earth's radius
% %test point
% x0=-73;
% y0=44;
% z0=15;
% maxdist=15; %defined the maximum distance of the neighbors (sphere radius).
goodpoints.idx=[];
goodpoints.loc=[];
goodpoints.val=[];

% n=1;
distx=deg2km(distance(y0,x0,y,x0));
clear idx00;
idx00=find(distx<=maxdistxy)';
if ~isempty(idx00)
    for j=1:length(idx00)
        disty=deg2km(distance(y0,x0,y(idx00(j)),x));
        clear idx0;
        idx0=find(disty<=maxdistxy)';
        if ~isempty(idx0)
            idx0=reshape(idx0,length(idx0),1);
            goodpoints.idx=[goodpoints.idx;idx0,ones(length(idx0),1)*idx00(j)];
            goodpoints.loc=[goodpoints.loc;reshape(x(idx0),length(idx0),1),ones(length(idx0),1)*y(idx00(j))];
%             j
%             idx00(j)
%             idx0
            goodpoints.val=[goodpoints.val;reshape(mv2d(idx0,idx00(j)),length(idx0),1)];
        end
    end
end
    %
% if ~isempty(goodpoints.idx)
%     plot(x0,y0,'r.','markersize',25);
%     hold on;
%     plot(goodpoints.loc(:,2),goodpoints.loc(:,1),'bo','markersize',5);
%     hold off;
%     grid on;
%     box on;
% end

return;
end