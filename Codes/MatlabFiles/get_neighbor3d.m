function goodpoints=get_neighbor3d(x,y,z,mv3d,x0,y0,z0,maxdistxy,maxdistz)
%get value from neighbors in 3d grid

R0=6371; %earth's radius
% %test point
% x0=-73;
% y0=44;
% z0=15;
% maxdist=15; %defined the maximum distance of the neighbors (sphere radius).
goodpoints=[];
% distgrid=nan(size(mv3d));
n=1;
for j=1:length(y)
    disty=R0*pi*distance(y0,x0,y(j),x0)/180;
    if disty <= maxdistxy
        for i=1:length(x)
            distx=R0*pi*distance(y0,x0,y0,x(i))/180;
            if distx <= maxdistxy
                distxy=R0*pi*distance(y0,x0,y(j),x(i))/180;
                if distxy <=maxdistxy
                    for k=1:length(z)
                        if abs(z0 - z(k)) <= maxdistz
%                             dist3d=sqrt(distxy^2 + (z0 - z(k))^2);
%                             if dist3d <= maxdistxy
                                goodpoints.idx(n,1:3)=[j,i,k];
                                goodpoints.loc(n,1:3)=[y(j),x(i),z(k)];
                                goodpoints.val(n)=mv3d(j,i,k);
                                n = n + 1;
%                             end
                        end
                    end
                end
            end
        end
    end
end
%
% plot3(x0,y0,z0,'r.','markersize',25);
% hold on;
% plot3(goodpoints.loc(:,2),goodpoints.loc(:,1),goodpoints.loc(:,3),'bo','markersize',5);
% hold off;
% grid on;
% box on;

return;
end