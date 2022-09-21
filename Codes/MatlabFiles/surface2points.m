function [z,xin,yin]=surface2points(xyzfile,x,y,dx,dy)
%Description: Extract point z values for given points [x, y] from surface defined by
%xyzfile.
%USAGE: z=surface2points(xyzfile,x,y,dx,dy)
%xyzfile: xyzfile used to define the surface;
%x,y: location of points;
%dx,dy: grid size in gridding the surface before extracting the points.
%
%Xiaotao Yang @ Indiana University   10/12/2015
%
% modification
%   1. Sep 6, 2017. Use boundary() to get the outline enclosing the original points. For
%   points outside the boundary, set them to NaN to avoid interpolating
%   non-existing values.
%   2. Sep 7, 2017. Improve efficiency by running griddata() within the
%   area enclosing the (x,y) points, instead of the whole data set area.
%   This is extremly fast when (x,y) is within a very small area.

data=load(xyzfile);

bd=boundary(data(:,1),data(:,2));

%get data average increment.
datadx=range(data(:,1))/length(unique(data(:,1)));
datady=range(data(:,2))/length(unique(data(:,2)));

%find indices of the data points within the search area enclosing (x,y) points. 
%Here I extend the search area by 3 times the average data increments. Some points in x
%and y may be within smaller area, there is no need to use the whole xyz
%area.
% max([dx,datadx])
% max([dy,datady])
% min(x)
% min(y)
% max(x)
% max(y)
idx=find(data(:,1) >= (min(x)-3*max([dx,datadx])) & data(:,1) <= (max(x)+3*max([dx,datadx])) & ...
    data(:,2) >= (min(y)-3*max([dy,datady])) & data(:,2) <= (max(y)+3*max([dy,datady])));
% idx=1:size(data,1);
x0=data(idx,1);y0=data(idx,2);z0=data(idx,3);
% min(x0)
% min(y0)
% max(x0)
% max(y0)
% plot for debugging. check points locations.
% plot(x0,y0,'ro'); hold on;
% plot(data(bd,1),data(bd,2),'b-'); 
% plot(x,y,'k+');
% hold off;
x_min=min(x0);x_max=max(x0);
y_min=min(y0);y_max=max(y0);
%min(z0)/1000
%grid points
xx=x_min:dx:x_max;
yy=y_min:dy:y_max;
zz=griddata(x0',y0,z0,xx',yy);
%min(min(zz))/1000
z=nan(length(x),1);

for i=1:length(x)
    %if(x(i) < x_min || x(i) > x_max || y(i) < y_min || y(i) > y_max )
%     inpolygon(x(i),y(i),x0(bd),y0(bd))
    if inpolygon(x(i),y(i),data(bd,1),data(bd,2))==0
        warning('!!Point is outside the coverage of the surface, set z value to nan.');
    else
        [p1,p2,p3,p4]=myfindneighbors(xx,yy,zz,x(i),y(i));
        if(numel(p1)==3 && numel(p2)==3 && numel(p3)==3 && numel(p4)==3 )
            z(i,1)=nanmean([p1(3),p2(3),p3(3),p4(3)]);
        else
            warning('!!Point is outside the coverage of the surface, set z value to nan.');
        end
    end
    %z(i,1)
end

%min_z=min(z)/1000
%max_z=max(z)/1000
xin=x;
yin=y;
return;
end