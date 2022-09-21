function [x,y,z,d]=surface2profile(xyzfile,x1,y1,x2,y2,dx,dy,npoints,geo,ztype)
%plot profiles.
%okay for small regions. approximate when computing points between two end
%points.
%USAGE: [x,y,z,d]=surface2profile(xyzfile,x1,y1,x2,y2,dx,dy,npoints,geo,ztype)
%
%PARAMETERS:
%dx,dy: grid size when griding data.
%npoints: number of points along the profile.
%geo: yes or no string.
%ztype: elevation (+ above sealevel) or depth (+ below sealevel), or other
%(will not be subtracted from Earth's radius when calculating distance.)
%xyzfile: text file with x,y,z values in each column;
%x1,y1 and x2,y2 are the two end points.
%
%RETURN VALUES: x,y,z are the point values along the profile. d is the
%distance in data unit if geo is no, otherwise in great circle distance
%(km).
%Xiaotao Yang @ Indiana University, September 21, 2015
%
% History:
%   1. May 18, 2018: For geo='yes', use great circle path to get the points between the
%   starting and ending points.
%
if(nargin<8)
    error('ERROR: not enough arguments.');
elseif(nargin==8)
    geo='yes';
    ztype='elevation';
elseif(nargin==9)
    ztype='elevation';
end

display(['Working on ',xyzfile,' ...']);
%grid data.
data=load(xyzfile);

grid_size_x=dx;
grid_size_y=dy;

x0=data(:,1);y0=data(:,2);z0=data(:,3);
x_min=min(x0);x_max=max(x0);
y_min=min(y0);y_max=max(y0);

%grid points
xx=x_min:grid_size_x:x_max;
yy=y_min:grid_size_y:y_max;
zz=griddata(x0',y0,z0,xx',yy);

%get points along the profile.
if strcmp(geo,'yes') %get great circle path points.
    [y,x]=track2(y1,x1,y2,x2,wgs84Ellipsoid,'degrees',npoints);
else
    if(abs(x2-x1)>0)
        ddx=(x2-x1)/npoints;
        x=x1:ddx:x2;
        if(y2==y1) y=y1*ones(length(x),1);
        else
           slope=(y2-y1)/(x2-x1);
           intercept=-x2*slope + y2;
           y=slope*x + intercept;
        end
    elseif(abs(y2-y1)>0)
        ddy=(y2-y1)/npoints;
        y=y1:ddy:y2;
        x=x1*ones(length(y));
    else
        error('Error: The given two points are the same!');
    end
end
z=nan(length(x),1);

for i=1:length(x)
    if(x(i) >=x_min && x(i) <=x_max && y(i) >=y_min && y(i) <=y_max)
        [p1,p2,p3,p4]=myfindneighbors(xx,yy,zz,x(i),y(i));

        z(i)=interpolate2d(p1,p2,p3,p4,x0,y0,'average');
    end
end

if(strcmp(geo,'yes'))
    
    d=zeros(length(x),2);
    for j=2:length(x)
       if(strcmp(ztype,'elevation'))
           [d(j,1),d(j,2)]=geodistance(y(j),x(j),y(1),x(1),-1.0*nanmean(z)); 
       elseif(strcmp(ztype,'depth'))
           [d(j,1),d(j,2)]=geodistance(y(j),x(j),y(1),x(1),nanmean(z));
       elseif(strcmp(ztype,'other'))
           [d(j,1),d(j,2)]=geodistance(y(j),x(j),y(1),x(1),0.0);
       else
           error('ERROR: unknow ztype, can only be one of elevation, depth, or other.');
       end
    end
end

return;
end

function z0=interpolate2d(p1,p2,p3,p4,x0,y0,method)
%method: average, weightsum

if(strcmp(method,'average'))
    if(numel(p1) ==3 && numel(p2) ==3 && numel(p3) ==3 && numel(p4) ==3)
        z0=nanmean([p1(3),p2(3),p3(3),p4(3)]);
    else
        z0=nan;
    end
end

return;
end