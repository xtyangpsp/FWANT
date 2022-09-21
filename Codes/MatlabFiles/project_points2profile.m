function [dist0,dd,xx,yy,zz,mm,cc]=project_points2profile(xyzfile,x1,y1,x2,y2,dmax,geo,ztype)
%project xyz points onto surface with given limit dmax as the maximum
%distance away from the profile.
%USAGE: [dist,dd,xx,yy,zz,mm]=project_points2profile(xyzfile,x1,y1,x2,y2,dmax,geo,ztype)
%xyzfile: file containing xyz points OR a N by 3 matrix containing the X,
%Y, Z data.
%x1,y1 and x2,y2 are the two end points.
%geo: yes or no string.
%ztype: elevation (+ above sealevel) or depth (+ below sealevel), or other
%dmax: stated above, in kilometers if geo is yes.
%
%RETURNS:
%distance: distance away from the start point, km when geo is yes.
%xx: x value of the output points
%yy: y value of the output points
%zz: z value of the points satisfy the cutoff;
%dd: distance of the points to the profile.
%mm: the 4th column, usually magnitue for earthquakes.
%
% External functions:
% [lat,lon] = gcxgc(lat1,lon1,az1,lat2,lon2,az2)%find intersections of the
% two great circles.
% distance()
%
% History:
%   May 24, 2018 by Xiaotao Yang: use great circle for geographical points. this is very
%   important for high lattitude regions. Change outputs to: [dist,dd,xx,yy,zz,mm]
%
%Xiaotao Yang @ Indiana University, October 8, 2015

if(nargin<6)
    error('ERROR: not enough arguments.');
elseif(nargin==6)
    geo='yes';
    ztype='elevation';
elseif(nargin==7)
    ztype='elevation';
end

R=6371.0;
if ischar(xyzfile)
    data=load(xyzfile);
else
    data=xyzfile; %xyzfile could be an array.
end
x=data(:,1); y=data(:,2);z=data(:,3);
column4=0;
if size(data,2) >=4
    column4=1;
    m=data(:,4);
end
column5=0;
if size(data,2) >=5
    column5=1;
    c=data(:,5);
end
%if geographical points, use track2() to get the points along great circle
% otherwise, use cartisian interpolation.
%
intersection_x=nan(length(x),1);
intersection_y=nan(length(x),1);
d=nan(length(x),1);
if(strcmp(geo,'yes'))
    az1=azimuth(y1,x1,y2,x2);
    for i=1:length(x)
        clear interlat interlon;
        [interlat,interlon]=gcxgc(y1,x1,az1,y(i),x(i),az1+90);
%         interlat
%         y(i)
        clear idxlat idxlon;
        if x1==x2 %along maridian
            idxlat=find(interlat>=min(y1,y2) & interlat <=max(y1,y2));
            if ~isempty(idxlat)
                intersection_x(i)=interlon(idxlat);
                intersection_y(i)=interlat(idxlat);
                d(i)=distance(y(i),x(i),intersection_y(i),intersection_x(i),wgs84Ellipsoid)/1000;
            end
        else
            idxlat=find(interlat>=min(y1,y2) & interlat <=max(y1,y2));
            idxlon=find(interlon>=min(x1,x2) & interlon<=max(x1,x2));
            if ~isempty(idxlat) && ~isempty(idxlon)
                intersection_x(i)=interlon(idxlon);
                intersection_y(i)=interlat(idxlat);
                d(i)=distance(y(i),x(i),intersection_y(i),intersection_x(i),wgs84Ellipsoid)/1000;
            end
        end
    end
%     d=distance(y,x,intersection_y,intersection_x,wgs84Ellipsoid)/1000;
else
    if(x1==x2)  %vertical line
        intersection_x=x1*ones(length(x),1);
        intersection_y=y;
        d=abs(x-x1);
    elseif(y1==y2) %horizontal line
        intersection_x=x;
        intersection_y=y1*ones(length(y),1);
        d=abs(y1-y);
    else
        slope=(y2-y1)/(x2-x1);
        intercept=-x2*slope + y2;

        slope_points=-1.0/slope;
        intercept_points=y-slope_points*x;

        intersection_x=(intercept_points-intercept)./(slope - slope_points);
        intersection_y=slope*(intercept_points-intercept)./(slope - slope_points)+intercept;
        d=sqrt((x-intersection_x)^2 + (y-intersection_y)^2);
    end
end
    %sometimes, there is nothing projected. assign a value for the first
%element is neccessary. 
% distance(1)=nan; %distance of each point from the starting point.
% zz(1)=nan;
% dd(1)=nan; %distance array of the finally saved points
% j=1;
clear idx0;
xx=[];yy=[];zz=[];dd=[];dist0=[];mm=[];cc=[];
idx0=find(d<=dmax);
if ~isempty(idx0)
    xx=intersection_x(idx0);
    yy=intersection_y(idx0);
    zz=z(idx0);
    dd=d(idx0);
    if column4
        mm=m(idx0);
    else
        mm=nan(length(idx0),1);    
    end
    if column5
        cc=c(idx0);
    else
        cc=nan(length(idx0),1);    
    end
    if(strcmp(ztype,'elevation'))
        dist0=geodistance(intersection_y,intersection_x,y1,x1,-1.0*z);
    elseif(strcmp(ztype,'depth'))
        dist0=geodistance(intersection_y,intersection_x,y1,x1,z);
    elseif(strcmp(ztype,'other'))
        dist0=geodistance(intersection_y,intersection_x,y1,x1,0.0);
    end
end

% if ~column4; mm=nan(length(xx),1);end
return;  
end