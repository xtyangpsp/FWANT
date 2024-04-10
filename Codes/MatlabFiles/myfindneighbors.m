function [p1,p2,p3,p4]=myfindneighbors(x,y,z,x0,y0)
%Description: Find the 4 neighbors of the given points in a x,y,z dataset.
%USAGE: [p1,p2,p3,p4]=myfindneighbors(x,y,z,x0,y0)
%   p1 (UL), p2(UR), p3 (LL), P4(LR). They all include 3 elements (x,y,z)
%   x,y,z are the three vectors containing the original data.
%   x0,y0 are the targetting point, which is the point we are looking for
%   heighbors for.
%
%AUTHOR: Xiaotao Yang @ Indiana University, September 21, 2015
%

%find x values.
as=find(x<=x0);

[m,c]=max(x(as));

x1=m; x3=x1;cx13=as(c);


as1=find(x>=x0);
[m,c]=min(x(as1));
x2=m; x4=x2;cx24=as1(c);

%find y values.
%find the colum of data with same x values.
%asx13=find(x==x1);
%y_x13=y(asx13); %z_x13=z(asx13); %subset y, and z for the same x values.
%asx24=find(x==x2);
%y_x24=y(asx24)
%z_x24=z(asx24)


as2=find(y>=y0);
[m,c]=min(y(as2));
y1=m;y2=y1;cy12=as2(c);
z1=z(cy12,cx13);
z2=z(cy12,cx24);

as3=find(y<=y0);
[m,c]=max(y(as3));
y3=m;y4=y3;cy34=as3(c);
z3=z(cy34,cx13);
z4=z(cy34,cx24);
%c=find(y_x24==y2);
%z2=z_x24(c);

% as3=find(y_x13<=y0);
% [m,c]=max(y_x13(as3));
% y3=m;y4=y3;cy34=as3(c);
% z3=z_x13(c);
% c=find(y_x24==y4);
% z4=z_x24(c);
% p1=nan(1,3)
% p2=nan(1,3)
% p3=nan(1,3)
% p4=nan(1,3)
% 

p1=[x1,y1,z1]; 
p2=[x2,y2,z2];
p3=[x3,y3,z3];
p4=[x4,y4,z4];

return;
end