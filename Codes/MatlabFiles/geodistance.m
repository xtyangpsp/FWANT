function [arclen,az]=geodistance(lat1,lon1,lat2,lon2,depth)
R=6371.0; 
lat1 = lat1*pi/180.0;
lat2 = lat2*pi/180.0;
lon1 = lon1*pi/180.0;
lon2 = lon2*pi/180.0;
  
arclen = acos(cos(lat1).*cos(lat2).*cos(lon1 - lon2) + sin(lat1).*sin(lat2));
  
  az = abs(arclen*180.0/pi);  
  arclen=arclen.*(R-depth);
  
  if( az > 180.0 ) 
      az = 360.0 - az;
      
  end
  
end