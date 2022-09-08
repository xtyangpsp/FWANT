% input CUB2.0 model
function [outlat, outlon, outdep, outvs]=read_CUB_revised(minlon,maxlon,minlat,maxlat)
define_latlon;
% SDT: diffraction tomography
% SRT: ray tomography.  See README in the directory
file='./vmodels/CUB2/CU_SDT1.0.mod';
outdep=0:4:396; % 4 km per grid from 0 to 396 km depth in CUB2.0

fid=fopen(file,'r');
data=fscanf(fid,'%f\n');
fclose(fid);
npoint=length(data)/(2+ length(outdep)*10); %(176/2+1)*180 = 16020 point at 2x2 degrees (88S to 88N)

% initialize the output grid dimensions
ilonlat=0;ilat=0;ilon=0;

% format at each point
% header: lat lon
% 100 lines and 10 columns (depth vs vsv vsh vsmin vsvmin vshmin vsmax vsvmax vshmax)
% 	              radial anisotropy not used

for i=1:npoint
  vel(i).lat=data((i-1)*1002+1);
  vel(i).lon=data((i-1)*1002+2);
  for j=1:length(outdep)
    vel(i).depth(j)=data((i-1)*1002+2+(j-1)*10+1);
    vel(i).vs(j)=data((i-1)*1002+2+(j-1)*10+2);
  end
  % determine if the point is within the study area
  if(vel(i).lon>=minlon && vel(i).lon<=maxlon && ...
     vel(i).lat>=minlat && vel(i).lat<=maxlat )
     
    if(ilat==0)
      ilat=ilat+1;
      outlat(ilat)=vel(i).lat;
    else
      if(vel(i).lat~=outlat(ilat))
        ilat=ilat+1;
        outlat(ilat)=vel(i).lat;
        ilon=0;
      end
    end
    
    ilon=ilon+1;
    outlon(ilon)=vel(i).lon;

    ilonlat=ilonlat+1;
    outvs(ilat,ilon,:)=vel(i).vs;
  end
end

end

