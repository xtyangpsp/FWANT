

function cmp = smoothing_hori(XY,index)

pnt = find(isnan(XY));
XY(pnt) = 0;
temp = zeros(size(XY));
[xx,yy] = size(XY);


switch index
    
  case 1,  % nearest 4 points smoothing
    XYtmp = zeros(xx+2,yy+2);
    XYtmp(2:end-1,2:end-1)=XY;
    
    wts = [0.4 0.15*ones(1,4)]';
    for i = 2:xx+1
    for j = 2:yy+1   
        cmps = [XYtmp(i,j) XYtmp(i,j-1) XYtmp(i,j+1) XYtmp(i+1,j) XYtmp(i-1,j)];
        cc=find(cmps==0); 
        if ~isempty(cc), wts(cc)=0; end
          
        sum_cmps = cmps*wts;
        sum_wts = sum(wts);                                             
        temp(i-1,j-1) = sum_cmps/sum_wts;
    end
    end
  
  case 2,  % nearest 8 points smoothing   
    XYtmp = zeros(xx+2,yy+2);
    XYtmp(2:end-1,2:end-1)=XY;
    
    wts = [0.4 0.1*ones(1,4) 0.05*ones(1,4)]';
    for i = 2:xx+1
    for j = 2:yy+1   
        cmps = [XYtmp(i,j) XYtmp(i-1,j) XYtmp(i,j-1) XYtmp(i,j+1) XYtmp(i+1,j) ...
                XYtmp(i-1,j-1) XYtmp(i+1,j-1) XYtmp(i-1,j+1) XYtmp(i+1,j+1)];
        cc=find(cmps==0); 
        if ~isempty(cc), wts(cc)=0; end
          
        sum_cmps = cmps*wts;
        sum_wts = sum(wts);                                             
        temp(i-1,j-1) = sum_cmps/sum_wts;
    end
    end  
    
  case 3,  % nearest 24 points smoothing   
    XYtmp = zeros(xx+4,yy+4);
    XYtmp(3:end-2,3:end-2)=XY;
    
    wts = [0.3 0.1*ones(1,4) 0.035*ones(1,4) 0.01*ones(1,16)]';
    for i = 3:xx+2
    for j = 3:yy+2   
        cmps = [XYtmp(i,j) XYtmp(i-1,j) XYtmp(i,j-1) XYtmp(i,j+1) XYtmp(i+1,j) ...
                XYtmp(i-1,j-1) XYtmp(i+1,j-1) XYtmp(i-1,j+1) XYtmp(i+1,j+1) ...
                XYtmp(i-2,j-2) XYtmp(i-2,j-1) XYtmp(i-2,j) XYtmp(i-2,j+1) XYtmp(i-2,j+2) ...
                XYtmp(i-1,j-2) XYtmp(i-1,j+2) XYtmp(i,j-2) XYtmp(i,j+2) XYtmp(i+1,j-2) XYtmp(i+1,j+2) ...
                XYtmp(i+2,j-2) XYtmp(i+2,j-1) XYtmp(i+2,j) XYtmp(i+2,j+1) XYtmp(i+2,j+2)];
        cc=find(cmps==0); 
        if ~isempty(cc), wts(cc)=0; end
          
        sum_cmps = cmps*wts;
        sum_wts = sum(wts);                                             
        temp(i-2,j-2) = sum_cmps/sum_wts;
    end
    end 
    
end
      
    
        
cmp = temp;

pnt1 = find(cmp==0);
cmp(pnt1) = NaN;

return
    
 