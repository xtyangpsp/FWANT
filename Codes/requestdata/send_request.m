

% this program is to send data request to Breqfast

clear all;

networks={'US'};

for idx = 1:length(networks)

nw = networks{idx};

%eval([ 'cd  ./' nw '_request' ]);

eval([ 'cd  ./' nw '_2000_2016' ]);

for yy = 2000:2017
for mm = 1:12
    
fname = [nw '_Small_' num2str(yy) '_' num2str(mm) ];
fid = fopen(fname,'r');
if fid<0, disp('no file found!'),continue, end
fclose(fid);

unix( [' mail breq_fast@iris.washington.edu < ' fname ]);  
%unix( [' mail  miniseed@iris.washington.edu < ' fname ]);  un@iris.washington.edu
display(['request ' fname ' is sent']);

pause(2);   % sleep the request for half an hour

end
end

cd ../

end


