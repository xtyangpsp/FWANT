
% this program is to generate request files

clear all;

nw={'Z4'};

%[network, sta,slat,slon, chan] = textread([char(nw) '_2014_2015_SIO.txt' ], '%s %s %*s %f %f %*f %*s %*s %s');
%[sta,slat,slon, chan] = textread([char(nw) '_EXZ_2015.txt' ], '%s %*s %f %f %*f %*s %*s %s');
%[network, sta, slon, slat] = textread([char(nw) '_formatted.txt' ], '%s %s %f %f %*f');
 [network,sta,slat,slon,chan] = textread([char(nw) '.txt' ], '%s %s %f %f %s');

for yy = 2010:2010
for mm = 1:12
    
fname = [char(nw) '_Small_' num2str(yy) '_' num2str(mm) ]
fid = fopen(fname,'w');

fprintf(fid, '%s \n', '.NAME Cong Li ');
fprintf(fid, '%s \n', '.INST UMASS AMHERST ');
fprintf(fid, '%s \n', '.EMAIL qingcongll@163.com ');
fprintf(fid, '%s \n', ['.LABEL ' char(nw) '_Small_'  num2str(yy) '_' num2str(mm) ]);

fprintf(fid, '%s \n', '.PHONE 413-404-6763 ');
fprintf(fid, '%s \n', '.MEDIA FTP ');
fprintf(fid, '%s \n', '.ALTERNATE MEDIA ');
fprintf(fid, '%s \n', '.ALTERNATE MEDIA ');
fprintf(fid, '%s \n', '.FROM_JWEED ');
fprintf(fid, '%s \n', '.END ');
% 
for ite = 1:length(sta)
    if mm<=11
    fprintf(fid, ' %s %s %d %d %s %s %s %s %d %d %s %s %s %s %s %s %s \n', char(sta(ite)), char(nw), yy, mm, '01', '00', '00', '00', yy, mm+1, '01', '00', '00', '00', '1', char(chan(ite)), '*' );
    elseif mm==12
    fprintf(fid, ' %s %s %d %d %s %s %s %s %d %d %s %s %s %s %s %s %s \n', char(sta(ite)), char(nw), yy, mm, '01', '00', '00', '00', yy+1, 1, '01', '00', '00', '00', '1', char(chan(ite)), '*' );
    end
end

% for ite = 1:length(sta)
% fprintf(fid, ' %s %s %d %s %s %s %s %s %d %s %s %s %s %s %s %s %s \n', char(sta(ite)), char(network(ite)), yy, '01', '01', '00', '00', '00', yy+1, '01', '01', '00', '00', '00', '1', 'L*', '*' );
% end
unix(['mkdir ' char(nw) '_2000_2016/']);
unix([ ' mv  ' fname ' ./' char(nw) '_2000_2016/' ]);
%unix([ ' mv  ' fname ' ./CRB/' char(nw) '_request' ]);
fclose(fid);

end
end
