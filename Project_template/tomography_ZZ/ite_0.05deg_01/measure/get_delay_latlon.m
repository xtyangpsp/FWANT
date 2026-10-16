function get_delay_latlon(stationinfofile,rawdelaydatafilelist)
%%%%
%%%%
numinput=numel(rawdelaydatafilelist);

for k=1:numinput
    rawdelaydatafile=rawdelaydatafilelist{k};
    tempdelaydatafile=strcat(rawdelaydatafile, '_latlon');
    %readin station info
    fidsta=fopen(stationinfofile,'r');
    stations=textscan(fidsta,'%s %s %f %f');
    fclose(fidsta);
    [nw,stn,lon,lat]=stations{1:4};

    %readin raw data
    fiddata=fopen(rawdelaydatafile,'r');
    rawdata=textscan(fiddata,'%s %f %*d %*s %*f %*f %s %f %f %f');
    fclose(fiddata);
    clear files dt fb err xcoe xsnr;
    [files,dt,fb,err,xcoe,xsnr]=rawdata{1:6};
    clear stations rawdata;
    %
    disp(strcat('Writing to: [ ',tempdelaydatafile,' ] ...'));
    fidtemp = fopen(tempdelaydatafile,'w');

    nfiles = length(dt);
    for ii = 1:nfiles
        tmpfile = char(files{ii});
        idxs = strfind(tmpfile,'/');
        src = tmpfile(1:idxs(1)-1);
        idxsd = strfind(src,'.');
        srcntwk = src(1:idxsd-1);
        srcstn = src(idxsd+1:end);

        idxr = strfind(tmpfile,'_');
        rcv = tmpfile(idxr(2)+1:idxr(3)-1);
        idxrd = strfind(rcv,'.');
        rcvntwk = rcv(1:idxrd-1);
        rcvstn = rcv(idxrd+1:end);

        clear ids idr
        ids = find(strcmp(srcntwk,nw) & strcmp(srcstn,stn));
        idr = find(strcmp(rcvntwk,nw) & strcmp(rcvstn,stn));
        if ~isempty(ids) && ~isempty(idr)
            lon1 = lon(ids);
            lat1 = lat(ids);
            lon2 = lon(idr);
            lat2 = lat(idr);

            fprintf(fidtemp,'%s %s %f %f %f %f %f %f %s %f %f\n',char(src), char(rcv), ...
                lat1, lon1, lat2, lon2, dt(ii), err(ii),char(fb{ii}),xcoe(ii),xsnr(ii));   
        end
    end
    fclose(fidtemp);
    % end of writing to temp data file.
end