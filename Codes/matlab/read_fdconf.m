function [confinfo]=read_fdconf(fnm_conf,id)
% Read major configuration parameters from SeisFD3D.conf
%
% Xiaotao Yang

% check
if ~ exist(fnm_conf,'file')
   error([mfilename ': file ' fnm_conf ' does not exist']);
end

str_cm='#';
stag=['snap_' num2str(id,'%3.3i')];
fid=fopen(fnm_conf);
conf=textscan(fid,'%s','delimiter','\n','whitespace','');
fclose(fid);
nline=size(conf{1});
confinfo=struct();
confinfo.stag=stag;
for n=1:nline

    str=conf{1}{n};
    if isempty(str)
       continue
    end
    npit=strfind(str,str_cm);
    if ~ isempty(npit)
       str(npit(1):end)=[];
    end
    
    str=regexprep(str,{'=','\|'},' ');
    [tag,s]=strtok(str);
    if isempty(tag)
       continue
    end

    switch tag
    case 'dims'
        confinfo.dims=sscanf(s,'%f',3)';
    case 'ni'
        confinfo.ni=sscanf(s,'%f',1);
    case 'nj'
        confinfo.nj=sscanf(s,'%f',1);
    case 'nk'
        confinfo.nk=sscanf(s,'%f',1);
    case 'nt'
        confinfo.nt=sscanf(s,'%f',1);
    case 'GRID_ROOT'
        confinfo.pnm_grid=sscanf(s,'%s',1);
    case 'MEDIA_ROOT'
        confinfo.pnm_media=sscanf(s,'%s',1);
    case 'SOURCE_ROOT'
        confinfo.pnm_src=sscanf(s,'%s',1);
    case 'OUTPUT_ROOT'
        confinfo.pnm_output=sscanf(s,'%s',1);
    case 'number_of_recv'
        confinfo.num_recv=sscanf(s,'%f',1);
    case 'number_of_inline'
        confinfo.num_line=sscanf(s,'%f',1);
    case 'number_of_snap'
        confinfo.num_snap=sscanf(s,'%f',1);
    case stag
        [a,~,~,npt]=sscanf(s,'%f',11);
        snap_subs=a(1:3)';
        snap_subc=a(4:6)';
        snap_subt=a(7:9)';
        snap_tinv=a(10);
        snap_tcnt=a(11);
        snap_stress=0;
        b=sscanf(s(npt:end),'%s');
        if ~ isempty(b)
            if strcmpi(b,'t')
               snap_stress=1;
            end
        end
    end  % select

end

confinfo.ngi=confinfo.ni*confinfo.dims(1);
confinfo.ngj=confinfo.nj*confinfo.dims(2);
confinfo.ngk=confinfo.nk*confinfo.dims(3);
confinfo.ngijk=[confinfo.ngi,confinfo.ngj,confinfo.ngk];

if id==0
    snap_subs=[  1, 1, 1 ];
    snap_subc=[ -1,-1,-1 ];
    snap_subt=[  1, 1, 1 ];
    snap_tinv=1;
    snap_tcnt=1;
end

confinfo.snap_subs=snap_subs;
confinfo.snap_subc=snap_subc;
confinfo.snap_subt=snap_subt;
confinfo.snap_tinv=snap_tinv;
confinfo.snap_tcnt=snap_tcnt;
if id>0; confinfo.snap_stress=snap_stress;end

if ~ exist('snap_subs','var')
    error([ 'id=',num2str(id),' doesn''t exist']);
end

end