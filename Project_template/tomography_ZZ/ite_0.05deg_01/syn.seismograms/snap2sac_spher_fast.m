% snap2sac_spher_fast: read the surface velocity snapshots of the fz run of every
% virtual source in sim.station and write SAC seismograms at all receivers.
% Same output as snap2sac_spher.m (fz run only, ZZ / Rayleigh), but faster:
%   1) block read: many time layers per NetCDF open instead of one gather_snap call
%      per layer (gather_snap opens every subdomain file twice per layer);
%   2) the receiver sampling uses precomputed bilinear weights instead of interp2 at
%      every layer (same result as interp2 'linear'; NaN outside the grid);
%   3) the receiver coordinates are column vectors (the original xrecv=nan(nrecv)
%      made an nrecv x nrecv matrix, so interp2 evaluated nrecv^2 points per layer).
%
% Files written for each pair into syn.seismograms/<src>/ as
%     <src>.to.<rec>.fz.Vz.SAC    velocity      (save_vel)
%     <src>.to.<rec>.fz.Uz.SAC    displacement  (save_disp)
% The self pair (source = receiver) is written, as in the original script.
%
% Needs the same functions as the original script (locate_snap, gather_coord,
% gather_snap, bsac, ch, wsac, distance, nc_varget, nc_getdiminfo, get_fnm_snap).

addpath(genpath('/depot/xtyang/data/codes/mexcdf'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise'))
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))

%% ===================== project =====================
PROJHOME = ['/depot/xtyang/data/projects/xtyang/craton'];
ite = ['/ite_0.025deg_01'];
ite_num = [PROJHOME ite];
sim_dir = [PROJHOME '/' ite '/sim.station'];
syn_dir = [ite_num '/syn.seismograms'];
cd (syn_dir);

MFILE_ROOT=[PROJHOME '/mfiles'];
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);

stnfile = [ PROJHOME '/STinfo/station_riogrande_withdata_formatok.txt'];

%% ===================== options =====================
fcomp = 'fz';            % force run folder
recvcomp = 'Vz';         % particle velocity component in the snapshots
recvcompdisp = 'Uz';     % displacement name in the file names

save_vel   = true;       % write the velocity SAC   (<fcomp>.Vz)
save_disp  = true;       % write the displacement SAC (<fcomp>.Uz, integrated velocity)

fast_read  = true;       % block read (see above); the first layer of every source is checked against
                         % gather_snap. false = original one-layer-at-a-time reader.
tblk       = 100;        % max number of time layers per block read; memory ~ tblk*nx*ny*8 bytes per worker

%% ===================== station info =====================
fid = fopen(stnfile);
recv = textscan(fid,'%s%s%f%f%f'); %ntwk stnm lon lat elevation
fclose(fid);
nrecv=length(recv{1});
ntwk=recv{1}; stnm=recv{2};
slon=recv{3}; slat=recv{4};
xrecv=nan(nrecv,1); yrecv=nan(nrecv,1);
for n=1:nrecv
    xrecv(n)=(90-slat(n))/180*pi;      % colatitude (rad)
    yrecv(n)=slon(n)/180*pi;           % longitude (rad)
end

sim_input = [ ite_num '/sim.input/' ];
pnm_metric = sim_input;
% ----------------------- parameter -----------------------
% snap_002 is the velocity output on the free surface
id=2;subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
n1=1;dn=1;
n2=4000; % steps saved every dn steps. n2 should be <= saved time steps.
layers=n1:dn:n2;
nlayers=length(layers);

nfile_per_pair = double(save_vel) + double(save_disp);
nfile_expected = nfile_per_pair*nrecv;          % the original check: 2 files per receiver
fprintf('%d SAC files per station pair, %d per virtual source\n', nfile_per_pair, nfile_expected);

%% ===================== extract =====================
parfor nstsrc=1:nrecv %
    event_name=[char(ntwk(nstsrc)) '.' char(stnm(nstsrc))];
    OUTPUT_ROOT=[ syn_dir '/' event_name ];
    if ~exist(OUTPUT_ROOT,'dir')
        mkdir(OUTPUT_ROOT);
    end
    %check if already have files
    if length(dir([OUTPUT_ROOT '/*.SAC']))<nfile_expected
        disp([ '--> extracting waveforms for virtual source  ' event_name ]);
    else
        disp([event_name ' already finished. Skip!'])
        continue
    end

    sim_stn = [ sim_dir '/' event_name ];
    sim_fcmp = [ sim_stn '/' fcomp];
    INPUT_ROOT=[ sim_fcmp ];
    fnm_conf=[INPUT_ROOT '/SeisFD3D.conf'];
    pnm_out =[INPUT_ROOT '/output/'];

    % -------------------- load coord --------------------------
    [snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);

    %-- get coord data
    [x,y,z]=gather_coord(snapinfo,'coorddir',pnm_metric);
    xs = squeeze(x);  ys = squeeze(y);
    nx = size(xs,1);  ny = size(xs,2);

    % bilinear sampling weights at the receivers (same result as interp2 'linear', NaN outside)
    xv = xs(:,1);  yv = ys(1,:)';
    rect = max(max(abs(xs - xv))) <= 1e-9 && max(max(abs(ys - yv'))) <= 1e-9 && all(diff(xv) > 0) && all(diff(yv) > 0);
    if rect
        fi = interp1(xv,(1:nx)',xrecv);
        fj = interp1(yv,(1:ny)',yrecv);
        ok = ~isnan(fi) & ~isnan(fj);
        i0 = min(floor(fi(ok)),nx-1);   j0 = min(floor(fj(ok)),ny-1);
        wx = fi(ok) - i0;               wy = fj(ok) - j0;
        l00 = i0 + (j0-1)*nx;  l10 = l00 + 1;  l01 = l00 + nx;  l11 = l00 + nx + 1;
        w00 = (1-wx).*(1-wy);  w10 = wx.*(1-wy);  w01 = (1-wx).*wy;  w11 = wx.*wy;
    else
        warning('%s: surface grid is not rectilinear/increasing; using interp2 per layer (slow).', event_name);
        ok = true(nrecv,1);  l00 = []; l10 = []; l01 = []; l11 = [];
        w00 = []; w10 = []; w01 = []; w11 = [];
    end

    % --------------------- load data --------------------------
    Wz = nan(nlayers,nrecv);  T = nan(nlayers,1);
    tload = tic;
    if fast_read && rect && n1 == 1 && dn == 1
        % ---- block read: several time layers per file open ----
        tcnt = snapinfo(1).tcnt;                       % time layers stored per NetCDF file
        fprintf('  %s: tcnt = %d layers per file (from SeisFD3D.conf), tblk = %d, block = %d layers\n', event_name, tcnt, tblk, min(tblk,tcnt));
        L = 1;
        while L <= nlayers
            nofnc  = fix((L-1)/tcnt) + 1;
            ninnc0 = mod(L-1,tcnt);
            nt = min([tblk, tcnt-ninnc0, nlayers-L+1]);
            [B,tv] = read_snap_block(snapinfo,id,nofnc,ninnc0,nt,{recvcomp},pnm_out,nx,ny);
            if L == 1                                  % check against the original reader
                [Vchk,tchk] = gather_snap(snapinfo,id,1,recvcomp,'outdir',pnm_out);
                Vchk = squeeze(Vchk);
                Vnew = reshape(B(1,:,:,1),nx,ny);
                if max(abs(Vchk(:)-Vnew(:))) > 1e-6*max(abs(Vchk(:))) || abs(tchk-tv(1)) > 1e-9
                    error('%s: block read differs from gather_snap; set fast_read=false.', event_name);
                end
            end
            for jj = 1:nt
                V2 = reshape(B(jj,:,:,1),nx,ny);
                vq = nan(nrecv,1);
                vq(ok) = w00.*V2(l00) + w10.*V2(l10) + w01.*V2(l01) + w11.*V2(l11);
                Wz(L+jj-1,:) = vq';
                T(L+jj-1) = tv(jj);
            end
            L = L + nt;
            fprintf('  %s: %d/%d layers read (%.0f s)\n', event_name, L-1, nlayers, toc(tload));
        end
    else
        for j=1:nlayers
            ly=layers(j);
            t = NaN;                          % set inside the parfor iteration (silences the temporary-variable warning)
            [Vsnap,t]=gather_snap(snapinfo,id,ly,recvcomp,'outdir',pnm_out);
            V2 = squeeze(Vsnap);
            if rect
                vq = nan(nrecv,1);
                vq(ok) = w00.*V2(l00) + w10.*V2(l10) + w01.*V2(l01) + w11.*V2(l11);
            else
                vq = interp2(ys,xs,V2,yrecv,xrecv,'linear');
            end
            Wz(j,:) = vq(1:nrecv)';
            T(j) = t;
            if mod(j,500) == 0 || j == nlayers
                fprintf('  %s: %d/%d layers read (%.0f s)\n', event_name, j, nlayers, toc(tload));
            end
        end  % end of nlayer
    end

    tt=[0; T];dt=tt(2)-tt(1);
    for nstrec=1:nrecv
        vzrecv=[0 Wz(:,nstrec)'];
        [disttemp,~]=distance(slat(nstsrc),slon(nstsrc),slat(nstrec),slon(nstrec));
        recvname = [ char(ntwk(nstrec)) '.' char(stnm(nstrec)) ];
        if save_vel
            Szrecv=bsac(tt,vzrecv);
            Szrecv=ch(Szrecv,'STLA',slat(nstrec),'STLO',slon(nstrec),'EVLA',slat(nstsrc),'EVLO',slon(nstsrc),...
                'DIST',disttemp);
            wsac([OUTPUT_ROOT '/' event_name '.to.' recvname '.' fcomp '.' recvcomp '.SAC'],Szrecv);
        end
        if save_disp
            % for displacement
            uzrecv=cumtrapz(vzrecv)*dt;
            Szrecv=bsac(tt,uzrecv);
            Szrecv=ch(Szrecv,'STLA',slat(nstrec),'STLO',slon(nstrec),'EVLA',slat(nstsrc),'EVLO',slon(nstsrc));
            wsac([OUTPUT_ROOT '/' event_name '.to.' recvname '.' fcomp '.' recvcompdisp '.SAC'],Szrecv);
        end
    end  % end of receiver

end  % end of source

%% ===================== local functions =====================
function [B,tv] = read_snap_block(snapinfo,id,nofnc,ninnc0,nt,comps,pnm_out,nx,ny)
% Surface velocity of nt consecutive time layers (layers ninnc0 .. ninnc0+nt-1 of file number nofnc)
% from all subdomain files, assembled on the (nx,ny) grid:  B(layer,i,j,component), plus times tv.
% Same file layout and indexing as gather_snap, but one open per file and component for nt layers
% instead of one per layer.  (Free-surface snapshot: the vertical dimension must be 1.)
nthd = length(snapinfo);
ncomp = numel(comps);
B = zeros(nt,nx,ny,ncomp);
tv = [];
for n = 1:nthd
    n_i=snapinfo(n).thisid(1); n_j=snapinfo(n).thisid(2); n_k=snapinfo(n).thisid(3);
    i1=snapinfo(n).indxs(1); j1=snapinfo(n).indxs(2); k1=snapinfo(n).indxs(3);
    i2=snapinfo(n).indxe(1); j2=snapinfo(n).indxe(2); k2=snapinfo(n).indxe(3);
    subs=snapinfo(n).subs; subc=snapinfo(n).subc; subt=snapinfo(n).subt;
    if k2 ~= k1, error('read_snap_block: expected a free-surface snapshot (one vertical index), set fast_read=false.'); end
    fnm_snap = get_fnm_snap(pnm_out,'vel_',n_i,n_j,n_k,id,nofnc);
    if ~exist(fnm_snap,'file'), error('read_snap_block: file %s does not exist',fnm_snap); end
    tdim = nc_getdiminfo(fnm_snap,'time');
    if tdim.Length < ninnc0+nt
        error('read_snap_block: layers %d-%d are beyond the time dim (%d) in %s',ninnc0+1,ninnc0+nt,tdim.Length,fnm_snap);
    end
    ni = i2-i1+1;  nj = j2-j1+1;
    for c = 1:ncomp
        d = nc_varget(fnm_snap,comps{c},[ninnc0,fliplr(subs)-1],[nt,fliplr(subc)],[1,fliplr(subt)]);
        d = reshape(d,[nt,nj,ni]);                 % (time, j, i); the vertical dimension is 1
        B(:,i1:i2,j1:j2,c) = permute(d,[1 3 2]);   % (time, i, j)
    end
    if n == 1, tv = nc_varget(fnm_snap,'time',ninnc0,nt); tv = tv(:); end
end
end
