% snap2sac_sphere_fxfy: read the surface velocity snapshots of the fx and fy
% runs of every virtual source in sim.station and write SAC seismograms at all
% receivers (station pairs).  Modified from snap2sac_spher.m (fz run only).
%
% Layout expected (as created by skel2station.sh with COMPS="fx fy"):
%     sim.station/<NET.STA>/fx/{SeisFD3D.conf,SeisSource.conf,output/}
%     sim.station/<NET.STA>/fy/{...}
%
% Code frame: x = colatitude direction (SOUTH), y = EAST, z = up.
% For a source A and receiver B let alpha = azimuth A->B and beta = azimuth B->A
% (degrees clockwise from north).  Unit vectors in the (x,y) frame:
%     source side   R_s = (-cos a,  sin a)   T_s = ( sin a,  cos a)
%     receiver side R_r = ( cos b, -sin b)   T_r = (-sin b, -cos b)
% A force along unit vector f gives u(f) = f(1)*u_fx + f(2)*u_fy.  Each run is
% first scaled to "force = f0 of the fx run" so amplitudes match the old fz
% runs (the kernel code assumes SGT_M0 = f0).
%
% Files written for each pair into syn.seismograms/<src>/ as
%     <src>.to.<rec>.<fcomp>.<comp>.SAC
%
%   fT.VT / fT.UT   TT : T_r . u(T_s)                    velocity / displacement  (Love, for phase delays)
%   fR.VR / fR.UR   RR : R_r . u(R_s)                    (optional, save_RR)
%   fT.Vx fT.Vy fT.Vz   response to the SOURCE-SIDE TRANSVERSE FORCE, NOT rotated
%                   at the receiver, in the code frame x,y,z.  These are the
%                   three synthetic components the kernel code reads (it rotates
%                   them itself with rot_mat), so no extra rotation back is needed.
%                   Vx is the colatitude (south) component, not north.
%   (do_rotate = false, or save_raw = true)  fx.Vx fx.Vy fx.Vz fy.Vx fy.Vy fy.Vz
%                   the unrotated responses of the two runs.
%
% TT and (RT = R_r.u(T_s)) are equivalent to fT.Vx/Vy up to a rotation by beta,
% so storing x,y is the simplest choice; it is also what the kernel needs.
%
% The self pair (source = receiver) has no azimuth and is skipped.
% An <src>.azimuth.txt per source lists alpha, beta and distance for every receiver
% (also AZ/BAZ in the SAC headers) for the kernel-configuration step.
%
% Needs MATLAB R2016b or newer.

addpath(genpath('/depot/xtyang/data/codes/mexcdf'))
addpath(genpath('/depot/xtyang/data/codes/MatNoise'))
addpath(genpath('/depot/xtyang/data/codes/FWANT/Codes/MatlabFiles'))

%% ===================== project (CHANGE for the TT project) =====================
PROJHOME = ['/depot/xtyang/data/projects/xtyang/RioGrande/FWANT_TT'];
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
run_x = 'fx';            % run folder with the unit force toward +x (south)
run_y = 'fy';            % run folder with the unit force toward +y (east)

% Defaults = only what the Love (TT) imaging workflow needs:
%   fT.UT         TT displacement      -> phase-delay measurement (synthetic)
%   fT.Vx, fT.Vy  x,y velocity for T_s -> kernel step (kern_synthetic4kernel_fxfy.m)
%   <src>.azimuth.txt                  -> kernel configuration (always written)
do_rotate     = true;    % true: build TT/RR/transverse-force responses; false: save raw fx and fy only
save_disp_TT  = true;    % TT (and RR) displacement SAC: fT.UT   [needed]
save_XYZ      = true;    % fT.Vx, fT.Vy   [needed by the kernel step; only used if do_rotate]
save_vel_TT   = false;   % TT (and RR) velocity SAC: fT.VT  (not used downstream)
save_RR       = false;   % also build RR (fR.UR / fR.VR) for quality checks only
save_Z        = false;   % fT.Vz (only the original SI_ker_sta.F90 reads it; SI_ker_sta_fxyz.F90 does not)
save_raw      = false;   % additionally save the unrotated fx/fy responses (always on if do_rotate = false)
save_disp_XYZ = false;   % displacement SAC for fT.x/y/z and raw (kern_synthetic4kernel integrates velocity itself)

fast_read     = true;    % read many time layers per NetCDF open (block read) instead of one gather_snap call per layer and component;
                         % the first layer of every source is checked against gather_snap. false = original slow path.
tblk          = 100;     % max number of time layers per block read (memory ~ tblk*nx*ny*ncomp doubles per worker)

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
n2=4000; % steps saved every dn steps. n2 should be <= saved time steps.  (must equal what your fz run used: the snapshot reading time scales with it)
layers=n1:dn:n2;
nlayers=length(layers);

% azimuth (deg clockwise from north) at point 1 toward point 2 on a sphere
sph_az = @(la1,lo1,la2,lo2) mod(atan2d(sind(lo2-lo1).*cosd(la2), ...
         cosd(la1).*sind(la2)-sind(la1).*cosd(la2).*cosd(lo2-lo1)),360);

%% ===================== output specification =====================
if ~do_rotate, save_raw = true; end
needZ = save_Z || save_raw;
if needZ, ncomp = 3; else, ncomp = 2; end

% rows: { fcomp , velocity name , displacement name ('' = none) , series id }
%   id 1 = TT, 2 = RR, 3/4/5 = x/y/z response to T_s force, 6..11 = raw (6 + (run-1)*3 + comp-1)
spec = cell(0,4);
if do_rotate
    spec(end+1,:) = {'fT',dname('VT',save_vel_TT),dname('UT',save_disp_TT),1};
    if save_RR,  spec(end+1,:) = {'fR',dname('VR',save_vel_TT),dname('UR',save_disp_TT),2}; end
    if save_XYZ
        spec(end+1,:) = {'fT','Vx',dname('Ux',save_disp_XYZ),3};
        spec(end+1,:) = {'fT','Vy',dname('Uy',save_disp_XYZ),4};
        if save_Z, spec(end+1,:) = {'fT','Vz',dname('Uz',save_disp_XYZ),5}; end
    end
end
if save_raw
    rn = {run_x, run_y};  cn = {'x','y','z'};
    for r = 1:2
        for c = 1:ncomp
            spec(end+1,:) = {rn{r}, ['V' cn{c}], dname(['U' cn{c}],save_disp_XYZ), 6+(r-1)*3+(c-1)};
        end
    end
end
nspec = size(spec,1);
nfile_per_pair = sum(~cellfun(@isempty,spec(:,2))) + sum(~cellfun(@isempty,spec(:,3)));
nfile_expected = nfile_per_pair*(nrecv-1);      % self pair is skipped
fprintf('%d SAC files per station pair, %d per virtual source\n', nfile_per_pair, nfile_expected);

%% ===================== pre-pass: check the sources of both runs =====================
% (serial, fails early if a run folder is missing or set up incorrectly)
scale_x = ones(nrecv,1);  scale_y = ones(nrecv,1);
for n = 1:nrecv
    stn = [sim_dir '/' char(ntwk(n)) '.' char(stnm(n))];
    [cs1,ls1,f01,fv1] = read_anchor_force([stn '/' run_x '/SeisSource.conf']);
    [cs2,ls2,f02,fv2] = read_anchor_force([stn '/' run_y '/SeisSource.conf']);
    if abs(cs1-cs2) > 1e-6 || abs(ls1-ls2) > 1e-6
        error('%s: fx and fy sources are at different locations.', stn);
    end
    if any(fv1(2:3) ~= 0) || fv1(1) == 0 || any(fv2([1 3]) ~= 0) || fv2(2) == 0
        error('%s: expected force directions (a,0,0) in %s and (0,b,0) in %s; found [%g %g %g] and [%g %g %g].', ...
              stn, run_x, run_y, fv1, fv2);
    end
    if abs(cs1-(90-slat(n))) > 0.05 || abs(mod(ls1-slon(n)+180,360)-180) > 0.05
        warning('%s: force location (%.3f,%.3f) differs from the station (colat %.3f, lon %.3f).', ...
                stn, cs1, ls1, 90-slat(n), slon(n));
    end
    % responses are scaled to a force of magnitude f0 of the fx run along each axis
    scale_x(n) = f01/(f01*fv1(1));
    scale_y(n) = f01/(f02*fv2(2));
end

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
    dirs = {[sim_stn '/' run_x], [sim_stn '/' run_y]};
    comps = {'Vx','Vy','Vz'};

    % -------------------- load coord --------------------------
    snapinfo = cell(1,2);  pnm_out = cell(1,2);
    for r = 1:2
        fnm_conf = [dirs{r} '/SeisFD3D.conf'];
        pnm_out{r} = [dirs{r} '/output/'];
        snapinfo{r} = locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
    end

    %-- get coord data
    [x,y,z]=gather_coord(snapinfo{1},'coorddir',pnm_metric);
    xs = squeeze(x);  ys = squeeze(y);
    nx = size(xs,1);  ny = size(xs,2);
    xv = xs(:,1);  yv = ys(1,:)';
    if max(max(abs(xs - xv))) > 1e-9 || max(max(abs(ys - yv'))) > 1e-9 || any(diff(xv) <= 0) || any(diff(yv) <= 0)
        error('Surface grid is not rectilinear/increasing; use interp2 as in snap2sac_spher.m.');
    end

    % bilinear sampling weights at the receivers (same result as interp2 'linear', NaN outside)
    fi = interp1(xv,(1:nx)',xrecv);
    fj = interp1(yv,(1:ny)',yrecv);
    ok = ~isnan(fi) & ~isnan(fj);
    i0 = min(floor(fi(ok)),nx-1);   j0 = min(floor(fj(ok)),ny-1);
    wx = fi(ok) - i0;               wy = fj(ok) - j0;
    l00 = i0 + (j0-1)*nx;  l10 = l00 + 1;  l01 = l00 + nx;  l11 = l00 + nx + 1;
    w00 = (1-wx).*(1-wy);  w10 = wx.*(1-wy);  w01 = (1-wx).*wy;  w11 = wx.*wy;

    % --------------------- load data --------------------------
    % A(layer, receiver, component, run); T(layer, run)
    A = nan(nlayers,nrecv,ncomp,2,'single');
    T = nan(nlayers,2);
    tload = tic;
    if fast_read && n1 == 1 && dn == 1
        % ---- block read: several time layers per file open ----
        for r = 1:2
            tcnt = snapinfo{r}(1).tcnt;               % time layers stored per NetCDF file
            if r == 1, fprintf('  %s: tcnt = %d layers per file (from SeisFD3D.conf), tblk = %d, block = %d layers\n', event_name, tcnt, tblk, min(tblk,tcnt)); end
            L = 1;
            while L <= nlayers
                nofnc  = fix((L-1)/tcnt) + 1;
                ninnc0 = mod(L-1,tcnt);
                nt = min([tblk, tcnt-ninnc0, nlayers-L+1]);
                [B,tv] = read_snap_block(snapinfo{r},id,nofnc,ninnc0,nt,comps(1:ncomp),pnm_out{r},nx,ny);
                if L == 1                              % check against the original reader
                    [Vchk,tchk] = gather_snap(snapinfo{r},id,1,comps{1},'outdir',pnm_out{r});
                    Vchk = squeeze(Vchk);
                    Vnew = reshape(B(1,:,:,1),nx,ny);
                    if max(abs(Vchk(:)-Vnew(:))) > 1e-6*max(abs(Vchk(:))) || abs(tchk-tv(1)) > 1e-9
                        error('%s: block read differs from gather_snap (run %d); set fast_read=false.', event_name, r);
                    end
                end
                for jj = 1:nt
                    for c = 1:ncomp
                        V2 = reshape(B(jj,:,:,c),nx,ny);
                        vq = nan(nrecv,1);
                        vq(ok) = w00.*V2(l00) + w10.*V2(l10) + w01.*V2(l01) + w11.*V2(l11);
                        A(L+jj-1,:,c,r) = single(vq');
                    end
                    T(L+jj-1,r) = tv(jj);
                end
                L = L + nt;
                if r == 1 || L > nlayers
                    fprintf('  %s: run %d, %d/%d layers read (%.0f s)\n', event_name, r, L-1, nlayers, toc(tload));
                end
            end
        end
    else
        for j=1:nlayers
            ly=layers(j);
            for r = 1:2
                t = NaN;                      % set inside the parfor iteration (silences the temporary-variable warning)
                for c = 1:ncomp
                    [Vsnap,t]=gather_snap(snapinfo{r},id,ly,comps{c},'outdir',pnm_out{r});
                    V2 = squeeze(Vsnap);
                    vq = nan(nrecv,1);
                    vq(ok) = w00.*V2(l00) + w10.*V2(l10) + w01.*V2(l01) + w11.*V2(l11);
                    A(j,:,c,r) = single(vq');
                end
                T(j,r) = t;
            end
            if mod(j,500) == 0 || j == nlayers
                fprintf('  %s: %d/%d layers read (%.0f s)\n', event_name, j, nlayers, toc(tload));
            end
        end  % end of nlayer
    end
    if max(abs(T(:,1)-T(:,2))) > 1e-6
        error('%s: fx and fy snapshots have different time axes.', event_name);
    end

    tt=[0; T(:,1)];dt=tt(2)-tt(1);

    % -------------------- geometry --------------------------
    alpha = sph_az(slat(nstsrc),slon(nstsrc),slat,slon);   % source -> receiver
    beta  = sph_az(slat,slon,slat(nstsrc),slon(nstsrc));   % receiver -> source
    sc = [scale_x(nstsrc) scale_y(nstsrc)];

    fidaz = fopen([OUTPUT_ROOT '/' event_name '.azimuth.txt'],'w');
    fprintf(fidaz,'# receiver  azimuth_src_to_rec(deg)  back_azimuth(deg)  dist(deg)\n');

    for nstrec=1:nrecv
        if nstrec == nstsrc, continue; end        % no azimuth for the self pair
        recvname = [ char(ntwk(nstrec)) '.' char(stnm(nstrec)) ];
        a = alpha(nstrec);  b = beta(nstrec);
        [disttemp,~]=distance(slat(nstsrc),slon(nstsrc),slat(nstrec),slon(nstrec));
        fprintf(fidaz,'%s %10.4f %10.4f %10.4f\n', recvname, a, b, disttemp);

        % responses of the two runs at this receiver (per unit-amplitude force)
        u1 = zeros(nlayers,3);  u2 = zeros(nlayers,3);      % run x, run y ; columns x,y,z
        for c = 1:ncomp
            u1(:,c) = sc(1)*double(A(:,nstrec,c,1));
            u2(:,c) = sc(2)*double(A(:,nstrec,c,2));
        end

        series = cell(11,1);
        if do_rotate
            ts = [ sind(a)  cosd(a)];   rs = [-cosd(a)  sind(a)];     % force directions at the source
            tr = [-sind(b) -cosd(b)];   rr = [ cosd(b) -sind(b)];     % T and R at the receiver
            uxT = ts(1)*u1(:,1) + ts(2)*u2(:,1);                      % response to a force along T_s
            uyT = ts(1)*u1(:,2) + ts(2)*u2(:,2);
            series{1} = tr(1)*uxT + tr(2)*uyT;                        % TT
            series{3} = uxT;  series{4} = uyT;
            if needZ, series{5} = ts(1)*u1(:,3) + ts(2)*u2(:,3); end
            if save_RR
                uxR = rs(1)*u1(:,1) + rs(2)*u2(:,1);                  % response to a force along R_s
                uyR = rs(1)*u1(:,2) + rs(2)*u2(:,2);
                series{2} = rr(1)*uxR + rr(2)*uyR;                    % RR
            end
        end
        if save_raw
            for c = 1:ncomp
                series{6+(c-1)}   = u1(:,c);
                series{6+3+(c-1)} = u2(:,c);
            end
        end

        for q = 1:nspec
            fcomp = spec{q,1};  vname = spec{q,2};  uname = spec{q,3};
            v = [0; series{spec{q,4}}]';                              % row vector, as in the original script
            S = bsac(tt,v);
            S = ch(S,'STLA',slat(nstrec),'STLO',slon(nstrec),'EVLA',slat(nstsrc),'EVLO',slon(nstsrc),...
                'DIST',disttemp);
            try, S = ch(S,'AZ',a,'BAZ',b); catch, end
            fbase = [OUTPUT_ROOT '/' event_name '.to.' recvname '.' fcomp '.'];
            if ~isempty(vname), wsac([fbase vname '.SAC'],S); end     % velocity
            if ~isempty(uname)                                         % displacement
                u = cumtrapz(v)*dt;
                S = bsac(tt,u);
                S = ch(S,'STLA',slat(nstrec),'STLO',slon(nstrec),'EVLA',slat(nstsrc),'EVLO',slon(nstsrc),...
                    'DIST',disttemp);
                try, S = ch(S,'AZ',a,'BAZ',b); catch, end
                wsac([fbase uname '.SAC'],S);
            end
        end  % end of output spec
    end  % end of receiver
    fclose(fidaz);

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

function s = dname(name, flag)
if flag, s = name; else, s = ''; end
end

function [colat,lon,f0,fvec] = read_anchor_force(fnm)
% First force line after <anchor_force>: x y z start f0 fx fy fz  (x = colatitude, y = longitude, degrees)
fid = fopen(fnm,'r');
if fid < 0, error('Cannot open %s', fnm); end
found = false;  colat = NaN;  lon = NaN;  f0 = NaN;  fvec = [NaN NaN NaN];
while true
    s = fgetl(fid);
    if ~ischar(s), break; end
    s = strtrim(s);
    if ~found
        found = strncmp(s,'<anchor_force>',14);
        continue
    end
    ic = find(s == '#', 1);  if ~isempty(ic), s = strtrim(s(1:ic-1)); end
    if isempty(s), continue; end
    v = sscanf(s,'%f');
    if numel(v) < 8, error('Unexpected force line in %s: %s', fnm, s); end
    colat = v(1);  lon = v(2);  f0 = v(5);  fvec = v(6:8)';
    break
end
fclose(fid);
if isnan(colat), error('No force line found after <anchor_force> in %s', fnm); end
end
