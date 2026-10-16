function stats = plt_data_dist_multiple_test(delayfilelist, fband, subplotpar, maxdt, maxy, cutoffmethod, labels, opts)
%PLT_DATA_DIST_MULTIPLE  Plot delay-vs-distance for several datasets, one panel per frequency band.
%
% USAGE
%   stats = plt_data_dist_multiple(delayfilelist, fband)
%   stats = plt_data_dist_multiple(delayfilelist, fband, subplotpar, maxdt, maxy, cutoffmethod, labels, opts)
%
% INPUTS
%   delayfilelist : cell array of delay file names (a single char also works).
%   fband         : [nfb x 2] frequency bands [fmin fmax] in Hz.
%   subplotpar    : [nrow ncol] panel layout.            (default: auto from nfb)
%   maxdt         : delay threshold (s) for 'flat'.      (default: 7.5)
%   maxy          : y-axis limit (s), axis is +-maxy.    (default: 1.5*maxdt)
%   cutoffmethod  : 'no' (plot only), 'flat', 'norm', 'nfalc', 'draw', 'dv'. (default: 'no')
%                   flat  : single threshold +-maxdt for all bands.
%                   norm  : mean +- nsigma*std of the delays in each band.
%                   nfalc : normal fit after linear (delay vs distance) correction.
%                   draw  : click two points for the lower and two for the upper bounding line.
%                   dv    : type dv_max for each band; lines are +-0.95*dv_max/cmax.
%   labels        : cell array / string array of dataset names for the legend,
%                   one per file. If empty or omitted: data1, data2, ...
%   opts          : optional struct, any of the fields below (defaults in brackets)
%                   .nsigma [1.96]    confidence multiplier for norm / nfalc
%                   .cmax [3]         used by 'dv'
%                   .colors []        [n x 3] RGB; default is generated automatically
%                   .markersize [8]
%                   .fontsize [11]
%                   .shownum [true]   print N (points per band) in each panel
%                   .save [true]      write the figure to file (see .format)
%                   .format ['png']   figure format: 'png','pdf','eps','svg','jpg','tif'
%                                     (also 'epsc'/'jpeg'/'tiff'); sets the file extension
%                                     and the matching print device (-dpng, -dpdf, ...).
%                   .outfile []       file name; default = [delayfilelist{end} '_plot.' format].
%                                     If it has no extension, '.format' is appended.
%                   .distunit ['km']  distance unit on the x axis: 'km', 'm', 'mi' or 'deg'
%                                     (deg uses 111.195 km/deg). Data are always read in km;
%                                     only the plotted values and axis label change, and the
%                                     drop lists / logs stay in km.
%                   .bandunit ['s']   unit of the panel titles: 's' (period, default),
%                                     'Hz' or 'mHz' (frequency). 'seconds' is also accepted.
%
% OUTPUT
%   stats : struct array (one element per dataset) with fields
%           label, file, n (nfb x 1), rms (nfb x 1), nremoved (nfb x 1).
%
% Cutoffs are only written to drop-list text files (see exclude_measurements.sh);
% no data are removed here.
%
% Modified from plt_data_dist.m by Haiying Gao.
% Xiaotao Yang @ UMASS, 12/13/2016, 1/24/2017.
% Revised: dataset labels as input, tiled layout, automatic colors,
% unified cutoff handling, optional outputs.
%
% Requires MATLAB R2019b+ (tiledlayout) and the Mapping Toolbox (distance).

%% ---------- inputs and defaults ----------
if nargin < 2 || isempty(fband)
    error('fband (nfb x 2 array of frequency bands) is required.');
end
if ischar(delayfilelist) || isstring(delayfilelist)
    delayfilelist = cellstr(delayfilelist);
end
numinput = numel(delayfilelist);
nfb   = size(fband, 1);
pband = flip(1./fband, 2);                      % period bands [Tmin Tmax]

if nargin < 3 || isempty(subplotpar)
    nr = floor(sqrt(nfb));  subplotpar = [nr, ceil(nfb/nr)];
end
if nfb > prod(subplotpar)
    error('subplotpar [%d %d] has fewer panels than the %d frequency bands.', subplotpar, nfb);
end
if nargin < 4 || isempty(maxdt),  maxdt = 7.5;         end
if nargin < 5 || isempty(maxy),   maxy  = 1.5*maxdt;   end
if nargin < 6 || isempty(cutoffmethod), cutoffmethod = 'no'; end
cutoffmethod = validatestring(cutoffmethod, {'no','flat','norm','nfalc','draw','dv'});

% dataset labels
if nargin < 7 || isempty(labels)
    labels = compose('data%d', 1:numinput);
else
    labels = cellstr(labels);  labels = labels(:)';
    if numel(labels) ~= numinput
        warning('%d labels given for %d datasets; missing ones are filled with dataN.', numel(labels), numinput);
        labels(end+1:numinput) = compose('data%d', numel(labels)+1:numinput);
        labels = labels(1:numinput);
    end
end

% options
o = struct('nsigma',1.96, 'cmax',3, 'colors',[], 'markersize',8, 'fontsize',11, ...
           'shownum',true, 'save',true, 'outfile',[], ...
           'distunit','km', 'bandunit','s', 'format','png');
if nargin >= 8 && ~isempty(opts)
    fn = fieldnames(opts);
    for n = 1:numel(fn), o.(fn{n}) = opts.(fn{n}); end
end
if isempty(o.colors), o.colors = default_colors(numinput); end
o.colors = o.colors(mod(0:numinput-1, size(o.colors,1)) + 1, :);   % cycle if too few

% units for the x axis and the panel titles
[dscale, dlabel] = dist_unit(o.distunit);           % km -> plotted unit
[bv, bu] = band_values(o.bandunit, fband, pband);   % band values shown in titles
bandstr = cell(nfb, 1);
for i = 1:nfb, bandstr{i} = sprintf('%.4g-%.4g %s', bv(i,1), bv(i,2), bu); end

docutoff    = ~strcmp(cutoffmethod, 'no');
interactive = any(strcmp(cutoffmethod, {'draw','dv'}));

%% ---------- figure and axes ----------
nr = subplotpar(1);  nc = subplotpar(2);
fig = figure('Color','w', 'Position',[100 100 250*nc+80 200*nr+110]);
tl  = tiledlayout(fig, nr, nc, 'TileSpacing','compact', 'Padding','compact');
ax  = gobjects(nfb, 1);
for i = 1:nfb
    ax(i) = nexttile(tl);
    hold(ax(i), 'on');  box(ax(i), 'on');
    ylim(ax(i), [-maxy maxy]);
end

maxdist = zeros(nfb, 1);
mindist = inf(nfb, 1);
stats   = repmat(struct('label','','file','','n',[],'rms',[],'nremoved',[]), numinput, 1);

%% ---------- loop over datasets ----------
for k = 1:numinput
    tempdelayfile = delayfilelist{k};
    col = o.colors(k, :);
    fprintf('\n[%s] %s\n', labels{k}, tempdelayfile);
    bands = read_delay_file(tempdelayfile, nfb);

    stats(k).label = labels{k};  stats(k).file = tempdelayfile;
    stats(k).n = zeros(nfb,1);   stats(k).rms = nan(nfb,1);
    pl_all = nan(2, nfb);  pu_all = nan(2, nfb);       % [slope; intercept] of drawn / dv lines

    for i = 1:nfb
        b = bands(i);
        n = numel(b.delay);
        stats(k).n(i) = n;
        if n == 0, continue; end
        stats(k).rms(i) = sqrt(mean(b.delay.^2));

        % bounds that need no interaction
        if docutoff && ~interactive
            b = auto_bounds(b, cutoffmethod, maxdt, o.nsigma);
        end

        plot(ax(i), b.dist*dscale, b.delay, '.', 'Color', col, 'MarkerSize', o.markersize, ...
             'HandleVisibility', 'off');
        maxdist(i) = max(maxdist(i), max(b.dist)*dscale);
        mindist(i) = min(mindist(i), min(b.dist)*dscale);

        % bounds that need interaction (data are already on screen)
        if interactive
            drawnow;
            fprintf('Band %d (%s):\n', i, bandstr{i});
            if strcmp(cutoffmethod, 'draw')
                [b, pl_all(:,i), pu_all(:,i)] = pick_bounds(ax(i), b, col, dscale);
            else
                dv_max = input(sprintf('--> Type dv max (0-1) for f%d: ', i));
                pl = [-0.95*dv_max/o.cmax, 0];  pu = [0.95*dv_max/o.cmax, 0];
                b.lo = polyval(pl, b.dist);     b.hi = polyval(pu, b.dist);
                pl_all(:,i) = pl(:);  pu_all(:,i) = pu(:);
            end
        end

        % draw the cutoff lines in the dataset color (dashed)
        if docutoff
            [xs, ix] = sort(b.dist);
            plot(ax(i), xs*dscale, b.lo(ix), '--', 'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
            plot(ax(i), xs*dscale, b.hi(ix), '--', 'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        end
        bands(i) = b;

        if o.shownum
            text(ax(i), 0.98, 0.97 - 0.085*(k-1), sprintf('N=%d', n), 'Units','normalized', ...
                 'HorizontalAlignment','right', 'VerticalAlignment','top', 'Color', col, ...
                 'FontSize', o.fontsize-2, 'FontWeight','bold');
        end
    end

    % cutoff procedures: write drop list + log, no data removed
    if docutoff
        [cleanupscript, logpar] = droplist_names(tempdelayfile, cutoffmethod, maxdt, o.nsigma);
        fprintf('Generating clean-up list: [ %s ] ...\n', cleanupscript);
        stats(k).nremoved = write_droplist(cleanupscript, cutoffmethod, bands, o.cmax);
        write_log([tempdelayfile '_cleanup_log.txt'], tempdelayfile, cutoffmethod, logpar, pl_all, pu_all);
        fprintf('---> Need to remove [ %d ] data points\n', sum(stats(k).nremoved));
    end
end

%% ---------- axes decoration ----------
for i = 1:nfb
    a = ax(i);
    yline(a, 0, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.75, 'HandleVisibility', 'off');
    if maxdist(i) > mindist(i)                      % x range = data range over all datasets
        xlim(a, [mindist(i) maxdist(i)]);
    elseif maxdist(i) == 0
        text(a, 0.5, 0.5, 'no data', 'Units','normalized', 'HorizontalAlignment','center', ...
             'Color',[0.5 0.5 0.5], 'FontSize', o.fontsize);
    end
    ylim(a, [-maxy maxy]);
    title(a, sprintf('\\bf%s.\\rm  %s', char('a'+i-1), bandstr{i}), ...
          'FontSize', o.fontsize+1, 'FontWeight','normal');
    set(a, 'TitleHorizontalAlignment','left', 'FontSize', o.fontsize, 'TickDir','out', 'Layer','top', 'LineWidth', 0.75, ...
           'XGrid','on', 'YGrid','on', 'GridAlpha', 0.15);
end
xlabel(tl, sprintf('Distance (%s)', dlabel), 'FontSize', o.fontsize+1);
ylabel(tl, 'Delay (s)',     'FontSize', o.fontsize+1);

% one shared legend below the panels (dummy NaN handles so legend markers are large)
hd = gobjects(numinput, 1);
for k = 1:numinput
    hd(k) = plot(ax(1), NaN, NaN, '.', 'Color', o.colors(k,:), 'MarkerSize', o.markersize*2.5, ...
                 'DisplayName', labels{k});
end
lgd = legend(ax(1), hd, 'Orientation','vertical', 'Box','off', 'FontSize', o.fontsize+1);
lgd.Layout.Tile = 'east';

%% ---------- save ----------
if o.save
    [ext, dev, isvec] = figure_format(o.format);
    outfile = o.outfile;
    if isempty(outfile)
        outfile = [delayfilelist{end} '_plot.' ext];
    else
        [~, ~, e] = fileparts(outfile);
        if isempty(e), outfile = [outfile '.' ext]; end
    end
    set(fig, 'PaperPositionMode', 'auto');
    if isvec
        print(fig, outfile, dev, '-vector');
    else
        print(fig, outfile, dev, '-r300');
    end
    fprintf('\nSaved figure: %s\n', outfile);
end

print_summary(stats, bv, bu);

end


%% =====================================================================
%  local functions
%% =====================================================================
function bands = read_delay_file(fname, nfb)
% Read one delay file and split it by frequency tag (f1, f2, ...).
fid = fopen(fname, 'r');
if fid < 0, error('Cannot open delay file: %s', fname); end
cl = onCleanup(@() fclose(fid));
d = textscan(fid, '%s %s %f %f %f %f %f %*f %s %*f %*f');
[src, rcv, slat, slon, rlat, rlon, delay, fb] = d{1:8};
dist = distance(slat, slon, rlat, rlon, [6378.14 0.0818]);      % km, ellipsoidal

bands = repmat(struct('delay',[], 'dist',[], 'src',{{}}, 'rcv',{{}}, 'lo',[], 'hi',[]), nfb, 1);
for i = 1:nfb
    idx = strcmp(sprintf('f%d', i), fb);
    bands(i).delay = delay(idx);
    bands(i).dist  = dist(idx);
    bands(i).src   = src(idx);
    bands(i).rcv   = rcv(idx);
    bands(i).lo    = nan(nnz(idx), 1);
    bands(i).hi    = nan(nnz(idx), 1);
end
end


function b = auto_bounds(b, method, delaycutoff, nsigma)
% Lower/upper cutoff for every measurement (non-interactive methods).
n = numel(b.delay);
switch method
    case 'flat'
        b.lo = -delaycutoff*ones(n,1);   b.hi = delaycutoff*ones(n,1);
    case 'norm'
        mu = mean(b.delay);  s = std(b.delay);
        b.lo = (mu - nsigma*s)*ones(n,1);  b.hi = (mu + nsigma*s)*ones(n,1);
    case 'nfalc'
        if n < 2, return; end
        p   = polyfit(b.dist, b.delay, 1);
        fit = polyval(p, b.dist);
        res = b.delay - fit;
        mu  = mean(res);  s = std(res);
        b.lo = fit + (mu - nsigma*s);    b.hi = fit + (mu + nsigma*s);
end
end


function [b, pl, pu] = pick_bounds(a, b, col, dscale)
% Click two points for the lower line and two for the upper line.
% Clicks are in plotted units; fits are kept per km so drop lists/logs are unchanged.
fig = ancestor(a, 'figure');
fig.CurrentAxes = a;
pl = pick_line(a, 'lower', col, dscale);
b.lo = polyval(pl, b.dist);
plot(a, sort(b.dist)*dscale, polyval(pl, sort(b.dist)), '-', 'Color', col, 'HandleVisibility','off');
pu = pick_line(a, 'upper', col, dscale);
b.hi = polyval(pu, b.dist);
pl = pl(:);  pu = pu(:);
end


function p = pick_line(a, name, col, dscale)
fprintf('  choose two points for the %s bound line:\n', name);
xy = zeros(2, 2);
for m = 1:2
    [x, y] = ginput(1);
    xy(m,:) = [x/dscale y];                 % back to km
    plot(a, x, y, 'o', 'Color', col, 'MarkerFaceColor', col, 'HandleVisibility','off');
end
p = polyfit(xy(:,1), xy(:,2), 1);
end


function [fname, par] = droplist_names(tempfile, method, delaycutoff, nsigma)
switch method
    case 'flat',  fname = [tempfile '_droplist_flat' num2str(delaycutoff) '.txt']; par = delaycutoff;
    case 'norm',  fname = [tempfile '_droplist_norm' num2str(nsigma) '.txt'];      par = nsigma;
    case 'nfalc', fname = [tempfile '_droplist_nfalc' num2str(nsigma) '.txt'];     par = nsigma;
    case 'draw',  fname = [tempfile '_droplist_draw.txt'];                         par = [];
    case 'dv',    fname = [tempfile '_droplist_dv.txt'];                           par = [];
end
end


function ncut = write_droplist(fname, method, bands, cmax)
% One line per measurement outside [lo, hi]: "src_rcv  f<band>"
fid = fopen(fname, 'w');
if fid < 0, error('Cannot write %s', fname); end
cl = onCleanup(@() fclose(fid));
if strcmp(method, 'dv'), fprintf(fid, 'cmax: %g\n', cmax); end
ncut = zeros(numel(bands), 1);
for i = 1:numel(bands)
    b   = bands(i);
    idx = find(b.delay > b.hi | b.delay < b.lo);
    ncut(i) = numel(idx);
    for j = idx(:)'
        fprintf(fid, '%s_%s  f%d\n', b.src{j}, b.rcv{j}, i);
    end
end
end


function write_log(fname, tempfile, method, par, pl_all, pu_all)
fid = fopen(fname, 'w');
if fid < 0, error('Cannot write %s', fname); end
cl = onCleanup(@() fclose(fid));
fprintf(fid, 'original file: %s\n', tempfile);
fprintf(fid, 'method: %s\n', method);
if any(strcmp(method, {'draw','dv'}))
    for i = 1:size(pl_all, 2)
        if any(isnan([pl_all(:,i); pu_all(:,i)])), continue; end
        fprintf(fid, 'parameter f%d: lower line [%f %f], upper line [%f %f]\n', ...
                i, pl_all(1,i), pl_all(2,i), pu_all(1,i), pu_all(2,i));
    end
else
    fprintf(fid, 'parameter: %f\n', par);
end
end


function c = default_colors(n)
% First color is neutral gray (reference), then a colorblind-friendly set;
% beyond 8 datasets fall back to a perceptual colormap.
base = [0.35 0.35 0.35;   % gray
        0.84 0.19 0.15;   % red
        0.00 0.45 0.70;   % blue
        0.00 0.62 0.45;   % green
        0.90 0.62 0.00;   % orange
        0.80 0.47 0.65;   % pink
        0.34 0.71 0.91;   % sky
        0.55 0.34 0.29];  % brown
if n <= size(base, 1)
    c = base(1:n, :);
elseif exist('turbo', 'file')
    c = turbo(n);
else
    c = parula(n);
end
end


function print_summary(stats, pband, bu)
fprintf('\n%-10s', sprintf('band (%s)', bu));
for k = 1:numel(stats), fprintf(' | %-24s', stats(k).label); end
fprintf('\n%-10s', '');
for k = 1:numel(stats), fprintf(' | %8s %8s %6s', 'N', 'RMS(s)', 'cut'); end
fprintf('\n');
for i = 1:size(pband, 1)
    fprintf('%-10s', sprintf('%.4g-%.4g', pband(i,1), pband(i,2)));
    for k = 1:numel(stats)
        nrm = 0;
        if ~isempty(stats(k).nremoved), nrm = stats(k).nremoved(i); end
        fprintf(' | %8d %8.3f %6d', stats(k).n(i), stats(k).rms(i), nrm);
    end
    fprintf('\n');
end
end

function [scale, label] = dist_unit(u)
% Conversion factor from km (data unit) to the requested plotting unit.
switch lower(char(u))
    case 'km',                     scale = 1;         label = 'km';
    case {'m','meter','meters'},   scale = 1e3;       label = 'm';
    case {'mi','mile','miles'},    scale = 0.621371;  label = 'mi';
    case {'deg','degree','degrees'}, scale = 1/111.195; label = 'deg';
    otherwise, error('opts.distunit must be ''km'', ''m'', ''mi'' or ''deg'' (got ''%s'').', char(u));
end
end


function [bv, bu] = band_values(u, fband, pband)
% Band values and unit string for the panel titles.
switch lower(char(u))
    case {'s','sec','second','seconds'}, bv = pband;      bu = 's';
    case 'hz',                           bv = fband;      bu = 'Hz';
    case 'mhz',                          bv = fband*1e3;  bu = 'mHz';
    otherwise, error('opts.bandunit must be ''s'', ''Hz'' or ''mHz'' (got ''%s'').', char(u));
end
end


function [ext, dev, isvec] = figure_format(f)
% File extension, print device and vector flag for opts.format.
switch lower(regexprep(char(f), '^\.', ''))
    case 'png',          ext = 'png'; dev = '-dpng';   isvec = false;
    case {'jpg','jpeg'}, ext = 'jpg'; dev = '-djpeg';  isvec = false;
    case {'tif','tiff'}, ext = 'tif'; dev = '-dtiff';  isvec = false;
    case 'pdf',          ext = 'pdf'; dev = '-dpdf';   isvec = true;
    case {'eps','epsc'}, ext = 'eps'; dev = '-depsc';  isvec = true;
    case 'svg',          ext = 'svg'; dev = '-dsvg';   isvec = true;
    otherwise
        error('opts.format must be png, jpg, tif, pdf, eps or svg (got ''%s'').', char(f));
end
end
