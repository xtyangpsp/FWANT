function stats = plt_raypath(delayfile,stationfile,fband,ite_tag,ray_min,gridstep,cmax,opts)
% PLT_RAYPATH  Plot ray-path coverage (rays per cell) for each period band.
%
%   stats = plt_raypath(delayfile,stationfile,fband,ite_tag,ray_min,gridstep)
%   stats = plt_raypath(...,gridstep,cmax,opts)
%
% delayfile   : phase delay file (output of get_delay_latlon()).
% stationfile : text file with station info (cols: name net lon lat elev).
% fband       : frequency bands [f1min f1max; f2min f2max; ...]
% ite_tag     : iteration tag (used in the figure title and outline file name).
% ray_min     : minimum ray count for a cell to be drawn / counted as covered.
% gridstep    : cell size in deg. Scalar = same size for lat and lon (square cells);
%               [dlat dlon] = separate latitude and longitude cell sizes.
% cmax        : (optional) maximum colorbar value. [] = auto (max over all bands).
% opts        : (optional) struct, any of
%                 .saveOutline  save boundary of the covered area (.mat)  (false)
%                 .nrows        number of panel rows; columns follow      (2)
%                 .margin       map margin around the stations, deg       (0.05)
%                 .savefig      write the PNG                             (true)
%                 .cmap         colormap, Nx3                             (parula)
%
% A ray is counted once in every cell it crosses.
%
% stats       : struct array, one element per band, with fields
%               band, period, nrays, maxcount, ncells, status
%               status = 'ok' | 'no data' | 'below ray_min'
%
% Bands with no rays, or with no cell reaching ray_min, are drawn as empty
% panels with an explanatory label instead of failing.

if nargin < 6 || ~isnumeric(gridstep) || ~any(numel(gridstep)==[1 2]) || any(~(gridstep>0))
    error('plt_raypath:gridstep','gridstep (positive scalar, or [dlat dlon], in deg) is required.');
end
if nargin < 7, cmax = []; end
if nargin < 8 || isempty(opts), opts = struct(); end
if isempty(ite_tag), ite_tag = ''; end
def = struct('saveOutline',false,'nrows',2,'margin',0.05,'savefig',true,'cmap',parula(256));
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts,fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end
end

wrap180 = @(x) mod(x+180,360)-180;   % accept either 0-360 or -180-180 input

%% ---- read inputs ---------------------------------------------------------
fid = fopen(stationfile,'r');
if fid < 0, error('plt_raypath:file','Cannot open station file: %s',stationfile); end
C = textscan(fid,'%*s %*s %f %f %*f'); fclose(fid);
stationlon = wrap180(C{1}); stationlat = C{2};

fid = fopen(delayfile,'r');
if fid < 0, error('plt_raypath:file','Cannot open delay file: %s',delayfile); end
C = textscan(fid,'%*s %*s %f %f %f %f %*f %*f %s'); fclose(fid);
[slat,slon,rlat,rlon,fb] = C{1:5};
slon = wrap180(slon); rlon = wrap180(rlon);
if isempty(slat)
    warning('plt_raypath:empty','Delay file contains no measurements: %s',delayfile);
end

%% ---- map area and grid ---------------------------------------------------
if isscalar(gridstep), dlat = gridstep; dlon = gridstep;   % same size in lat and lon
else,                  dlat = gridstep(1); dlon = gridstep(2);
end
maparea.lon = [min(stationlon) max(stationlon)] + opts.margin*[-1 1];
maparea.lat = [min(stationlat) max(stationlat)] + opts.margin*[-1 1];

lonEdges = min(stationlon) : dlon : min(stationlon) + ceil((max(stationlon)-min(stationlon))/dlon)*dlon;
latEdges = min(stationlat) : dlat : min(stationlat) + ceil((max(stationlat)-min(stationlat))/dlat)*dlat;
nlon = numel(lonEdges)-1;  nlat = numel(latEdges)-1;
if nlon < 1 || nlat < 1
    error('plt_raypath:grid','gridstep (dlat=%g, dlon=%g) is too large for the station extent.',dlat,dlon);
end
lonc = lonEdges(1:end-1) + dlon/2;           % cell centres
latc = latEdges(1:end-1) + dlat/2;

state = [];
if exist('us_states.mat','file')
    S = load('us_states.mat');
    if isfield(S,'Alaska')
        state(1).polygon = [S.Alaska.lon(:) S.Alaska.lat(:)];
    end
else
    warning('plt_raypath:nostates','us_states.mat not found; state outlines skipped.');
end

%% ---- pass 1: ray counts for every band ----------------------------------
nb = size(fband,1);
cc = cell(nb,1);
stats = struct('band',num2cell((1:nb)'),'period',[],'nrays',0,'maxcount',0,'ncells',0,'status','');
for k = 1:nb
    stats(k).period = [1/fband(k,2) 1/fband(k,1)];
    fm = find(strcmp(fb, ['f' num2str(k)]));
    x1 = slon(fm); y1 = slat(fm); x2 = rlon(fm); y2 = rlat(fm);
    good = isfinite(x1) & isfinite(y1) & isfinite(x2) & isfinite(y2);
    x1 = x1(good); y1 = y1(good); x2 = x2(good); y2 = y2(good);
    stats(k).nrays = numel(x1);

    c = zeros(nlat,nlon);
    for idx = 1:numel(x1)
        dx = x2(idx)-x1(idx);  dy = y2(idx)-y1(idx);
        n  = max(2, ceil(max(abs(dx)/(dlon/2), abs(dy)/(dlat/2))) + 1);   % >= 2 samples per cell
        t  = linspace(0,1,n);
        ix = floor((x1(idx)+t*dx - lonEdges(1))/dlon) + 1;
        iy = floor((y1(idx)+t*dy - latEdges(1))/dlat) + 1;
        ok = ix>=1 & ix<=nlon & iy>=1 & iy<=nlat;
        lin = unique(sub2ind([nlat nlon], iy(ok), ix(ok)));   % once per cell
        c(lin) = c(lin) + 1;
    end
    cc{k} = c;
    stats(k).maxcount = max(c(:));
    stats(k).ncells   = nnz(c >= ray_min);
    if stats(k).nrays == 0
        stats(k).status = 'no data';
    elseif stats(k).ncells == 0
        stats(k).status = 'below ray_min';
    else
        stats(k).status = 'ok';
    end
end

if isempty(cmax), cmax = max([stats.maxcount]); end        % common scale across bands
if cmax <= ray_min, cmax = ray_min + 1; end                % CLim must be increasing

%% ---- pass 2: plotting ---------------------------------------------------
nrows = max(1, min(round(opts.nrows), nb));  ncols = ceil(nb/nrows);   % rows first, columns follow
fig = figure('Color','w','Position',[100 100 270*ncols+120 250*nrows+90]);
tl  = tiledlayout(fig,nrows,ncols,'TileSpacing','compact','Padding','compact');
if isempty(ite_tag), ttl = sprintf('Ray path coverage (%.3g\\circ x %.3g\\circ cells)',dlon,dlat);
else,                ttl = sprintf('Ray path coverage, %s (%.3g\\circ x %.3g\\circ cells)',strrep(ite_tag,'_','\_'),dlon,dlat);
end
title(tl,ttl,'FontSize',14,'FontWeight','bold');
xlabel(tl,'Longitude (\circ)','FontSize',12);
ylabel(tl,'Latitude (\circ)','FontSize',12);

axs = gobjects(nb,1);
for k = 1:nb
    ax = nexttile(tl); axs(k) = ax;
    hold(ax,'on'); box(ax,'on');
    c = cc{k};
    has = strcmp(stats(k).status,'ok');

    set(ax,'Color',[0.94 0.94 0.94],'Layer','top','TickDir','out','FontSize',11, ...
        'YDir','normal','CLim',[ray_min cmax]);
    colormap(ax,opts.cmap);

    if has
        hi = image(ax,lonc,latc,c,'CDataMapping','scaled');
        hi.AlphaData = double(c >= ray_min);
    end
    for sb = 1:numel(state)
        plot(ax,state(sb).polygon(:,1),state(sb).polygon(:,2),'Color',[.2 .2 .2],'LineWidth',.5);
    end
    plot(ax,stationlon,stationlat,'^','MarkerSize',4,'LineWidth',0.5, ...
        'MarkerEdgeColor',[.25 .25 .25],'MarkerFaceColor',[.85 .85 .85]);

    axis(ax,[maparea.lon maparea.lat]);
    daspect(ax,[1 cosd(mean(maparea.lat)) 1]);
    title(ax,sprintf('%c. %g-%g s  (%d rays)',char('a'+k-1),stats(k).period,stats(k).nrays), ...
        'FontSize',14,'FontWeight','normal');
    set(gca, 'TitleHorizontalAlignment','left');

    % keep tick labels only on the outer edges of the tile grid
    if k + ncols <= nb, ax.XTickLabel = []; end
    if mod(k-1,ncols) ~= 0, ax.YTickLabel = []; end

    if ~has   % label empty panels so a blank map is not mistaken for a plotting bug
        if strcmp(stats(k).status,'no data')
            msg = 'No data';
        else
            msg = sprintf('No cell with\n\\geq %g rays\n(max %d)',ray_min,stats(k).maxcount);
        end
        text(ax,mean(maparea.lon),mean(maparea.lat),msg,'HorizontalAlignment','center', ...
            'FontSize',12,'Color',[.3 .3 .3],'BackgroundColor','w','EdgeColor',[.7 .7 .7],'Margin',4);
    end

    if has && opts.saveOutline
        [jj,ii] = find(c >= ray_min);
        if numel(ii) >= 3
            try
                pts = [lonc(ii(:))' latc(jj(:))'];
                kk = boundary(pts(:,1),pts(:,2),1);
                raycover.data = pts(kk,:); raycover.cutoff = ray_min; %#ok<STRNU>
                save(sprintf('RayCoverOutline_%s_%g-%gs_cutoff%g.mat',ite_tag,stats(k).period,ray_min), ...
                    'raycover');
            catch ME
                warning('plt_raypath:boundary','Outline failed for band %d: %s',k,ME.message);
            end
        end
    end
    hold(ax,'off');
end

% one shared colorbar (all panels use the same colour limits)
cb = colorbar(axs(end),'eastoutside');
cb.Layout.Tile = 'east';
cb.Label.String = sprintf('Rays per cell (shown where \\geq %g)',ray_min);
cb.TickDirection = 'out'; cb.FontSize = 11;
drawnow

%% ---- summary and output --------------------------------------------------
fprintf('%-5s %-12s %8s %9s %8s  %s\n','band','period(s)','rays','maxcount','cells','status');
for k = 1:nb
    fprintf('%-5d %5g-%-6g %8d %9d %8d  %s\n',k,stats(k).period,stats(k).nrays, ...
        stats(k).maxcount,stats(k).ncells,stats(k).status);
end

if opts.savefig
    outfile = [delayfile '_pathcoverage' num2str(ray_min) '.png'];
    if exist('exportgraphics','file')
        exportgraphics(fig,outfile,'Resolution',300);
    else
        set(fig,'PaperPositionMode','auto');
        print(fig,outfile,'-dpng','-painters');
    end
end
end
