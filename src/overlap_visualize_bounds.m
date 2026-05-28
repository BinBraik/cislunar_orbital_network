function B = overlap_visualize_bounds(V, SA, SB, O, cfg, outdir, tag)
%OVERLAP_VISUALIZE_BOUNDS  Plot voxel-wise DVtotal bounds (vectorized).
%
% Bounds (m/s) per overlap voxel:
%   LB    = min(dv_turn_A) + min(dv_turn_B)
%   proxy = LB + DVpatch_ub   (DVpatch_ub uses voxel-center speed with CJ*)
%   UB    = max(dv_turn_A) + DVpatch_ub + max(dv_turn_B)
%
% V must be the flat-array struct from overlap_extract_voxel_info.
% All computations are fully vectorized — no per-voxel loop.

if nargin < 6 || isempty(tag),    tag    = 'bounds'; end
if nargin < 5 || isempty(outdir), outdir = pwd;      end
if ~exist(outdir,'dir'), mkdir(outdir); end

if ~isstruct(V) || ~isfield(V,'x') || isempty(V.x)
    warning('[overlap] Bounds visualization skipped: empty or old-format voxel metadata.');
    B = struct('dv_lb',[],'dv_ub',[],'dv_patch_ub',[],'dv_proxy',[], ...
               'min_dvproxy',NaN,'imin',NaN,'x_at_min',NaN,'y_at_min',NaN);
    return;
end

doLB    = local_plot_enabled(cfg, 'plot.overlap.bounds_lb',    true);
doUB    = local_plot_enabled(cfg, 'plot.overlap.bounds_ub',    false);
doProxy = local_plot_enabled(cfg, 'plot.overlap.bounds_proxy', true);

grid3   = SA.grid3;
safeTag = sanitize_fname(tag);
VU_mps  = local_cfg_get(cfg, 'units.VU_mps', 1.0);
CJstar  = min(SA.CJ, SB.CJ);
Nx = numel(grid3.x_centers);
Ny = numel(grid3.y_centers);
Nt = numel(grid3.th_centers);

% ---------------- FULL sets (BRS = R(FRS), new schema from PR #25) --------
% rows_FRS_upper + rows_FRS_lower are directly stored from forward integration.
% BRS_full = R(FRS_full): BRS_upper=R(FRS_lower), BRS_lower=R(FRS_upper).
% Gracefully degrade for old-format caches missing rows_FRS_lower.
hasFullA = isfield(SA,'Step4') && isfield(SA.Step4,'rows_FRS_upper') && isfield(SA.Step4,'rows_FRS_lower');
hasFullB = isfield(SB,'Step4') && isfield(SB.Step4,'rows_FRS_upper') && isfield(SB.Step4,'rows_FRS_lower');

Ax = []; Ay = []; Bx = []; By = [];
idsO = O.ids(:);
[Ox, Oy, ~] = local_ids_to_centers(idsO, grid3, Ny, Nx, Nt);

if hasFullA && hasFullB
    % A.FRS_full = FRS_upper + FRS_lower (both directly stored)
    rowsA_F = local_rows_cat(SA.Step4.rows_FRS_upper, SA.Step4.rows_FRS_lower);

    % B.BRS_full = R(B.FRS_full): mirror both halves
    rowsB_B_u = atlas_rows_mirror_lower(SB.Step4.rows_FRS_lower, grid3, 2);
    rowsB_B_l = atlas_rows_mirror_lower(SB.Step4.rows_FRS_upper, grid3, 2);
    rowsB_B   = local_rows_cat(rowsB_B_u, rowsB_B_l);

    idsA = unique(local_rows_to_vid(rowsA_F, Ny, Nx, Nt));
    idsB = unique(local_rows_to_vid(rowsB_B, Ny, Nx, Nt));

    idsA_only = setdiff(idsA, idsO);
    idsB_only = setdiff(idsB, idsO);

    [Ax, Ay, ~] = local_ids_to_centers(idsA_only, grid3, Ny, Nx, Nt);
    [Bx, By, ~] = local_ids_to_centers(idsB_only, grid3, Ny, Nx, Nt);
else
    warning('[overlap] Atlas missing rows_FRS_lower — non-overlap background omitted from proxy plot. Rebuild cache to include it.');
end

% ---------- flat arrays from V ----------
x  = V.x(:);
y  = V.y(:);

% ---------- vectorized voxel-center potential (one call) ----------
pot    = cr3bp_potential(x, y, SA.mu);
v_box  = sqrt(max(2*pot.U - CJstar, 0));
dv_patch_ub = 2 * v_box * sin(abs(grid3.dtheta)/2) * VU_mps;

% ---------- bounds ----------
dv_lb = V.dv_turn_mps_min_A + V.dv_turn_mps_min_B;

if doUB
    dv_ub = V.dv_turn_mps_max_A + dv_patch_ub + V.dv_turn_mps_max_B;
else
    dv_ub = [];   % not needed
end

dv_proxy = dv_lb + dv_patch_ub;

% ---------- best voxel ----------
valid = isfinite(dv_proxy);
if any(valid)
    idxValid = find(valid);
    [~, iLocal] = min(dv_proxy(valid));
    iMin = idxValid(iLocal);
else
    iMin = NaN;
end

B = struct();
B.dv_lb       = dv_lb;
B.dv_ub       = dv_ub;
B.dv_patch_ub = dv_patch_ub;
B.dv_proxy    = dv_proxy;
B.min_dvproxy = double(dv_proxy(iMin));
B.imin        = iMin;
if isfinite(iMin) && iMin >= 1 && iMin <= numel(x)
    B.x_at_min = x(iMin);
    B.y_at_min = y(iMin);
else
    B.x_at_min = NaN;
    B.y_at_min = NaN;
end

% ---------- plots ----------
CJbg = min(SA.CJ, SB.CJ);

if doLB
    fig = figure('Color','w', 'Name',['DVtotal Lower Bound XY ' tag], ...
                 'Visible', local_fig_visible(cfg));
    ax = gca;
    cr3bp_plot_background(CJbg, SA.mu, ax);
    set(ax.Children,'HandleVisibility','off');
    hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
    scatter(ax, x, y, 18, dv_lb, 'filled', 'MarkerFaceAlpha', 0.92, 'MarkerEdgeAlpha', 0.25);
    cb = colorbar(ax); ylabel(cb, 'DVtotal lower bound (m/s)');
    title(ax, sprintf('Voxel DVtotal lower bound | %s', tag), 'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'y');
    local_apply_zoom(cfg, ax);
    io_save_figure(fig, outdir, ['overlap_' safeTag '_bounds_lb_xy'], cfg, 'Resolution', local_fig_res(cfg));
    local_close_if_hidden(cfg, fig);
end

if doUB && ~isempty(dv_ub)
    fig = figure('Color','w', 'Name',['DVtotal Upper Bound XY ' tag], ...
                 'Visible', local_fig_visible(cfg));
    ax = gca;
    cr3bp_plot_background(CJbg, SA.mu, ax);
    set(ax.Children,'HandleVisibility','off');
    hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
    scatter(ax, x, y, 18, dv_ub, 'filled', 'MarkerFaceAlpha', 0.92, 'MarkerEdgeAlpha', 0.25);
    cb = colorbar(ax); ylabel(cb, 'DVtotal upper bound (m/s)');
    title(ax, sprintf('Voxel DVtotal upper bound | %s', tag), 'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'y');
    local_apply_zoom(cfg, ax);
    io_save_figure(fig, outdir, ['overlap_' safeTag '_bounds_ub_xy'], cfg, 'Resolution', local_fig_res(cfg));
    local_close_if_hidden(cfg, fig);
end

if doProxy
    fig = figure('Color','w', 'Name',['DVproxy XY ' tag], ...
                 'Visible', local_fig_visible(cfg));
    ax = gca;
    cr3bp_plot_background(CJbg, SA.mu, ax);
    set(ax.Children,'HandleVisibility','off');
    hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
        % Faint non-overlap points (only when full row data is available)
    if ~isempty(Ax)
        scatter(ax, Ax, Ay, 10, [0.20 0.45 0.95], '.', ...
            'MarkerEdgeAlpha', 0.12, 'MarkerFaceAlpha', 0.12);
    end
    if ~isempty(Bx)
        scatter(ax, Bx, By, 10, [0.92 0.35 0.15], '.', ...
            'MarkerEdgeAlpha', 0.12, 'MarkerFaceAlpha', 0.12);
    end
    scatter(ax, x, y, 20, dv_proxy, 'filled', 'MarkerFaceAlpha', 0.94, 'MarkerEdgeAlpha', 0.30);
    cb = colorbar(ax); ylabel(cb, 'DVproxy = DVlb + DVpatch_{ub} (m/s)', 'Interpreter','tex');
    if isfinite(iMin) && iMin >= 1 && iMin <= numel(x)
        plot(ax, x(iMin), y(iMin), 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    end
    title(ax, sprintf('Voxel DVproxy heatmap | %s', tag), 'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'y');
    local_apply_zoom(cfg, ax);
    io_save_figure(fig, outdir, ['overlap_' safeTag '_bounds_proxy_xy'], cfg, 'Resolution', local_fig_res(cfg));
    local_close_if_hidden(cfg, fig);
end

end

% ===== helpers =====
function v = local_fig_visible(cfg)
v = 'off';
if isfield(cfg,'io') && isfield(cfg.io,'fig_visible') && ~isempty(cfg.io.fig_visible)
    v = cfg.io.fig_visible;
end
end

function res = local_fig_res(cfg)
res = 220;
if isfield(cfg,'io') && isfield(cfg.io,'fig_resolution') && ~isempty(cfg.io.fig_resolution)
    res = cfg.io.fig_resolution;
end
end

function local_close_if_hidden(cfg, fig)
v = local_fig_visible(cfg);
if ischar(v) || isstring(v)
    if strcmpi(char(v),'off')
        close(fig);
    end
end
end

function local_apply_zoom(cfg, ax)
if isfield(cfg,'diag') && isfield(cfg.diag,'zoom') && isfield(cfg.diag.zoom,'enable') && cfg.diag.zoom.enable
    if isfield(cfg.diag.zoom,'xlim') && ~isempty(cfg.diag.zoom.xlim), xlim(ax, cfg.diag.zoom.xlim); end
    if isfield(cfg.diag.zoom,'ylim') && ~isempty(cfg.diag.zoom.ylim), ylim(ax, cfg.diag.zoom.ylim); end
end
end

function v = local_cfg_get(cfg, path, defaultVal)
v = defaultVal;
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        k = parts{i};
        if ~isstruct(cur) || ~isfield(cur, k), return; end
        cur = cur.(k);
    end
    if ~isempty(cur), v = cur; end
catch
    v = defaultVal;
end
end

function tf = local_plot_enabled(cfg, path, defaultVal)
tf = logical(defaultVal);
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        k = parts{i};
        if ~isstruct(cur) || ~isfield(cur, k), return; end
        cur = cur.(k);
    end
    tf = logical(cur);
catch
    tf = logical(defaultVal);
end
end

function r = local_rows_cat(a, b)
if isstruct(a) && isstruct(b)
    nA = double(a.n); nB = double(b.n);
    if nA==0, r=b; return; end
    if nB==0, r=a; return; end
    r = atlas_rows_empty(nA+nB);
    r.n = uint32(nA+nB);
    r.iSeed(1:nA)    = a.iSeed(1:nA);    r.iSeed(nA+1:end)    = b.iSeed(1:nB);
    r.iHead(1:nA)    = a.iHead(1:nA);    r.iHead(nA+1:end)    = b.iHead(1:nB);
    r.leg(1:nA)      = a.leg(1:nA);      r.leg(nA+1:end)      = b.leg(1:nB);
    r.halfFlag(1:nA) = a.halfFlag(1:nA); r.halfFlag(nA+1:end) = b.halfFlag(1:nB);
    r.t(1:nA)        = a.t(1:nA);        r.t(nA+1:end)        = b.t(1:nB);
    r.ix(1:nA)       = a.ix(1:nA);       r.ix(nA+1:end)       = b.ix(1:nB);
    r.iy(1:nA)       = a.iy(1:nA);       r.iy(nA+1:end)       = b.iy(1:nB);
    r.it(1:nA)       = a.it(1:nA);       r.it(nA+1:end)       = b.it(1:nB);
    return;
end
if isempty(a), r=b; return; end
if isempty(b), r=a; return; end
r = [a; b];
end

function ids = local_rows_to_vid(rows, Ny, Nx, Nt)
if isempty(rows), ids = zeros(0,1); return; end
if isstruct(rows)
    n = double(rows.n);
    if n==0, ids = zeros(0,1); return; end
    iy = double(rows.iy(1:n));
    ix = double(rows.ix(1:n));
    it = double(rows.it(1:n));
else
    ix = double(rows(:,5)); iy = double(rows(:,6)); it = double(rows(:,7));
end
ix = max(1, min(Nx, ix));
iy = max(1, min(Ny, iy));
it = max(1, min(Nt, it));
ids = sub2ind([Ny, Nx, Nt], iy, ix, it);
end

function [x,y,th] = local_ids_to_centers(ids, grid3, Ny, Nx, Nt)
if isempty(ids)
    x = []; y = []; th = [];
    return;
end
[iy, ix, it] = ind2sub([Ny, Nx, Nt], ids);
x  = grid3.x_centers(ix);
y  = grid3.y_centers(iy);
th = grid3.th_centers(it);
end