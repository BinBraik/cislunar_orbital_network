function rs4_overlap_visualize(O, SA, SB, cfg, outdir, tag)
%RS4_OVERLAP_VISUALIZE  Optional 2D/3D plots for overlap voxels.
% Shows FRS (red, faint), BRS (blue, faint), and overlap (green, solid).

if nargin < 6 || isempty(tag), tag = 'overlap'; end
if nargin < 5 || isempty(outdir), outdir = pwd; end
if ~exist(outdir,'dir'), mkdir(outdir); end

grid3 = SA.grid3;
safeTag = rs3_sanitize_fname(tag);

Nx  = numel(grid3.x_centers);
Ny_g = numel(grid3.y_centers);
Nt  = numel(grid3.th_centers);

% --- extract FRS, BRS, and overlap voxel centres ---
rowsA_F   = local_rows_cat(SA.Step4.rows_FRS_upper, SA.Step4.rows_FRS_lower);
rowsB_B_u = rs3_rows_mirror_lower(SB.Step4.rows_FRS_lower, grid3, 2);
rowsB_B_l = rs3_rows_mirror_lower(SB.Step4.rows_FRS_upper, grid3, 2);
rowsB_B   = local_rows_cat(rowsB_B_u, rowsB_B_l);

idsA = unique(local_rows_to_vid(rowsA_F, Ny_g, Nx, Nt));
idsB = unique(local_rows_to_vid(rowsB_B, Ny_g, Nx, Nt));
idsO = O.ids(:);

idsA_only = setdiff(idsA, idsO);
idsB_only = setdiff(idsB, idsO);

[Ax, Ay, Ath] = local_ids_to_centers(idsA_only, grid3, Ny_g, Nx, Nt);
[Bx, By, Bth] = local_ids_to_centers(idsB_only, grid3, Ny_g, Nx, Nt);
[Ox, Oy, Oth] = local_ids_to_centers(idsO,      grid3, Ny_g, Nx, Nt);

doXY  = local_plot_enabled(cfg, 'plot.rs4.overlap_xy', true);
doXYZ = local_plot_enabled(cfg, 'plot.rs4.overlap_xyz', true);

% ---------- 2D XY projection ----------
if doXY
    fig = figure('Color','w', 'Name',['Overlap XY ' tag], 'Visible', local_fig_visible(cfg));
    ax = gca;

    CJbg = min(SA.CJ, SB.CJ);
    rs3_core_plot_cislunar_background(CJbg, SA.mu, ax);
    set(ax.Children,'HandleVisibility','off');
    hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');

    % Faint FRS (red) and BRS (blue)
    hA = scatter(ax, Ax, Ay, 10, 'r', '.', ...
        'MarkerEdgeAlpha', 0.12, 'MarkerFaceAlpha', 0.12);
    hB = scatter(ax, Bx, By, 10, 'b', '.', ...
        'MarkerEdgeAlpha', 0.12, 'MarkerFaceAlpha', 0.12);

    % Overlap on top (green, solid)
    hO = scatter(ax, Ox, Oy, 16, [0.10 0.65 0.35], 'filled', ...
        'MarkerFaceAlpha', 0.95, 'MarkerEdgeAlpha', 0.20);

    title(ax, sprintf('Overlap XY projection | %s', tag), 'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'y');
    legend(ax, [hA hB hO], {'FRS (red)', 'BRS (blue)', 'Overlap'}, 'Location','northwest');

    local_apply_zoom(cfg, ax);
    rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_overlap_xy'], cfg, 'Resolution', local_fig_res(cfg));
    local_close_if_hidden(cfg, fig);
    
end

% ---------- 3D (x,y,theta) ----------
if doXYZ
    fig = figure('Color','w', 'Name',['Overlap 3D ' tag], 'Visible', local_fig_visible(cfg));
    ax = gca;

    hA = scatter3(ax, Ax, Ay, Ath, 8, 'r', '.', 'MarkerEdgeAlpha', 0.20);
    hold(ax,'on');
    hB = scatter3(ax, Bx, By, Bth, 8, 'b', '.', 'MarkerEdgeAlpha', 0.20);
    hO = scatter3(ax, Ox, Oy, Oth, 20, [0.10 0.65 0.35], 'o', 'filled', ...
        'MarkerFaceAlpha', 0.95, 'MarkerEdgeAlpha', 0.3);

    grid(ax,'on'); xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'\theta (rad)');
    title(ax, sprintf('Overlap in (x,y,\\theta) | %s', tag), 'Interpreter','none');
    legend(ax, [hA hB hO], {'FRS (red)', 'BRS (blue)', 'Overlap'}, 'Location','northwest');
    view(ax, 35, 20);

    rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_overlap_xyz'], cfg, 'Resolution', local_fig_res(cfg));
    local_close_if_hidden(cfg, fig);
end

end

% ===== helpers =====

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

function r = local_rows_cat(a, b)
if isstruct(a) && isstruct(b)
    nA = double(a.n); nB = double(b.n);
    if nA==0, r=b; return; end
    if nB==0, r=a; return; end
    r = rs3_rows_empty(nA+nB);
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

function tf = local_plot_enabled(cfg, path, defaultVal)
tf = logical(defaultVal);
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        k = parts{i};
        if ~isstruct(cur) || ~isfield(cur, k)
            return;
        end
        cur = cur.(k);
    end
    tf = logical(cur);
catch
    tf = logical(defaultVal);
end
end
