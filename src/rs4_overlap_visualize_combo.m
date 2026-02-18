function rs4_overlap_visualize_combo(SA, SB, O, cfg, outdir, tag)
%RS4_OVERLAP_VISUALIZE_COMBO
% Plot A.FRS_full, B.BRS_full, and their overlap together.
% Outputs:
%   - 2D XY combined plot on cislunar background
%   - 3D (x,y,theta) combined scatter

if nargin < 6 || isempty(tag), tag = 'combo'; end
if nargin < 5 || isempty(outdir), outdir = pwd; end
if ~exist(outdir,'dir'), mkdir(outdir); end

grid3 = SA.grid3;

% Robust dims
Nx = numel(grid3.x_centers);
Ny = numel(grid3.y_centers);
Nt = numel(grid3.th_centers);

safeTag = rs3_sanitize_fname(tag);

% ---------------- FULL sets ----------------
% A.FRS_full
rowsA_F_u = SA.Step4.rows_FRS_upper;
rowsA_B_u = SA.Step4.rows_BRS_upper;
rowsA_F_l = rs3_rows_mirror_lower(rowsA_B_u, grid3, 1);
rowsA_F   = local_rows_cat(rowsA_F_u, rowsA_F_l);

% B.BRS_full
rowsB_B_u = SB.Step4.rows_BRS_upper;
rowsB_F_u = SB.Step4.rows_FRS_upper;
rowsB_B_l = rs3_rows_mirror_lower(rowsB_F_u, grid3, 2);
rowsB_B   = local_rows_cat(rowsB_B_u, rowsB_B_l);

idsA = unique(local_rows_to_vid(rowsA_F, Ny, Nx, Nt));
idsB = unique(local_rows_to_vid(rowsB_B, Ny, Nx, Nt));
idsO = O.ids(:);  % already unique from overlap

% exclusive parts (nice visual)
idsA_only = setdiff(idsA, idsO);
idsB_only = setdiff(idsB, idsO);

% Convert to (x,y,theta) centers for plotting
[Ax, Ay, Ath] = local_ids_to_centers(idsA_only, grid3, Ny, Nx, Nt);
[Bx, By, Bth] = local_ids_to_centers(idsB_only, grid3, Ny, Nx, Nt);
[Ox, Oy, Oth] = local_ids_to_centers(idsO,      grid3, Ny, Nx, Nt);

% ---------------- 2D XY combined ----------------
fig = figure('Color','w', 'Name',['FRS+BRS+Overlap XY ' tag], 'Visible', local_fig_visible(cfg));
ax = gca;

CJbg = min(SA.CJ, SB.CJ);
rs3_core_plot_cislunar_background(CJbg, SA.mu, ax);
set(ax.Children,'HandleVisibility','off');
hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');

% Plot exclusive parts faint-ish by marker style
hA = plot(ax, Ax, Ay, '.', 'MarkerSize', 6);       % A only
hB = plot(ax, Bx, By, '.', 'MarkerSize', 6);       % B only
hO = plot(ax, Ox, Oy, 'o', 'MarkerSize', 4);       % overlap highlighted

title(ax, sprintf('A.FRS (full) + B.BRS (full) + Overlap | %s', tag), 'Interpreter','none');
xlabel(ax,'x'); ylabel(ax,'y');
legend(ax, [hA hB hO], {'A.FRS only', 'B.BRS only', 'Overlap'}, 'Location','best');

local_apply_zoom(cfg, ax);
rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_combo_xy'], cfg, 'Resolution', local_fig_res(cfg));
local_close_if_hidden(cfg, fig);

% ---------------- 2D XY voxel cells + DVpatch heat ----------------
fig = figure('Color','w', 'Name',['FRS+BRS cells + DVpatch heat ' tag], 'Visible', local_fig_visible(cfg));
ax = gca;

CJbg = min(SA.CJ, SB.CJ);
rs3_core_plot_cislunar_background(CJbg, SA.mu, ax);
set(ax.Children,'HandleVisibility','off');
hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');

% Build XY occupancy masks using actual voxel cells (not center scatters)
[iyA, ixA, ~] = ind2sub([Ny, Nx, Nt], idsA);
[iyB, ixB, ~] = ind2sub([Ny, Nx, Nt], idsB);

maskA = false(Ny, Nx);
maskB = false(Ny, Nx);
maskA(sub2ind([Ny, Nx], iyA, ixA)) = true;
maskB(sub2ind([Ny, Nx], iyB, ixB)) = true;

hA = local_draw_xy_cells(ax, grid3.x_edges, grid3.y_edges, maskA, [0.20 0.45 0.95], 0.14, 0.18);
hB = local_draw_xy_cells(ax, grid3.x_edges, grid3.y_edges, maskB, [0.92 0.35 0.15], 0.14, 0.18);

% DVpatch upper bound heat per overlap voxel using box speed at center
% CJ* = min(CJ_A, CJ_B), dv_patch_ub = 2*v_box*sin(dtheta/2), reported in m/s.
idsO = idsO(:);
dv_ub = zeros(numel(idsO),1);
[iyO, ixO, ~] = ind2sub([Ny, Nx, Nt], idsO);
xO = grid3.x_centers(ixO);
yO = grid3.y_centers(iyO);

CJstar = min(SA.CJ, SB.CJ);
VU_mps = local_cfg_get(cfg, 'units.VU_mps', 1.0);
for k = 1:numel(idsO)
    pot = rs3_core_cr3bp_U_and_derivs(xO(k), yO(k), SA.mu);
    v_box = sqrt(max(2*pot.U - CJstar, 0));
    dv_ub(k) = 2*v_box*sin(abs(grid3.dtheta)/2) * VU_mps;
end

dvXY = nan(Ny, Nx);
if ~isempty(dv_ub)
    linO = sub2ind([Ny, Nx], iyO, ixO);
    dvXYmax = accumarray(linO, dv_ub, [Ny*Nx, 1], @max, NaN);
    dvXY = reshape(dvXYmax, [Ny, Nx]);
end

maskOxy = ~isnan(dvXY);
hO = local_draw_xy_heat_cells(ax, grid3.x_edges, grid3.y_edges, maskOxy, dvXY, 0.85, 0.22);

if ~isempty(hO)
    cb = colorbar(ax);
    ylabel(cb, '\DeltaV_{patch,ub} (m/s)', 'Interpreter', 'tex');
end

title(ax, sprintf('XY voxel footprint + overlap \DeltaV_{patch,ub} heat | %s', tag), 'Interpreter','none');
xlabel(ax,'x'); ylabel(ax,'y');
legend(ax, [hA hB hO], {'A.FRS XY cells', 'B.BRS XY cells', 'Overlap DVpatch UB'}, 'Location','best');

local_apply_zoom(cfg, ax);
rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_combo_xy_dvpatch_heat'], cfg, 'Resolution', local_fig_res(cfg));
local_close_if_hidden(cfg, fig);

% ---------------- 3D (x,y,theta) combined ----------------
fig = figure('Color','w', 'Name',['FRS+BRS+Overlap 3D ' tag], 'Visible', local_fig_visible(cfg));
ax = gca;

% Use different marker types so it’s clear even without color assumptions
hA = scatter3(ax, Ax, Ay, Ath, 8, '.', 'MarkerEdgeAlpha', 0.7);
hold(ax,'on');
hB = scatter3(ax, Bx, By, Bth, 8, '.', 'MarkerEdgeAlpha', 0.7);
hO = scatter3(ax, Ox, Oy, Oth, 18, 'o', 'filled');  % overlap strong

grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'\theta (rad)');
title(ax, sprintf('Overlap in (x,y,\\theta) | %s', tag), 'Interpreter','none');
legend(ax, [hA hB hO], {'A.FRS only', 'B.BRS only', 'Overlap'}, 'Location','best');
view(ax, 35, 20);

rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_combo_xyz'], cfg, 'Resolution', local_fig_res(cfg));
local_close_if_hidden(cfg, fig);

end

% ================= helpers =================

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

function h = local_draw_xy_cells(ax, x_edges, y_edges, maskXY, faceColor, faceAlpha, edgeAlpha)
h = gobjects(0);
[Ny, Nx] = size(maskXY);
for iy = 1:Ny
    for ix = 1:Nx
        if ~maskXY(iy, ix), continue; end
        x1 = x_edges(ix); x2 = x_edges(ix+1);
        y1 = y_edges(iy); y2 = y_edges(iy+1);
        h(end+1) = patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], faceColor, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', faceColor, 'EdgeAlpha', edgeAlpha, ...
            'LineWidth', 0.5); %#ok<AGROW>
    end
end
if isempty(h)
    h = patch(ax, nan, nan, faceColor, 'FaceAlpha', faceAlpha, 'EdgeColor', faceColor, 'EdgeAlpha', edgeAlpha);
else
    h = h(1);
end
end

function h = local_draw_xy_heat_cells(ax, x_edges, y_edges, maskXY, valXY, faceAlpha, edgeAlpha)
h = gobjects(0);
[Ny, Nx] = size(maskXY);

vals = valXY(maskXY);
if isempty(vals)
    h = patch(ax, nan, nan, [0.4 0.4 0.4], 'FaceAlpha', faceAlpha, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', edgeAlpha);
    return;
end

vmin = min(vals);
vmax = max(vals);
if abs(vmax - vmin) < 1e-14
    vmax = vmin + 1e-14;
end

cmap = parula(256);
for iy = 1:Ny
    for ix = 1:Nx
        if ~maskXY(iy, ix), continue; end
        v = valXY(iy, ix);
        a = (v - vmin) / (vmax - vmin);
        a = max(0, min(1, a));
        ci = 1 + floor(a * (size(cmap,1)-1));
        c = cmap(ci,:);

        x1 = x_edges(ix); x2 = x_edges(ix+1);
        y1 = y_edges(iy); y2 = y_edges(iy+1);
        h(end+1) = patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], c, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', [0.12 0.12 0.12], ...
            'EdgeAlpha', edgeAlpha, 'LineWidth', 0.6); %#ok<AGROW>
    end
end

colormap(ax, cmap);
caxis(ax, [vmin, vmax]);

if isempty(h)
    h = patch(ax, nan, nan, [0.7 0.7 0.7], 'FaceAlpha', faceAlpha, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', edgeAlpha);
else
    h = h(1);
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
