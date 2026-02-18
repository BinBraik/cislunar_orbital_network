function rs3_hits_visual_validate(S, cfg, outdir)
%RS3_HITS_VISUAL_VALIDATE  Visual checks for Step 4 hit logging.
%
% Updates:
%   - Plots FRS and BRS on separate figures.
%   - Adds LOWER half via CR3BP symmetry (Upper generated, Lower derived).
%   - Uses rs3_io_save_figure for consistent PNG/FIG saving.
%
% Expected inputs:
%   S.Step4.rows_FRS_upper  (packed rows struct or [N x 8] matrix)
%   S.Step4.rows_BRS_upper  (packed rows struct or [N x 8] matrix)
%   S.grid3 has x_centers, y_centers, th_edges, th_centers, Nx, Ny

if ~isfield(cfg,'io') || ~isfield(cfg.io,'save_figs') || ~cfg.io.save_figs
    return;
end
if nargin < 3 || isempty(outdir)
    outdir = pwd;
end
if ~exist(outdir,'dir'), mkdir(outdir); end

grid3 = S.grid3;
Nx = grid3.Nx; Ny = grid3.Ny;

rowsF_u = S.Step4.rows_FRS_upper;
rowsB_u = S.Step4.rows_BRS_upper;

% ---- Symmetry-complete lower halves (for plotting) ----
[rowsF_l, rowsB_l] = local_symmetry_complete(rowsF_u, rowsB_u, grid3);

% Full (upper + lower) for histograms, etc.
rowsF_full = local_rows_cat(rowsF_u, rowsF_l);
rowsB_full = local_rows_cat(rowsB_u, rowsB_l);

% ---- XY visited masks (projected), split by half ----
maskF_u = local_rows_to_mask_xy(rowsF_u, Ny, Nx);
maskF_l = local_rows_to_mask_xy(rowsF_l, Ny, Nx);
maskB_u = local_rows_to_mask_xy(rowsB_u, Ny, Nx);
maskB_l = local_rows_to_mask_xy(rowsB_l, Ny, Nx);

[xc, yc] = meshgrid(grid3.x_centers, grid3.y_centers);
safeName = rs3_sanitize_fname(S.name);

% ===================== FIG 1: FRS (upper+lower) =====================
fig = figure('Color','w','Name',['rs3 FRS hits ' S.name], 'Visible', local_get_fig_visible(cfg));
ax = gca;
rs3_core_plot_cislunar_background(S.CJ, S.mu, ax);
set(ax.Children,'HandleVisibility','off');
hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');

% Optional: show the periodic orbit
if local_show_orbit(cfg, S)
    stride = local_po_stride(cfg);
    plot(ax, S.Xpo(1:stride:end,1), S.Xpo(1:stride:end,2), 'k-', ...
        'LineWidth', 1.0, 'HandleVisibility','off');
end

h1 = plot(ax, xc(maskF_u), yc(maskF_u), '.', 'MarkerSize', 6);
h2 = plot(ax, xc(maskF_l), yc(maskF_l), '.', 'MarkerSize', 6);

title(ax, sprintf('FRS (upper + lower via symmetry) | %s', S.name), 'Interpreter','none');
legend(ax, [h1 h2], {'FRS upper (generated)', 'FRS lower (sym)'}, 'Location','best');

local_apply_zoom(cfg, ax);
rs3_io_save_figure(fig, outdir, ['step4_' safeName '_FRS_xy'], cfg, ...
    'Resolution', local_get_fig_resolution(cfg));
local_close_if_hidden(cfg, fig);

% ===================== FIG 2: BRS (upper+lower) =====================
fig = figure('Color','w','Name',['rs3 BRS hits ' S.name], 'Visible', local_get_fig_visible(cfg));
ax = gca;
rs3_core_plot_cislunar_background(S.CJ, S.mu, ax);
set(ax.Children,'HandleVisibility','off');
hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');

% Optional: show the periodic orbit
if local_show_orbit(cfg, S)
    stride = local_po_stride(cfg);
    plot(ax, S.Xpo(1:stride:end,1), S.Xpo(1:stride:end,2), 'k-', ...
        'LineWidth', 1.0, 'HandleVisibility','off');
end

h1 = plot(ax, xc(maskB_u), yc(maskB_u), '.', 'MarkerSize', 6);
h2 = plot(ax, xc(maskB_l), yc(maskB_l), '.', 'MarkerSize', 6);

title(ax, sprintf('BRS (upper + lower via symmetry) | %s', S.name), 'Interpreter','none');
legend(ax, [h1 h2], {'BRS upper (generated)', 'BRS lower (sym)'}, 'Location','best');

local_apply_zoom(cfg, ax);
rs3_io_save_figure(fig, outdir, ['step4_' safeName '_BRS_xy'], cfg, ...
    'Resolution', local_get_fig_resolution(cfg));
local_close_if_hidden(cfg, fig);

% ===================== FIG 3: Theta histograms (full) =====================
fig = figure('Color','w','Name',['rs3 theta hist ' S.name], 'Visible', local_get_fig_visible(cfg));
ax = gca;

thF = local_rows_to_theta(rowsF_full, grid3);
thB = local_rows_to_theta(rowsB_full, grid3);

if ~isempty(thF)
    histogram(ax, thF, 40);
    hold(ax,'on');
end
if ~isempty(thB)
    histogram(ax, thB, 40);
end

grid(ax,'on');
xlabel(ax,'theta (rad)'); ylabel(ax,'count');
title(ax, sprintf('Visited theta bins (upper+lower) | %s', S.name), 'Interpreter','none');
legend(ax, {'FRS full', 'BRS full'}, 'Location','best');

rs3_io_save_figure(fig, outdir, ['step4_' safeName '_theta_hist_full'], cfg, ...
    'Resolution', local_get_fig_resolution(cfg));
local_close_if_hidden(cfg, fig);

end

% =======================================================================
% Helpers
% =======================================================================

function [rowsF_l, rowsB_l] = local_symmetry_complete(rowsF_u, rowsB_u, grid3)
% CR3BP symmetry completion:
%   FRS_lower = mirror(BRS_upper) with leg=1
%   BRS_lower = mirror(FRS_upper) with leg=2

if isstruct(rowsF_u) && isstruct(rowsB_u)
    rowsF_l = rs3_rows_mirror_lower(rowsB_u, grid3, 1);
    rowsB_l = rs3_rows_mirror_lower(rowsF_u, grid3, 2);
    return;
end

% Matrix fallback ([N x 8]): [iSeed iHead leg t ix iy it halfFlag]
rowsF_l = local_mirror_lower_matrix(rowsB_u, grid3, 1);
rowsB_l = local_mirror_lower_matrix(rowsF_u, grid3, 2);
end

function rL = local_mirror_lower_matrix(rU, grid3, target_leg)
if isempty(rU)
    rL = zeros(0, 8);
    return;
end

Ny = grid3.Ny;

rL = rU;
rL(:,3) = target_leg;     % leg
rL(:,4) = -rU(:,4);       % time flip
rL(:,6) = (Ny - rU(:,6) + 1);  % iy mirror
rL(:,8) = -1;             % halfFlag = -1 (lower)

% Theta mirror: wrapToPi(pi - th)
th  = grid3.th_centers(rU(:,7));
thm = rs3_wrapToPi(pi - th);
itm = discretize(thm, grid3.th_edges);

bad = isnan(itm);
if any(bad)
    rL = rL(~bad,:);
    itm = itm(~bad);
end
rL(:,7) = itm;
end

function r = local_rows_cat(a, b)
if isstruct(a) && isstruct(b)
    r = rs3_rows_vcat(a, b);
else
    if isempty(a), r = b; return; end
    if isempty(b), r = a; return; end
    r = [a; b];
end
end

function mask = local_rows_to_mask_xy(rows, Ny, Nx)
mask = false(Ny, Nx);
if isempty(rows)
    return;
end
if isstruct(rows)
    n = double(rows.n);
    if n == 0, return; end
    ind = sub2ind([Ny, Nx], double(rows.iy(1:n)), double(rows.ix(1:n)));
    mask(ind) = true;
else
    % Matrix: col5=ix, col6=iy
    ind = sub2ind([Ny, Nx], rows(:,6), rows(:,5));
    mask(ind) = true;
end
end

function th = local_rows_to_theta(rows, grid3)
th = [];
if isempty(rows)
    return;
end
if isstruct(rows)
    n = double(rows.n);
    if n == 0, return; end
    th = grid3.th_centers(double(rows.it(1:n)));
else
    th = grid3.th_centers(rows(:,7));
end
end

function v = local_get_fig_visible(cfg)
v = 'off';
if isfield(cfg,'io') && isfield(cfg.io,'fig_visible') && ~isempty(cfg.io.fig_visible)
    v = cfg.io.fig_visible;
end
end

function res = local_get_fig_resolution(cfg)
res = 220;
if isfield(cfg,'io') && isfield(cfg.io,'fig_resolution') && ~isempty(cfg.io.fig_resolution)
    res = cfg.io.fig_resolution;
end
end

function local_close_if_hidden(cfg, fig)
v = local_get_fig_visible(cfg);
if ischar(v) || isstring(v)
    if strcmpi(char(v),'off')
        close(fig);
    end
end
end

function local_apply_zoom(cfg, ax)
if isfield(cfg,'diag') && isfield(cfg.diag,'zoom') && isfield(cfg.diag.zoom,'enable') && cfg.diag.zoom.enable
    if isfield(cfg.diag.zoom,'xlim') && ~isempty(cfg.diag.zoom.xlim)
        xlim(ax, cfg.diag.zoom.xlim);
    end
    if isfield(cfg.diag.zoom,'ylim') && ~isempty(cfg.diag.zoom.ylim)
        ylim(ax, cfg.diag.zoom.ylim);
    end
end
end

function tf = local_show_orbit(cfg, S)
tf = false;
if isfield(S,'Xpo') && ~isempty(S.Xpo)
    if isfield(cfg,'diag') && isfield(cfg.diag,'show_orbits')
        tf = logical(cfg.diag.show_orbits);
    else
        tf = true;
    end
end
end

function stride = local_po_stride(cfg)
stride = 8;
if isfield(cfg,'diag') && isfield(cfg.diag,'po_stride') && ~isempty(cfg.diag.po_stride)
    stride = cfg.diag.po_stride;
end
stride = max(1, round(stride));
end
