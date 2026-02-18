function rs4_overlap_visualize(O, SA, SB, cfg, outdir, tag)
%RS4_OVERLAP_VISUALIZE  2D + 3D plots for overlap voxels.

if nargin < 6 || isempty(tag), tag = 'overlap'; end
if nargin < 5 || isempty(outdir), outdir = pwd; end
if ~exist(outdir,'dir'), mkdir(outdir); end

grid3 = SA.grid3;
safeTag = rs3_sanitize_fname(tag);

% ---------- 2D XY projection ----------
fig = figure('Color','w', 'Name',['Overlap XY ' tag], 'Visible', local_fig_visible(cfg));
ax = gca;

% background from A (assumes same CJ/mu conceptually; you can decide what to show)
CJbg = min(SA.CJ, SB.CJ);
rs3_core_plot_cislunar_background(CJbg, SA.mu, ax);
set(ax.Children,'HandleVisibility','off');
hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');

% projection: collapse theta, show counts as “intensity”
Ny=grid3.Ny; Nx=grid3.Nx;
countXY = accumarray([double(O.iy), double(O.ix)], 1, [Ny, Nx], @sum, 0);
mask = countXY > 0;

[xc, yc] = meshgrid(grid3.x_centers, grid3.y_centers);
plot(ax, xc(mask), yc(mask), '.', 'MarkerSize', 7);

title(ax, sprintf('Overlap XY projection | %s', tag), 'Interpreter','none');
xlabel(ax,'x'); ylabel(ax,'y');

local_apply_zoom(cfg, ax);
rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_overlap_xy'], cfg, 'Resolution', local_fig_res(cfg));
local_close_if_hidden(cfg, fig);

% ---------- 3D (x,y,theta) ----------
fig = figure('Color','w', 'Name',['Overlap 3D ' tag], 'Visible', local_fig_visible(cfg));
ax = gca;
scatter3(ax, O.x, O.y, O.th, 8, 'filled'); % no explicit colors
grid(ax,'on'); xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'\theta (rad)');
title(ax, sprintf('Overlap in (x,y,\\theta) | %s', tag), 'Interpreter','none');
view(ax, 35, 20);

rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_overlap_xyz'], cfg, 'Resolution', local_fig_res(cfg));
local_close_if_hidden(cfg, fig);

end

% ===== helpers
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
