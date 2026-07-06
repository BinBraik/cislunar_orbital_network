function overlap_visualize_dark(O, SA, SB, cfg, outdir, tag)
%OVERLAP_VISUALIZE_DARK  Dark-theme, presentation/paper-style overlap figure.
%
% Same FRS (A) / BRS (B) / overlap voxel projection as overlap_visualize.m,
% restyled for a dark background and exported as vector formats suitable
% for a paper or poster: .fig (MATLAB), .pdf, and .eps.
%
% Usage:
%   overlap_visualize_dark(O, SA, SB, cfg, outdir, tag)

if nargin < 6 || isempty(tag), tag = 'overlap'; end
if nargin < 5 || isempty(outdir), outdir = pwd; end
if ~exist(outdir,'dir'), mkdir(outdir); end

safeTag = sanitize_fname(tag);

BG      = [0.047 0.055 0.078];
FRS_COL = [1.00 0.42 0.20];   % warm bright orange-red
BRS_COL = [0.25 0.70 1.00];   % bright cyan-blue
OVL_COL = [0.70 0.98 0.25];   % bright lime-green
TXTC    = [0.92 0.94 0.99];
LEGBG   = [0.09 0.10 0.14];
LEGEDGE = [0.40 0.44 0.52];

[Ax, Ay, ~, Bx, By, ~, Ox, Oy, ~] = overlap_extract_xy_sets(SA, SB, O);

fig = figure('Color', BG, 'Name', ['Overlap XY (dark) ' tag], ...
    'Visible', local_fig_visible(cfg), 'InvertHardcopy', 'off', ...
    'Units', 'inches', 'Position', [0 0 7 5.4]);
ax = axes('Parent', fig);

CJbg = min(SA.CJ, SB.CJ);
cr3bp_plot_background_dark(CJbg, SA.mu, ax);
set(ax.Children, 'HandleVisibility', 'off');
hold(ax, 'on'); axis(ax, 'equal'); grid(ax, 'on');

hA = scatter(ax, Ax, Ay, 10, FRS_COL, '.', ...
    'MarkerEdgeAlpha', 0.30, 'MarkerFaceAlpha', 0.30);
hB = scatter(ax, Bx, By, 10, BRS_COL, '.', ...
    'MarkerEdgeAlpha', 0.30, 'MarkerFaceAlpha', 0.30);

% Overlap on top: soft glow halo + solid bright core
scatter(ax, Ox, Oy, 60, OVL_COL, 'filled', ...
    'MarkerFaceAlpha', 0.18, 'MarkerEdgeAlpha', 0, 'HandleVisibility', 'off');
hO = scatter(ax, Ox, Oy, 18, OVL_COL, 'filled', ...
    'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor', [1 1 1], 'MarkerEdgeAlpha', 0.35);

title(ax, sprintf('Overlap XY projection | %s', tag), 'Interpreter', 'none', 'Color', TXTC);
xlabel(ax, 'x', 'Color', TXTC); ylabel(ax, 'y', 'Color', TXTC);
lg = legend(ax, [hA hB hO], {'FRS (A)', 'BRS (B)', 'Overlap'}, 'Location', 'northwest');
set(lg, 'TextColor', TXTC, 'Color', LEGBG, 'EdgeColor', LEGEDGE);

local_apply_zoom(cfg, ax);

baseName = ['overlap_' safeTag '_overlap_xy_dark'];
local_save_dark(fig, outdir, baseName, BG);
local_close_if_hidden(cfg, fig);

end

% ===== helpers =====

function local_save_dark(fig, outdir, baseName, BG)
%LOCAL_SAVE_DARK  Save .fig, .pdf, .eps preserving the dark background.
set(fig, 'Color', BG, 'InvertHardcopy', 'off');

figPath = fullfile(outdir, [baseName '.fig']);
try
    savefig(fig, figPath);
catch ME
    warning('[overlap_visualize_dark] Failed to save FIG (%s): %s', figPath, ME.message);
end

pdfPath = fullfile(outdir, [baseName '.pdf']);
try
    exportgraphics(fig, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'current');
catch ME
    warning('[overlap_visualize_dark] Failed to save PDF (%s): %s', pdfPath, ME.message);
end

epsPath = fullfile(outdir, [baseName '.eps']);
try
    exportgraphics(fig, epsPath, 'ContentType', 'vector', 'BackgroundColor', 'current');
catch
    try
        print(fig, epsPath, '-depsc', '-painters');
    catch ME
        warning('[overlap_visualize_dark] Failed to save EPS (%s): %s', epsPath, ME.message);
    end
end
end

function v = local_fig_visible(cfg)
v = 'off';
if isfield(cfg,'io') && isfield(cfg.io,'fig_visible') && ~isempty(cfg.io.fig_visible)
    v = cfg.io.fig_visible;
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
