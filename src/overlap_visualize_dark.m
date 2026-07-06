function overlap_visualize_dark(O, SA, SB, cfg, outdir, tag, varargin)
%OVERLAP_VISUALIZE_DARK  Dark-theme, presentation/paper-style overlap figure.
%
% Same FRS (A) / BRS (B) / overlap voxel projection as overlap_visualize.m,
% restyled for a dark background and exported as .fig (MATLAB), .pdf, and
% .eps suitable for a paper or poster.
%
% Usage:
%   overlap_visualize_dark(O, SA, SB, cfg, outdir, tag)
%   overlap_visualize_dark(..., 'MaxPointsPerSet', 8000, 'Resolution', 300)
%
% Name-value options (all optional):
%   'MaxPointsPerSet'  cap on FRS/BRS/overlap points actually drawn (default 8000).
%                      Voxel clouds can run into the tens/hundreds of thousands;
%                      at this marker size a random subsample looks identical
%                      but renders and exports far faster and smaller. Set to
%                      Inf to plot every voxel.
%   'Resolution'       DPI used to rasterize the PDF/EPS export (default 300).
%                      Vector export of a huge scatter embeds one path per
%                      point and can balloon to hundreds of MB; rasterizing
%                      at a print-quality DPI keeps files small and fast to
%                      generate while remaining crisp.

if nargin < 6 || isempty(tag), tag = 'overlap'; end
if nargin < 5 || isempty(outdir), outdir = pwd; end
if ~exist(outdir,'dir'), mkdir(outdir); end

p = inputParser;
addParameter(p, 'MaxPointsPerSet', 8000);
addParameter(p, 'Resolution', 300);
parse(p, varargin{:});
opt = p.Results;

safeTag = sanitize_fname(tag);

BG      = [0.047 0.055 0.078];
FRS_COL = [1.00 0.42 0.20];   % warm bright orange-red
BRS_COL = [0.25 0.70 1.00];   % bright cyan-blue
OVL_COL = [0.70 0.98 0.25];   % bright lime-green
TXTC    = [0.92 0.94 0.99];
LEGBG   = [0.09 0.10 0.14];
LEGEDGE = [0.40 0.44 0.52];

[Ax, Ay, ~, Bx, By, ~, Ox, Oy, ~] = overlap_extract_xy_sets(SA, SB, O);

nA0 = numel(Ax); nB0 = numel(Bx); nO0 = numel(Ox);
[Ax, Ay] = overlap_subsample_xy(Ax, Ay, opt.MaxPointsPerSet, 1);
[Bx, By] = overlap_subsample_xy(Bx, By, opt.MaxPointsPerSet, 2);
[Ox, Oy] = overlap_subsample_xy(Ox, Oy, opt.MaxPointsPerSet, 3);
fprintf('[overlap_visualize_dark] plotting %d/%d FRS, %d/%d BRS, %d/%d overlap points\n', ...
    numel(Ax), nA0, numel(Bx), nB0, numel(Ox), nO0);

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
local_save_dark(fig, outdir, baseName, BG, opt.Resolution);
local_close_if_hidden(cfg, fig);

end

% ===== helpers =====

function local_save_dark(fig, outdir, baseName, BG, res)
%LOCAL_SAVE_DARK  Save .fig, .pdf, .eps preserving the dark background.
%
% PDF/EPS are rasterized (ContentType 'image') at `res` DPI rather than
% exported as vector: a vector export of a large scatter cloud embeds one
% drawing primitive per point and can reach hundreds of MB, whereas a
% rasterized embed at print DPI is a few MB and generates in a fraction of
% the time.
set(fig, 'Color', BG, 'InvertHardcopy', 'off');

figPath = fullfile(outdir, [baseName '.fig']);
try
    savefig(fig, figPath);
catch ME
    warning('[overlap_visualize_dark] Failed to save FIG (%s): %s', figPath, ME.message);
end

pdfPath = fullfile(outdir, [baseName '.pdf']);
try
    exportgraphics(fig, pdfPath, 'ContentType', 'image', 'Resolution', res, 'BackgroundColor', 'current');
catch ME
    warning('[overlap_visualize_dark] Failed to save PDF (%s): %s', pdfPath, ME.message);
end

epsPath = fullfile(outdir, [baseName '.eps']);
try
    exportgraphics(fig, epsPath, 'ContentType', 'image', 'Resolution', res, 'BackgroundColor', 'current');
catch
    try
        print(fig, epsPath, sprintf('-r%d', res), '-depsc');
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
