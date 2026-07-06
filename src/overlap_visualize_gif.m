function gifPath = overlap_visualize_gif(O, SA, SB, cfg, outdir, tag, varargin)
%OVERLAP_VISUALIZE_GIF  Animated dark-theme GIF of the overlap-pair story.
%
% Storyboard:
%   1) background (ZVC) + the two periodic orbits (POs)
%   2) family A's forward reachable set (FRS) grows in
%   3) family B's backward reachable set (BRS) grows in
%   4) the overlap region is highlighted on top
%
% Usage:
%   gifPath = overlap_visualize_gif(O, SA, SB, cfg, outdir, tag)
%   gifPath = overlap_visualize_gif(..., 'FPS', 20, 'HoldStart', 1.5, ...)
%
% Name-value options (all optional):
%   'FPS'            frames/sec for the build stages             (default 20)
%   'HoldStart'      seconds pausing on background+POs            (default 1.5)
%   'HoldFRS'        seconds pausing once FRS is fully shown      (default 1.0)
%   'HoldBRS'        seconds pausing once BRS is fully shown      (default 1.0)
%   'HoldEnd'        seconds pausing on the final highlighted frame (default 3.0)
%   'NBuildFRS'      number of reveal frames for the FRS stage    (default 12)
%   'NBuildBRS'      number of reveal frames for the BRS stage    (default 12)
%   'NBuildOverlap'  number of reveal frames for the overlap stage (default 10)
%   'FigSize'        [width height] in pixels                     (default [900 700])
%   'Seed'           RNG seed controlling the point reveal order  (default 1)
%   'MaxPointsPerSet' cap on FRS/BRS/overlap points drawn per frame
%                     (default 20000 — dense, matching the standard overlap
%                     figure look). Redrawing every voxel on every frame is
%                     the main cost driver for slow GIF generation, so lower
%                     this (e.g. a few thousand) if a particular family pair
%                     is too slow to render. Set to Inf to animate every voxel.
%
% Output:
%   gifPath          full path to the written .gif

if nargin < 6 || isempty(tag), tag = 'overlap'; end
if nargin < 5 || isempty(outdir), outdir = pwd; end
if ~exist(outdir,'dir'), mkdir(outdir); end
safeTag = sanitize_fname(tag);

p = inputParser;
addParameter(p, 'FPS', 20);
addParameter(p, 'HoldStart', 1.5);
addParameter(p, 'HoldFRS', 1.0);
addParameter(p, 'HoldBRS', 1.0);
addParameter(p, 'HoldEnd', 3.0);
addParameter(p, 'NBuildFRS', 12);
addParameter(p, 'NBuildBRS', 12);
addParameter(p, 'NBuildOverlap', 10);
addParameter(p, 'FigSize', [900 700]);
addParameter(p, 'Seed', 1);
addParameter(p, 'MaxPointsPerSet', 20000);
parse(p, varargin{:});
opt = p.Results;

BG       = [0.0   0.0   0.0  ];   % true black background
FRS_COL  = [1.00 0.42 0.20];
BRS_COL  = [0.25 0.70 1.00];
OVL_COL  = [0.70 0.98 0.25];
PO_A_COL = 0.55*FRS_COL + 0.45*[1 1 1];
PO_B_COL = 0.55*BRS_COL + 0.45*[1 1 1];
TXTC     = [0.92 0.94 0.99];

[Ax, Ay, ~, Bx, By, ~, Ox, Oy, ~] = overlap_extract_xy_sets(SA, SB, O);

nA0 = numel(Ax); nB0 = numel(Bx); nO0 = numel(Ox);
[Ax, Ay] = overlap_subsample_xy(Ax, Ay, opt.MaxPointsPerSet, 11);
[Bx, By] = overlap_subsample_xy(Bx, By, opt.MaxPointsPerSet, 12);
[Ox, Oy] = overlap_subsample_xy(Ox, Oy, opt.MaxPointsPerSet, 13);
fprintf('[overlap_visualize_gif] animating %d/%d FRS, %d/%d BRS, %d/%d overlap points\n', ...
    numel(Ax), nA0, numel(Bx), nB0, numel(Ox), nO0);

rng(opt.Seed);
ordA = randperm(max(numel(Ax),1)); ordA = ordA(1:numel(Ax));
ordB = randperm(max(numel(Bx),1)); ordB = ordB(1:numel(Bx));
ordO = randperm(max(numel(Ox),1)); ordO = ordO(1:numel(Ox));

fig = figure('Color', BG, 'Name', ['Overlap process ' tag], 'Visible', 'off', ...
    'Units', 'pixels', 'Position', [0 0 opt.FigSize(1) opt.FigSize(2)]);
ax = axes('Parent', fig);

CJbg = min(SA.CJ, SB.CJ);
cr3bp_plot_background_dark(CJbg, SA.mu, ax);
set(ax.Children, 'HandleVisibility', 'off');
hold(ax, 'on'); axis(ax, 'equal'); grid(ax, 'on');

% Periodic orbits (closed loops), drawn once — static throughout
% (dense Xpo is only present when cfg.cache.store_dense_po=true / on a
% fresh build; otherwise fall back to the lightweight cached PO_xy trace)
local_plot_po(ax, local_po_xy(SA), PO_A_COL);
local_plot_po(ax, local_po_xy(SB), PO_B_COL);

title(ax, tag, 'Interpreter', 'none', 'Color', TXTC);
xlabel(ax, 'x', 'Color', TXTC); ylabel(ax, 'y', 'Color', TXTC);

local_apply_zoom(cfg, ax);

% Reveal handles (start empty)
% Markers are sized up relative to the static figure's since the GIF canvas
% renders at a much lower pixel resolution (no export DPI to fall back on) —
% at the static figure's marker size, dots round down to sub-pixel here and
% the cloud reads as sparse even at full point density.
hA     = scatter(ax, NaN, NaN, 16, FRS_COL, 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0);
hB     = scatter(ax, NaN, NaN, 16, BRS_COL, 'filled', 'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0);
hOglow = scatter(ax, NaN, NaN, 90, OVL_COL, 'filled', 'MarkerFaceAlpha', 0.20, 'MarkerEdgeAlpha', 0);
hO     = scatter(ax, NaN, NaN, 24, OVL_COL, 'filled', 'MarkerFaceAlpha', 0.95, ...
    'MarkerEdgeColor', [1 1 1], 'MarkerEdgeAlpha', 0.35);

lg = legend(ax, [hA hB hO], {'FRS (A)', 'BRS (B)', 'Overlap'}, 'Location', 'northwest');
set(lg, 'TextColor', TXTC, 'Color', [0.09 0.10 0.14], 'EdgeColor', [0.40 0.44 0.52], 'AutoUpdate', 'off');

cap_h = annotation(fig, 'textbox', [0.02 0.955 0.96 0.04], ...
    'String', 'Background + periodic orbits', 'Color', TXTC, 'FontSize', 10, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

imgs   = {};
delays = [];
dt_build = 1/opt.FPS;

% ---- Stage 1: background + POs ----
drawnow;
imgs{end+1} = local_capture(fig); delays(end+1) = opt.HoldStart; %#ok<AGROW>

% ---- Stage 2: FRS grows in ----
set(cap_h, 'String', sprintf('%s -- forward reachable set (FRS)', SA.name));
[imgs, delays] = local_reveal(fig, imgs, delays, hA, Ax, Ay, ordA, opt.NBuildFRS, dt_build, opt.HoldFRS);

% ---- Stage 3: BRS grows in ----
set(cap_h, 'String', sprintf('%s -- backward reachable set (BRS)', SB.name));
[imgs, delays] = local_reveal(fig, imgs, delays, hB, Bx, By, ordB, opt.NBuildBRS, dt_build, opt.HoldBRS);

% ---- Stage 4: overlap highlighted ----
set(cap_h, 'String', 'Overlap region');
[imgs, delays] = local_reveal_pair(fig, imgs, delays, hOglow, hO, Ox, Oy, ordO, ...
    opt.NBuildOverlap, dt_build, opt.HoldEnd);

close(fig);

% ---- Write GIF with one shared (global) colormap to avoid per-frame flicker ----
gifPath = fullfile(outdir, ['overlap_' safeTag '_process.gif']);
[~, gmap] = rgb2ind(imgs{end}, 256);
for k = 1:numel(imgs)
    imind = rgb2ind(imgs{k}, gmap, 'nodither');
    if k == 1
        imwrite(imind, gmap, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', delays(k));
    else
        imwrite(imind, gmap, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', delays(k));
    end
end

fprintf('[overlap_visualize_gif] Wrote %d frames -> %s\n', numel(imgs), gifPath);

end

% ===== helpers =====

function im = local_capture(fig)
fr = getframe(fig);
im = fr.cdata;
end

function xy = local_po_xy(S)
%LOCAL_PO_XY  Dense Xpo if available, else the cached lightweight PO_xy trace.
if isfield(S, 'Xpo') && ~isempty(S.Xpo)
    xy = S.Xpo(:,1:2);
elseif isfield(S, 'PO_xy') && ~isempty(S.PO_xy)
    xy = S.PO_xy;
else
    xy = zeros(0,2);
end
end

function local_plot_po(ax, xy, rgb)
if isempty(xy), return; end
xy = [xy(:,1:2); xy(1,1:2)];
plot(ax, xy(:,1), xy(:,2), '--', 'Color', rgb, 'LineWidth', 1.3, 'HandleVisibility', 'off');
end

function local_apply_zoom(cfg, ax)
if isfield(cfg,'diag') && isfield(cfg.diag,'zoom') && isfield(cfg.diag.zoom,'enable') && cfg.diag.zoom.enable
    if isfield(cfg.diag.zoom,'xlim') && ~isempty(cfg.diag.zoom.xlim), xlim(ax, cfg.diag.zoom.xlim); end
    if isfield(cfg.diag.zoom,'ylim') && ~isempty(cfg.diag.zoom.ylim), ylim(ax, cfg.diag.zoom.ylim); end
end
end

function [imgs, delays] = local_reveal(fig, imgs, delays, h, X, Y, ord, nBuild, dt, dtHold)
%LOCAL_REVEAL  Cumulatively reveal points on handle h over nBuild frames.
n = numel(X);
if n == 0
    drawnow;
    imgs{end+1} = local_capture(fig); delays(end+1) = dtHold; %#ok<AGROW>
    return;
end
for k = 1:nBuild
    m = max(1, round(n * k / nBuild));
    idx = ord(1:m);
    set(h, 'XData', X(idx), 'YData', Y(idx));
    drawnow;
    d = dt;
    if k == nBuild, d = dtHold; end
    imgs{end+1} = local_capture(fig); delays(end+1) = d; %#ok<AGROW>
end
end

function [imgs, delays] = local_reveal_pair(fig, imgs, delays, h1, h2, X, Y, ord, nBuild, dt, dtHold)
%LOCAL_REVEAL_PAIR  Same as local_reveal but drives two handles (glow + core) together.
n = numel(X);
if n == 0
    drawnow;
    imgs{end+1} = local_capture(fig); delays(end+1) = dtHold; %#ok<AGROW>
    return;
end
for k = 1:nBuild
    m = max(1, round(n * k / nBuild));
    idx = ord(1:m);
    set(h1, 'XData', X(idx), 'YData', Y(idx));
    set(h2, 'XData', X(idx), 'YData', Y(idx));
    drawnow;
    d = dt;
    if k == nBuild, d = dtHold; end
    imgs{end+1} = local_capture(fig); delays(end+1) = d; %#ok<AGROW>
end
end
