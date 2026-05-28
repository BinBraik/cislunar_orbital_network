function io_save_figure(fig, outdir, baseName, cfg, varargin)
%IO_SAVE_FIGURE  Save figure as PNG and/or FIG based on cfg.
%
% Usage:
%   io_save_figure(fig, outdir, baseName, cfg)
%   io_save_figure(fig, outdir, baseName, cfg, 'Resolution', 220)
%
% Notes:
%   - Uses cfg.io.save_figs as master plot-save enable.
%   - Saves PNG by default when save_figs=true.
%   - Saves FIG when cfg.io.save_fig=true.
%   - Stores under outdir/[cfg.io.fig_subdir]/ if provided.

if nargin < 4 || isempty(cfg)
    cfg = struct();
end

% master switch
saveFigs = true;
if isfield(cfg,'io') && isfield(cfg.io,'save_figs')
    saveFigs = logical(cfg.io.save_figs);
end
if ~saveFigs
    return;
end

if nargin < 2 || isempty(outdir), outdir = pwd; end
if nargin < 3 || isempty(baseName), baseName = 'figure'; end

subdir = '';
if isfield(cfg,'io') && isfield(cfg.io,'fig_subdir') && ~isempty(cfg.io.fig_subdir)
    subdir = cfg.io.fig_subdir;
end

if ~isempty(subdir)
    outdir2 = fullfile(outdir, subdir);
else
    outdir2 = outdir;
end
if ~exist(outdir2,'dir'), mkdir(outdir2); end

% parse optional args
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'Resolution',220);
parse(p,varargin{:});
res = p.Results.Resolution;

% Always save PNG when save_figs enabled
pngPath = fullfile(outdir2, [baseName '.png']);
try
    if exist('exportgraphics','file')
        exportgraphics(fig, pngPath, 'Resolution', res);
    else
        saveas(fig, pngPath);
    end
catch ME
    warning('[atlas] Failed to save PNG (%s): %s', pngPath, ME.message);
end

% Save FIG if requested
saveFig = false;
if isfield(cfg,'io') && isfield(cfg.io,'save_fig')
    saveFig = logical(cfg.io.save_fig);
end
if saveFig
    figPath = fullfile(outdir2, [baseName '.fig']);
    try
        savefig(fig, figPath);
    catch ME
        warning('[atlas] Failed to save FIG (%s): %s', figPath, ME.message);
    end
end

% Save SVG if requested (default true)
saveSvg = true;
if isfield(cfg,'io') && isfield(cfg.io,'save_svg')
    saveSvg = logical(cfg.io.save_svg);
end
if saveSvg
    svgPath = fullfile(outdir2, [baseName '.svg']);
    try
        if exist('exportgraphics','file')
            exportgraphics(fig, svgPath, 'ContentType', 'vector');
        else
            print(fig, svgPath, '-dsvg');
        end
    catch ME
        warning('[atlas] Failed to save SVG (%s): %s', svgPath, ME.message);
    end
end

end
