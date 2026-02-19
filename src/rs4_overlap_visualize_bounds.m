function B = rs4_overlap_visualize_bounds(V, SA, SB, cfg, outdir, tag)
%RS4_OVERLAP_VISUALIZE_BOUNDS  Plot voxel-wise DVtotal bounds (vectorized).
%
% Bounds (m/s) per overlap voxel:
%   LB    = min(dv_turn_A) + min(dv_turn_B)
%   proxy = LB + DVpatch_ub   (DVpatch_ub uses voxel-center speed with CJ*)
%   UB    = max(dv_turn_A) + DVpatch_ub + max(dv_turn_B)
%
% V must be the flat-array struct from rs4_overlap_extract_voxel_info.
% All computations are fully vectorized — no per-voxel loop.

if nargin < 6 || isempty(tag),    tag    = 'bounds'; end
if nargin < 5 || isempty(outdir), outdir = pwd;      end
if ~exist(outdir,'dir'), mkdir(outdir); end

if ~isstruct(V) || ~isfield(V,'x') || isempty(V.x)
    warning('[rs4] Bounds visualization skipped: empty or old-format voxel metadata.');
    B = struct('dv_lb',[],'dv_ub',[],'dv_patch_ub',[],'dv_proxy',[], ...
               'min_dvproxy',NaN,'imin',NaN,'x_at_min',NaN,'y_at_min',NaN);
    return;
end

doLB    = local_plot_enabled(cfg, 'plot.rs4.bounds_lb',    true);
doUB    = local_plot_enabled(cfg, 'plot.rs4.bounds_ub',    false);
doProxy = local_plot_enabled(cfg, 'plot.rs4.bounds_proxy', true);

grid3   = SA.grid3;
safeTag = rs3_sanitize_fname(tag);
VU_mps  = local_cfg_get(cfg, 'units.VU_mps', 1.0);
CJstar  = min(SA.CJ, SB.CJ);

% ---------- flat arrays from V ----------
x  = V.x(:);
y  = V.y(:);

% ---------- vectorized voxel-center potential (one call) ----------
pot    = rs3_core_cr3bp_U_and_derivs(x, y, SA.mu);
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
    rs3_core_plot_cislunar_background(CJbg, SA.mu, ax);
    set(ax.Children,'HandleVisibility','off');
    hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
    scatter(ax, x, y, 18, dv_lb, 'filled', 'MarkerFaceAlpha', 0.92, 'MarkerEdgeAlpha', 0.25);
    cb = colorbar(ax); ylabel(cb, 'DVtotal lower bound (m/s)');
    title(ax, sprintf('Voxel DVtotal lower bound | %s', tag), 'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'y');
    local_apply_zoom(cfg, ax);
    rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_bounds_lb_xy'], cfg, 'Resolution', local_fig_res(cfg));
    local_close_if_hidden(cfg, fig);
end

if doUB && ~isempty(dv_ub)
    fig = figure('Color','w', 'Name',['DVtotal Upper Bound XY ' tag], ...
                 'Visible', local_fig_visible(cfg));
    ax = gca;
    rs3_core_plot_cislunar_background(CJbg, SA.mu, ax);
    set(ax.Children,'HandleVisibility','off');
    hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
    scatter(ax, x, y, 18, dv_ub, 'filled', 'MarkerFaceAlpha', 0.92, 'MarkerEdgeAlpha', 0.25);
    cb = colorbar(ax); ylabel(cb, 'DVtotal upper bound (m/s)');
    title(ax, sprintf('Voxel DVtotal upper bound | %s', tag), 'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'y');
    local_apply_zoom(cfg, ax);
    rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_bounds_ub_xy'], cfg, 'Resolution', local_fig_res(cfg));
    local_close_if_hidden(cfg, fig);
end

if doProxy
    fig = figure('Color','w', 'Name',['DVproxy XY ' tag], ...
                 'Visible', local_fig_visible(cfg));
    ax = gca;
    rs3_core_plot_cislunar_background(CJbg, SA.mu, ax);
    set(ax.Children,'HandleVisibility','off');
    hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
    scatter(ax, x, y, 20, dv_proxy, 'filled', 'MarkerFaceAlpha', 0.94, 'MarkerEdgeAlpha', 0.30);
    cb = colorbar(ax); ylabel(cb, 'DVproxy = DVlb + DVpatch_{ub} (m/s)', 'Interpreter','tex');
    if isfinite(iMin) && iMin >= 1 && iMin <= numel(x)
        plot(ax, x(iMin), y(iMin), 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    end
    title(ax, sprintf('Voxel DVproxy heatmap | %s', tag), 'Interpreter','none');
    xlabel(ax,'x'); ylabel(ax,'y');
    local_apply_zoom(cfg, ax);
    rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_bounds_proxy_xy'], cfg, 'Resolution', local_fig_res(cfg));
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
