function rs4_voxel_traj_visualize_compare(T, R, SA, SB, cfg, outdir, tag)
%RS4_VOXEL_TRAJ_VISUALIZE_COMPARE Compare baseline and DC-corrected trajectories.
%
% Shows baseline A/B arcs and closest-point markers from T, and overlays
% differential-corrector (DC) arcs and patch states from R. The figure
% annotates DV comparisons, residual statistics, and convergence state.
%
% Inputs
%   T      : output struct from rs4_voxel_traj_extract (baseline trajectory)
%   R      : output struct from rs4_voxel_dc_solve (DC solution)
%   SA, SB : family structs
%   cfg    : config struct
%   outdir : output directory for saving figure
%   tag    : string tag for filename and title

if nargin < 6 || isempty(outdir), outdir = pwd; end
if nargin < 7 || isempty(tag),    tag = 'traj'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end

grid3 = SA.grid3;
CJbg  = min(SA.CJ, SB.CJ);
mu    = SA.mu;

relTol = local_cfg_get(cfg, 'propag.relTol', 1e-9);
absTol = local_cfg_get(cfg, 'propag.absTol', 1e-9);
TU_days = local_cfg_get(cfg, 'units.TU_days', 1.0);

% ---- build DC-corrected arcs from solved decision vector ----
[XA_dc, xB_dc, yB_dc] = local_build_dc_arcs(T, R, SA, SB, relTol, absTol);

% matched point (minimum A/B separation among integrated points)
[~, ia, ib] = local_min_pair_distance(XA_dc(:,1), XA_dc(:,2), xB_dc, yB_dc);

% ---- figure ----
fig = figure('Color','w', 'Name', ['Transfer Trajectory Compare ' tag], ...
             'Visible', local_fig_visible(cfg), ...
             'Units','pixels', 'Position',[100 100 900 720]);
ax = axes('Parent', fig);
rs3_core_plot_cislunar_background(CJbg, mu, ax);
set(ax.Children, 'HandleVisibility', 'off');
hold(ax, 'on');
axis(ax, 'equal');
grid(ax, 'on');

% ---- periodic orbits ----
local_plot_po(ax, SA, [0.15 0.45 0.80], 'Origin PO');
local_plot_po(ax, SB, [0.80 0.20 0.15], 'Target PO');

% ---- baseline arcs ----
plot(ax, T.XA(:,1), T.XA(:,2), '-', ...
    'Color', [0.10 0.40 0.85], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Baseline A (%.1f d)', T.tof_A_days));
plot(ax, T.x_B, T.y_B, '-', ...
    'Color', [0.85 0.15 0.10], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Baseline B (%.1f d)', T.tof_B_days));

% baseline patch markers (closest-to-voxel points)
plot(ax, T.XA(T.i_star,1), T.XA(T.i_star,2), 'o', ...
    'Color', [0.10 0.40 0.85], 'MarkerFaceColor', [0.10 0.40 0.85], ...
    'MarkerSize', 6, 'DisplayName', 'Baseline A patch marker');
plot(ax, T.x_B(T.j_star), T.y_B(T.j_star), 'o', ...
    'Color', [0.85 0.15 0.10], 'MarkerFaceColor', [0.85 0.15 0.10], ...
    'MarkerSize', 6, 'DisplayName', 'Baseline B patch marker');

% ---- DC-corrected arcs ----
plot(ax, XA_dc(:,1), XA_dc(:,2), '--', ...
    'Color', [0.10 0.40 0.85], 'LineWidth', 1.8, ...
    'DisplayName', sprintf('DC A (%.1f d)', max(R.t_A,0) * TU_days));
plot(ax, xB_dc, yB_dc, '--', ...
    'Color', [0.85 0.15 0.10], 'LineWidth', 1.8, ...
    'DisplayName', sprintf('DC B (%.1f d)', max(R.t_B,0) * TU_days));

% matched point markers from DC arc pair
plot(ax, XA_dc(ia,1), XA_dc(ia,2), 'd', ...
    'Color', [0.10 0.40 0.85], 'MarkerFaceColor', [1 1 1], ...
    'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'DC matched A point');
plot(ax, xB_dc(ib), yB_dc(ib), 'd', ...
    'Color', [0.85 0.15 0.10], 'MarkerFaceColor', [1 1 1], ...
    'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'DC matched B point');

% patch states from DC residual evaluation
plot(ax, R.stateA_patch.x, R.stateA_patch.y, 'p', ...
    'Color', [0.10 0.40 0.85], 'MarkerFaceColor', [0.10 0.40 0.85], ...
    'MarkerSize', 9, 'DisplayName', 'DC patch state A');
plot(ax, R.stateB_patch.x, R.stateB_patch.y, 'p', ...
    'Color', [0.85 0.15 0.10], 'MarkerFaceColor', [0.85 0.15 0.10], ...
    'MarkerSize', 9, 'DisplayName', 'DC patch state B');

% voxel rectangle + center
rectangle('Parent', ax, ...
    'Position', [T.xc - grid3.dx/2, T.yc - grid3.dy/2, grid3.dx, grid3.dy], ...
    'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.2, 'FaceColor', 'none');
plot(ax, T.xc, T.yc, 'k+', 'MarkerSize', 6, 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

% ---- annotation ----
res = R.residual(:);
resMean = mean(res);
resStd = std(res);
resMax = max(abs(res));
convFlag = 'false';
if R.converged
    convFlag = 'true';
end
convState = sprintf('%s (%s)', convFlag, R.message);

dvText = sprintf(['DV_{patch,true}=%.2f m/s   DV_{patch,dc}=%.2f m/s\\n' ...
                  'DV_{total,true}=%.2f m/s   DV_{total,dc}=%.2f m/s'], ...
                 T.DV_patch_mps, R.DV_patch_dc_mps, ...
                 T.DV_total_true_mps, R.DV_total_dc_mps);

resText = sprintf(['Residual norm=%.3e   mean=%.3e   std=%.3e   max|r|=%.3e\\n' ...
                   'iters=%d   converged=%s'], ...
                  R.final_residual_norm, resMean, resStd, resMax, ...
                  R.iterations, convState);

text(ax, 0.02, 0.98, dvText, 'Units','normalized', ...
    'VerticalAlignment','top', 'HorizontalAlignment','left', ...
    'BackgroundColor', [1 1 1], 'Margin', 6, ...
    'FontSize', 8, 'Interpreter', 'tex');

text(ax, 0.02, 0.86, resText, 'Units','normalized', ...
    'VerticalAlignment','top', 'HorizontalAlignment','left', ...
    'BackgroundColor', [1 1 1], 'Margin', 6, ...
    'FontSize', 8, 'Interpreter', 'tex');

legend(ax, 'Location', 'best', 'FontSize', 8);

title(ax, {sprintf('%s  \rightarrow  %s  |  voxel #%d', SA.name, SB.name, T.vid), ...
           'Baseline vs. Differential-Corrected Trajectories'}, ...
      'Interpreter', 'tex', 'FontSize', 10);
xlabel(ax, 'x (nd)');
ylabel(ax, 'y (nd)');

% ---- Save ----
safeTag = regexprep(tag, '[^A-Za-z0-9_]', '_');
baseName = ['rs4_' safeTag '_traj_compare_dc'];
rs3_io_save_figure(fig, outdir, baseName, cfg);
fprintf('[rs4] Saved compare figure: %s\n', local_saved_png_path(outdir, baseName, cfg));
end

% =========================================================================
% Local helpers
% =========================================================================

function [XA_dc, xB_dc, yB_dc] = local_build_dc_arcs(T, R, SA, SB, relTol, absTol)
seedA = T.seed_A(:);
seedB = T.seed_B_frs(:);

thA = rs3_wrapToPi(seedA(3) + R.phi_A + R.delta_A);
thB = rs3_wrapToPi(seedB(3) + R.phi_B + R.delta_B);

IC_A = [seedA(1); seedA(2); thA];
IC_B = [seedB(1); seedB(2); thB];

odeOpts = odeset('RelTol', relTol, 'AbsTol', absTol);

[~, XA_dc] = ode113(@(tt,XX) rs3_core_reduced_cr3bp_model(tt, XX, SA.CJ, SA.mu, false), ...
                    [0, max(R.t_A,0)], IC_A, odeOpts);
[~, XB_frs_dc] = ode113(@(tt,XX) rs3_core_reduced_cr3bp_model(tt, XX, SB.CJ, SB.mu, false), ...
                        [0, max(R.t_B,0)], IC_B, odeOpts);

xB_dc = XB_frs_dc(:,1);
yB_dc = -XB_frs_dc(:,2);
end

function [dmin, ia, ib] = local_min_pair_distance(xA, yA, xB, yB)
NA = numel(xA);
NB = numel(xB);

best = inf;
ia = 1;
ib = 1;
for i = 1:NA
    dx = xB - xA(i);
    dy = yB - yA(i);
    [d2, j] = min(dx.^2 + dy.^2);
    if d2 < best
        best = d2;
        ia = i;
        ib = j;
    end
end

dmin = sqrt(best);
end

function local_plot_po(ax, S, rgb, dispName)
hasPOxy = isfield(S, 'PO_xy') && ~isempty(S.PO_xy);
hasXpo  = isfield(S, 'Xpo')   && ~isempty(S.Xpo);

if hasPOxy
    xy = S.PO_xy;
elseif hasXpo
    xy = S.Xpo(:, 1:2);
else
    if isfield(S, 'SeedsUpper') && ~isempty(S.SeedsUpper)
        plot(ax, S.SeedsUpper(:,1), S.SeedsUpper(:,2), '.', ...
            'Color', rgb, 'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    return;
end

if norm(xy(end,:) - xy(1,:)) > 1e-6
    xy = [xy; xy(1,:)];
end

plot(ax, xy(:,1), xy(:,2), '--', 'Color', rgb, 'LineWidth', 1.2, ...
    'DisplayName', dispName);
end

function v = local_fig_visible(cfg)
v = 'on';
try
    if isfield(cfg,'io') && isfield(cfg.io,'fig_visible')
        v = cfg.io.fig_visible;
    end
catch
end
end

function v = local_cfg_get(cfg, path, defaultVal)
v = defaultVal;
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        if ~isstruct(cur) || ~isfield(cur, parts{i})
            return;
        end
        cur = cur.(parts{i});
    end
    v = cur;
catch
    v = defaultVal;
end
end

function p = local_saved_png_path(outdir, baseName, cfg)
subdir = '';
if isfield(cfg,'io') && isfield(cfg.io,'fig_subdir') && ~isempty(cfg.io.fig_subdir)
    subdir = cfg.io.fig_subdir;
end

if isempty(subdir)
    p = fullfile(outdir, [baseName '.png']);
else
    p = fullfile(outdir, subdir, [baseName '.png']);
end
end
