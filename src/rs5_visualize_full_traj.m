function rs5_visualize_full_traj(T, traj, SA, SB, cfg, outdir, tag)
%RS5_VISUALIZE_FULL_TRAJ  Plot a complete transfer trajectory with directional arrows.
%
% Draws the full transfer from PO_A to the target orbit as a single coherent
% picture:
%   - Cislunar background (potential contours, primary bodies)
%   - Origin PO   (blue dashed)
%   - Target PO R-transformed (red dashed, y → -y)
%   - Arc A (blue solid) from PO_A to patch point
%   - Arc B reversed (red solid) from patch toward PO_B_R
%   - Departure arrow (quiver) at the departure seed on PO_A
%   - Arrival   arrow (quiver) at the end of the approach arc near PO_B_R
%   - Patch point marker (gold star)
%   - DV budget + TOF in title
%
% Inputs
%   T    : T or Tc struct (rs4_voxel_traj_extract / rs4_diffcorr)
%   traj : full-trajectory struct from rs5_build_full_traj(T)
%   SA   : origin family struct
%   SB   : target family struct
%   cfg  : config struct
%   outdir : output directory
%   tag    : filename / title tag

if nargin < 6 || isempty(outdir), outdir = pwd;    end
if nargin < 7 || isempty(tag),    tag    = 'traj'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end

relTol = local_cfg_get(cfg, 'propag.relTol', 1e-9);
absTol = local_cfg_get(cfg, 'propag.absTol', 1e-9);
SA = local_ensure_xpo(SA, relTol, absTol, 1001);
SB = local_ensure_xpo(SB, relTol, absTol, 1001);

mu   = SA.mu;
CJbg = min(SA.CJ, SB.CJ);

clr_A     = [0.10 0.40 0.85];   % blue  — departure / arc A
clr_B     = [0.85 0.15 0.10];   % red   — arrival  / arc B
clr_patch = [0.15 0.65 0.15];   % green — patch point

% ---- figure ----
fig = figure('Color', 'w', 'Visible', local_fig_visible(cfg), ...
             'Units', 'pixels', 'Position', [100 100 940 750]);
ax = axes('Parent', fig);
rs3_core_plot_cislunar_background(CJbg, mu, ax);
set(ax.Children, 'HandleVisibility', 'off');
hold(ax, 'on');
axis(ax, 'equal');
grid(ax, 'on');

% ---- periodic orbits ----
local_plot_po(ax, SA, clr_A, [SA.name ' (origin)'],  false);
local_plot_po(ax, SB, clr_B, [SB.name ' (target)'],  true);

% ---- arc A: departure to patch (blue) ----
nA = traj.n_A;
plot(ax, traj.x(1:nA), traj.y(1:nA), '-', ...
    'Color', clr_A, 'LineWidth', 2.2, ...
    'DisplayName', sprintf('Departure arc  (%.1f d)', traj.tof_A_days));

% ---- arc B reversed: patch toward target orbit (red) ----
plot(ax, traj.x(nA+1:end), traj.y(nA+1:end), '-', ...
    'Color', clr_B, 'LineWidth', 2.2, ...
    'DisplayName', sprintf('Arrival arc  (%.1f d)',   traj.tof_B_days));

% ---- patch point (gold star) ----
plot(ax, traj.patch_x, traj.patch_y, 'p', ...
    'Color', clr_patch, 'MarkerFaceColor', [1.0 0.85 0.0], 'MarkerSize', 14, ...
    'LineWidth', 1.2, ...
    'DisplayName', sprintf('Patch  (\\DeltaV_{patch}=%.1f m/s)', traj.DV_patch_mps));

% ---- departure seed marker (open blue square on PO_A) ----
plot(ax, traj.depart_x, traj.depart_y, 's', ...
    'Color', clr_A, 'MarkerFaceColor', 'w', 'MarkerSize', 10, ...
    'LineWidth', 1.8, 'HandleVisibility', 'off');

% ---- arrival point marker (open red square near PO_B_R) ----
plot(ax, traj.arrive_x, traj.arrive_y, 's', ...
    'Color', clr_B, 'MarkerFaceColor', 'w', 'MarkerSize', 10, ...
    'LineWidth', 1.8, 'HandleVisibility', 'off');

% ---- departure arrow (quiver at departure seed, pointing along IC_A(3)) ----
arrow_scale = 0.05;   % nondimensional length — adjust if plot looks cluttered
dep_u = arrow_scale * cos(traj.depart_th);
dep_v = arrow_scale * sin(traj.depart_th);
quiver(ax, traj.depart_x, traj.depart_y, dep_u, dep_v, 0, ...
    'Color', clr_A, 'LineWidth', 2.4, 'MaxHeadSize', 1.8, ...
    'DisplayName', 'Departure direction');

% ---- arrival arrow (quiver at arrival point, pointing in approach direction) ----
arr_u = arrow_scale * cos(traj.arrive_th);
arr_v = arrow_scale * sin(traj.arrive_th);
quiver(ax, traj.arrive_x, traj.arrive_y, arr_u, arr_v, 0, ...
    'Color', clr_B, 'LineWidth', 2.4, 'MaxHeadSize', 1.8, ...
    'DisplayName', 'Arrival direction');

% ---- convergence label (only for corrected Tc) ----
dc_str = '';
if isfield(T, 'converged')
    if T.converged
        dc_str = '  \color{green}[DC \checkmark]';
    else
        dc_str = '  \color{red}[DC \times]';
    end
end

% ---- legend, title, labels ----
legend(ax, 'Location', 'best', 'FontSize', 8);

title(ax, { ...
    sprintf('%s  \\rightarrow  %s%s', SA.name, SB.name, dc_str), ...
    sprintf('\\DeltaV_A=%.1f + \\DeltaV_{patch}=%.1f + \\DeltaV_B=%.1f = %.1f m/s', ...
        traj.DV_turn_A_mps, traj.DV_patch_mps, traj.DV_turn_B_mps, traj.DV_total_mps), ...
    sprintf('TOF = %.1f + %.1f = %.1f days', ...
        traj.tof_A_days, traj.tof_B_days, traj.tof_total_days)}, ...
    'Interpreter', 'tex', 'FontSize', 9);

xlabel(ax, 'x  (nondimensional)');
ylabel(ax, 'y  (nondimensional)');

% ---- save ----
safeTag = rs3_sanitize_fname(tag);
rs3_io_save_figure(fig, outdir, ['rs5_' safeTag], cfg);
end

% =========================================================================
% Local helpers  (mirrored from rs4_voxel_traj_visualize_single)
% =========================================================================

function local_plot_po(ax, S, rgb, dispName, y_flip)
if nargin < 5, y_flip = false; end
hasPOxy = isfield(S, 'PO_xy') && ~isempty(S.PO_xy);
hasXpo  = isfield(S, 'Xpo')   && ~isempty(S.Xpo);

if hasXpo
    xy = S.Xpo(:, 1:2);
elseif hasPOxy
    xy = S.PO_xy;
else
    if isfield(S, 'SeedsUpper') && ~isempty(S.SeedsUpper)
        y_vals = S.SeedsUpper(:,2);
        if y_flip, y_vals = -y_vals; end
        plot(ax, S.SeedsUpper(:,1), y_vals, '.', ...
            'Color', rgb, 'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    return;
end
if y_flip, xy(:,2) = -xy(:,2); end
if norm(xy(end,:) - xy(1,:)) > 1e-6, xy = [xy; xy(1,:)]; end
plot(ax, xy(:,1), xy(:,2), '--', 'Color', rgb, 'LineWidth', 1.2, ...
    'DisplayName', dispName);
end

% -------------------------------------------------------------------------

function S = local_ensure_xpo(S, relTol, absTol, N_po)
if isfield(S, 'Xpo') && ~isempty(S.Xpo) && ...
   isfield(S, 't_dense') && ~isempty(S.t_dense)
    return;
end
fprintf('[rs5_viz] Xpo not cached for "%s" — re-integrating PO ...\n', S.name);
opts    = odeset('RelTol', relTol, 'AbsTol', absTol);
solPO   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,S.CJ,S.mu,true), ...
                 [0, S.Tf_PO], S.X0, opts);
t_dense = linspace(0, S.Tf_PO, N_po)';
S.t_dense = t_dense;
S.Xpo     = deval(solPO, t_dense)';
end

% -------------------------------------------------------------------------

function v = local_fig_visible(cfg)
v = 'on';
try
    if isfield(cfg,'io') && isfield(cfg.io,'fig_visible')
        v = cfg.io.fig_visible;
    end
catch
end
end

% -------------------------------------------------------------------------

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
