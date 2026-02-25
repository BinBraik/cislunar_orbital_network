function rs4_voxel_traj_visualize_single(T, SA, SB, cfg, outdir, tag)
%RS4_VOXEL_TRAJ_VISUALIZE_SINGLE  Plot one transfer trajectory pair.
%
% Shows the cislunar background, both periodic orbits (dashed), the
% departure arc A (blue) and arrival arc B (red, after R-transform), the
% target voxel as a rectangle, and the miss-distance lines from each arc's
% closest point to the voxel center.
%
% Inputs
%   T      : output struct from rs4_voxel_traj_extract
%   SA, SB : family structs
%   cfg    : config struct
%   outdir : output directory for saving figure
%   tag    : string tag for filename and title

if nargin < 5 || isempty(outdir), outdir = pwd;    end
if nargin < 6 || isempty(tag),    tag    = 'traj'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end

grid3  = SA.grid3;
CJbg   = min(SA.CJ, SB.CJ);
mu     = SA.mu;

% ---- figure ----
fig = figure('Color','w', 'Name', ['Transfer Trajectory ' tag], ...
             'Visible', local_fig_visible(cfg), ...
             'Units','pixels', 'Position',[100 100 900 720]);
ax = axes('Parent', fig);
rs3_core_plot_cislunar_background(CJbg, mu, ax);
set(ax.Children, 'HandleVisibility', 'off');
hold(ax, 'on');
axis(ax, 'equal');
grid(ax, 'on');

% ---- Periodic orbit A (blue dashed) ----
local_plot_po(ax, SA, [0.15 0.45 0.80], 'Origin PO');

% ---- Periodic orbit B (red dashed, R-transformed to match BRS frame) ----
local_plot_po(ax, SB, [0.80 0.20 0.15], 'Target PO', true);  % y → -y

% ---- Departure arc A (blue solid) ----
hA = plot(ax, T.XA(:,1), T.XA(:,2), '-', ...
    'Color', [0.10 0.40 0.85], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Departure arc A  (%.1f days)', T.tof_A_days));

% ---- Arrival arc B (red solid, after R-transform) ----
hB = plot(ax, T.x_B, T.y_B, '-', ...
    'Color', [0.85 0.15 0.10], 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Arrival arc B  (%.1f days)', T.tof_B_days));

% ---- Voxel rectangle ----
dx = grid3.dx;
dy = grid3.dy;
rectangle('Parent', ax, ...
    'Position', [T.xc - dx/2, T.yc - dy/2, dx, dy], ...
    'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.2, 'FaceColor', 'none');
plot(ax, T.xc, T.yc, 'k+', 'MarkerSize', 6, 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

% ---- Closest point markers ----
% A closest point (filled blue circle)
plot(ax, T.XA(T.i_star, 1), T.XA(T.i_star, 2), 'o', ...
    'Color', [0.10 0.40 0.85], 'MarkerFaceColor', [0.10 0.40 0.85], ...
    'MarkerSize', 7, 'HandleVisibility', 'off');
% B closest point (filled red circle)
plot(ax, T.x_B(T.j_star), T.y_B(T.j_star), 'o', ...
    'Color', [0.85 0.15 0.10], 'MarkerFaceColor', [0.85 0.15 0.10], ...
    'MarkerSize', 7, 'HandleVisibility', 'off');

% ---- Miss-distance lines (dashed) ----
xA_star = T.XA(T.i_star, 1);  yA_star = T.XA(T.i_star, 2);
xB_star = T.x_B(T.j_star);    yB_star = T.y_B(T.j_star);

plot(ax, [xA_star, T.xc], [yA_star, T.yc], '--', ...
    'Color', [0.10 0.40 0.85], 'LineWidth', 0.9, 'HandleVisibility', 'off');
plot(ax, [xB_star, T.xc], [yB_star, T.yc], '--', ...
    'Color', [0.85 0.15 0.10], 'LineWidth', 0.9, 'HandleVisibility', 'off');

% ---- Miss-distance annotations ----
midA = [(xA_star + T.xc)/2, (yA_star + T.yc)/2];
midB = [(xB_star + T.xc)/2, (yB_star + T.yc)/2];
text(ax, midA(1), midA(2), sprintf('d_A=%.4f', T.dA_nd), ...
    'FontSize', 7, 'Color', [0.10 0.40 0.85], ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
text(ax, midB(1), midB(2), sprintf('d_B=%.4f', T.dB_nd), ...
    'FontSize', 7, 'Color', [0.85 0.15 0.10], ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% ---- Seed markers ----
plot(ax, T.seed_A(1), T.seed_A(2), 's', ...
    'Color', [0.10 0.40 0.85], 'MarkerFaceColor', 'w', ...
    'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');

% B seed is the R-transform of IC_B_frs start: (x, -y)
xB_seed = T.IC_B_frs(1);
yB_seed = -T.IC_B_frs(2);     % R-transform: y → -y
plot(ax, xB_seed, yB_seed, 's', ...
    'Color', [0.85 0.15 0.10], 'MarkerFaceColor', 'w', ...
    'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');

% ---- Legend ----
legend(ax, 'Location', 'best', 'FontSize', 8);

% ---- Title ----
famA = SA.name;
famB = SB.name;
title(ax, {sprintf('%s  \\rightarrow  %s  |  voxel #%d', famA, famB, T.vid), ...
           sprintf('DV_{turn,A}=%.1f  +  DV_{patch}=%.1f  +  DV_{turn,B}=%.1f  =  %.1f m/s', ...
               T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps, T.DV_total_true_mps), ...
           sprintf('\\Delta\\theta=%.3f\\circ   miss: d_A=%.4f nd, d_B=%.4f nd   (proxy was %.1f m/s)', ...
               rad2deg(T.delta_th_rad), T.dA_nd, T.dB_nd, T.DV_proxy_mps)}, ...
    'Interpreter', 'tex', 'FontSize', 9);

xlabel(ax, 'x (nd)');
ylabel(ax, 'y (nd)');

% ---- Save ----
safeTag = regexprep(tag, '[^A-Za-z0-9_]', '_');
rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_traj_single'], cfg);
end

% =========================================================================
% Local helpers
% =========================================================================

function local_plot_po(ax, S, rgb, dispName, y_flip)
% Plot the periodic orbit from S.PO_xy or S.Xpo (whichever is available).
% y_flip=true applies the R-transform (y → -y) to display in BRS frame.
if nargin < 5, y_flip = false; end
hasPOxy = isfield(S, 'PO_xy') && ~isempty(S.PO_xy);
hasXpo  = isfield(S, 'Xpo')   && ~isempty(S.Xpo);

if hasPOxy
    xy = S.PO_xy;
elseif hasXpo
    xy = S.Xpo(:, 1:2);
else
    % No orbit trace cached — just mark the seeds
    if isfield(S, 'SeedsUpper') && ~isempty(S.SeedsUpper)
        y_vals = S.SeedsUpper(:,2);
        if y_flip, y_vals = -y_vals; end
        plot(ax, S.SeedsUpper(:,1), y_vals, '.', ...
            'Color', rgb, 'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    return;
end

if y_flip
    xy(:,2) = -xy(:,2);
end

% Close the orbit if not already closed
if norm(xy(end,:) - xy(1,:)) > 1e-6
    xy = [xy; xy(1,:)];
end

plot(ax, xy(:,1), xy(:,2), '--', 'Color', rgb, 'LineWidth', 1.2, ...
    'DisplayName', dispName);
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
