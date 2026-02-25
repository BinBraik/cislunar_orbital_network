function rs4_diffcorr_visualize(T, Tc, SA, SB, cfg, outdir, tag)
%RS4_DIFFCORR_VISUALIZE  Side-by-side before/after differential correction.
%
% Creates a 1×2 figure:
%   Left  panel  — original arcs (T)  showing miss distances to voxel centre
%   Right panel  — corrected arcs (Tc) showing exact patch-point connection
%
% Inputs
%   T      : struct from rs4_voxel_traj_extract     (before DC)
%   Tc     : struct from rs4_diffcorr               (after  DC)
%   SA, SB : family structs
%   cfg    : config struct
%   outdir : output directory for saving figure
%   tag    : string tag for filename / title

if nargin < 6 || isempty(outdir), outdir = pwd;      end
if nargin < 7 || isempty(tag),    tag    = 'diffcorr'; end
if ~exist(outdir, 'dir'), mkdir(outdir); end

grid3 = SA.grid3;
CJbg  = min(SA.CJ, SB.CJ);
mu    = SA.mu;

col_A     = [0.10 0.40 0.85];   % blue  — arc A
col_B     = [0.85 0.15 0.10];   % red   — arc B
col_patch = [0.10 0.65 0.20];   % green — exact patch point

% =========================================================================
% Figure layout
% =========================================================================
fig = figure('Color','w', ...
             'Name',  ['Diff Correction Comparison ' tag], ...
             'Visible', local_fig_visible(cfg), ...
             'Units','pixels', 'Position',[60 60 1400 680]);

ax1 = subplot(1,2,1, 'Parent', fig);
ax2 = subplot(1,2,2, 'Parent', fig);

% =========================================================================
% Left panel — BEFORE
% =========================================================================
local_draw_background(ax1, CJbg, mu, SA, SB);

% --- departure arc A ---
plot(ax1, T.XA(:,1), T.XA(:,2), '-', ...
    'Color', col_A, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Arc A  (%.1f d)', T.tof_A_days));

% --- arrival arc B (R-transform already in T) ---
plot(ax1, T.x_B, T.y_B, '-', ...
    'Color', col_B, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Arc B  (%.1f d)', T.tof_B_days));

% --- departure seeds ---
local_plot_seed(ax1, T.seed_A(1),  T.seed_A(2),  col_A);
local_plot_seed(ax1, T.IC_B_frs(1), -T.IC_B_frs(2), col_B);

% --- voxel rectangle ---
rectangle('Parent', ax1, ...
    'Position', [T.xc - grid3.dx/2, T.yc - grid3.dy/2, grid3.dx, grid3.dy], ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2, 'FaceColor', 'none');
plot(ax1, T.xc, T.yc, 'k+', 'MarkerSize', 6, 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

% --- closest-point markers ---
xA_star = T.XA(T.i_star,1);  yA_star = T.XA(T.i_star,2);
xB_star = T.x_B(T.j_star);   yB_star = T.y_B(T.j_star);

plot(ax1, xA_star, yA_star, 'o', ...
    'Color', col_A, 'MarkerFaceColor', col_A, 'MarkerSize', 7, ...
    'HandleVisibility', 'off');
plot(ax1, xB_star, yB_star, 'o', ...
    'Color', col_B, 'MarkerFaceColor', col_B, 'MarkerSize', 7, ...
    'HandleVisibility', 'off');

% --- miss-distance dashed lines ---
plot(ax1, [xA_star, T.xc], [yA_star, T.yc], '--', ...
    'Color', col_A, 'LineWidth', 0.9, 'HandleVisibility', 'off');
plot(ax1, [xB_star, T.xc], [yB_star, T.yc], '--', ...
    'Color', col_B, 'LineWidth', 0.9, 'HandleVisibility', 'off');

% --- miss-distance labels ---
midA = [(xA_star+T.xc)/2, (yA_star+T.yc)/2];
midB = [(xB_star+T.xc)/2, (yB_star+T.yc)/2];
text(ax1, midA(1), midA(2), sprintf('d_A=%.4f', T.dA_nd), ...
    'FontSize', 7, 'Color', col_A, 'HorizontalAlignment', 'left');
text(ax1, midB(1), midB(2), sprintf('d_B=%.4f', T.dB_nd), ...
    'FontSize', 7, 'Color', col_B, 'HorizontalAlignment', 'left');

% --- DV annotation ---
local_dv_box(ax1, T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps, T.DV_total_true_mps);

legend(ax1, 'Location', 'best', 'FontSize', 7);
title(ax1, {'\bf Before differential correction', ...
    sprintf('DV_{turn,A}=%.1f + DV_{patch}=%.1f + DV_{turn,B}=%.1f = %.1f m/s', ...
        T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps, T.DV_total_true_mps), ...
    sprintf('\\Delta\\theta=%.3f\\circ   miss: d_A=%.4f, d_B=%.4f nd', ...
        rad2deg(T.delta_th_rad), T.dA_nd, T.dB_nd)}, ...
    'Interpreter','tex', 'FontSize', 8);
xlabel(ax1,'x (nd)');  ylabel(ax1,'y (nd)');

% =========================================================================
% Right panel — AFTER
% =========================================================================
local_draw_background(ax2, CJbg, mu, SA, SB);

% --- corrected arc A ---
plot(ax2, Tc.XA(:,1), Tc.XA(:,2), '-', ...
    'Color', col_A, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Arc A  (%.1f d)', Tc.tof_A_days));

% --- corrected arc B ---
plot(ax2, Tc.x_B, Tc.y_B, '-', ...
    'Color', col_B, 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Arc B  (%.1f d)', Tc.tof_B_days));

% --- departure seeds ---
local_plot_seed(ax2, Tc.seed_A(1),   Tc.seed_A(2),   col_A);
local_plot_seed(ax2, Tc.IC_B_frs(1), -Tc.IC_B_frs(2), col_B);

% --- exact patch point (green star — arcs meet here) ---
plot(ax2, Tc.xp, Tc.yp, 'p', ...
    'Color', col_patch, 'MarkerFaceColor', col_patch, 'MarkerSize', 14, ...
    'DisplayName', sprintf('Patch point  (DV_{patch}=%.2f m/s)', Tc.DV_patch_mps));

% --- residual vector (tiny, should be nearly invisible) ---
plot(ax2, [Tc.XA(end,1), Tc.x_B(end)], [Tc.XA(end,2), Tc.y_B(end)], ...
    'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% --- DV annotation ---
local_dv_box(ax2, Tc.DV_turn_A_mps, Tc.DV_patch_mps, Tc.DV_turn_B_mps, Tc.DV_total_mps);

legend(ax2, 'Location', 'best', 'FontSize', 7);
title(ax2, {'\bf After differential correction', ...
    sprintf('DV_{turn,A}=%.1f + DV_{patch}=%.1f + DV_{turn,B}=%.1f = %.1f m/s', ...
        Tc.DV_turn_A_mps, Tc.DV_patch_mps, Tc.DV_turn_B_mps, Tc.DV_total_mps), ...
    sprintf('\\Delta\\theta=%.4f\\circ   ||r||=%.2e   exitflag=%d', ...
        rad2deg(Tc.delta_th_rad), norm(Tc.r_final), Tc.exitflag)}, ...
    'Interpreter','tex', 'FontSize', 8);
xlabel(ax2,'x (nd)');  ylabel(ax2,'y (nd)');

% =========================================================================
% Match axis limits across both panels
% =========================================================================
local_match_axes(ax1, ax2);

% =========================================================================
% Super-title
% =========================================================================
sgtitle(fig, sprintf('%s  \\rightarrow  %s  |  DV: %.1f m/s  \\rightarrow  %.1f m/s  (\\Delta=%.1f m/s)', ...
    SA.name, SB.name, T.DV_total_true_mps, Tc.DV_total_mps, ...
    T.DV_total_true_mps - Tc.DV_total_mps), ...
    'Interpreter','tex', 'FontSize', 10, 'FontWeight','bold');

% =========================================================================
% Save
% =========================================================================
safeTag = regexprep(tag, '[^A-Za-z0-9_]', '_');
rs3_io_save_figure(fig, outdir, ['rs4_' safeTag '_diffcorr_compare'], cfg);
end

% =========================================================================
% Local helpers
% =========================================================================

function local_draw_background(ax, CJbg, mu, SA, SB)
% Draw cislunar background + both periodic orbits onto axes ax.
% PO_B is R-transformed (y → -y) to match the BRS frame in which arc B
% and its seed marker are plotted.
    rs3_core_plot_cislunar_background(CJbg, mu, ax);
    set(ax.Children, 'HandleVisibility', 'off');
    hold(ax, 'on');
    axis(ax, 'equal');
    grid(ax, 'on');
    local_plot_po(ax, SA, [0.10 0.40 0.85], 'PO_A');
    local_plot_po(ax, SB, [0.85 0.15 0.10], 'PO_B', true);  % R-transform: y → -y
end

% -------------------------------------------------------------------------

function local_plot_po(ax, S, rgb, dispName, y_flip)
    if nargin < 5, y_flip = false; end
    if isfield(S,'PO_xy') && ~isempty(S.PO_xy)
        xy = S.PO_xy;
    elseif isfield(S,'Xpo') && ~isempty(S.Xpo)
        xy = S.Xpo(:,1:2);
    else
        return;
    end
    if y_flip
        xy(:,2) = -xy(:,2);
    end
    if norm(xy(end,:) - xy(1,:)) > 1e-6
        xy = [xy; xy(1,:)];
    end
    plot(ax, xy(:,1), xy(:,2), '--', 'Color', rgb, 'LineWidth', 1.0, ...
        'DisplayName', dispName);
end

% -------------------------------------------------------------------------

function local_plot_seed(ax, x, y, rgb)
% White-filled square marker at departure point.
    plot(ax, x, y, 's', 'Color', rgb, 'MarkerFaceColor','w', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility','off');
end

% -------------------------------------------------------------------------

function local_dv_box(ax, dvA, dvP, dvB, dvTot)
% Annotation text box with DV breakdown in lower-left corner.
    str = sprintf('DV_{A}=%.1f\nDV_P=%.1f\nDV_{B}=%.1f\nTotal=%.1f m/s', ...
        dvA, dvP, dvB, dvTot);
    annotation_ax = ax;
    xl = xlim(annotation_ax);
    yl = ylim(annotation_ax);
    text(annotation_ax, xl(1) + 0.02*(xl(2)-xl(1)), ...
                        yl(1) + 0.04*(yl(2)-yl(1)), ...
         str, 'FontSize', 7, 'VerticalAlignment','bottom', ...
         'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.6 0.6 0.6], ...
         'HandleVisibility','off');
end

% -------------------------------------------------------------------------

function local_match_axes(ax1, ax2)
% Set both axes to the same limits (union of both).
    xl1 = xlim(ax1);  yl1 = ylim(ax1);
    xl2 = xlim(ax2);  yl2 = ylim(ax2);
    xl  = [min(xl1(1), xl2(1)), max(xl1(2), xl2(2))];
    yl  = [min(yl1(1), yl2(1)), max(yl1(2), yl2(2))];
    xlim(ax1, xl);  ylim(ax1, yl);
    xlim(ax2, xl);  ylim(ax2, yl);
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
