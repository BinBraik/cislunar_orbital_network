%% RUN_BETWEENNESS_DC
%
% Post-processing script: applies differential correction to the per-example
% *_data.mat files produced by run_betweenness_explainer.
%
% Does NOT require the reachable-set cache or any atlas/overlap computation.
% All inputs come from:
%   1. The *_data.mat files under rs3_betweenness_explainer/
%   2. rs3_core_family_ic  — mu, CJ, Tf_PO for each family (hardcoded constants)
%
% For each mat file found:
%   a. Builds minimal family structs from rs3_core_family_ic + PO_* arrays in mat
%   b. Runs rs4_diffcorr on T1 (origin→bridge) and T2 (bridge→dest)
%   c. Recomputes coast arc using corrected bridge endpoints
%   d. Prints raw vs corrected DV comparison table
%   e. Saves <ex_dir>/<tag>_dc.mat  (Tc1, Tc2, corrected DV scalars)
%   f. Saves corrected figure

clear; clc;

% ── USER KNOBS ───────────────────────────────────────────────────────────────

% Root folder produced by run_betweenness_explainer
EXPLAINER_DIR = fullfile(pwd, 'rs3_betweenness_explainer');

% Grid spacings used in the original run (must match)
GRID_DX     = 0.001;
GRID_DY     = 0.001;
GRID_DTHETA = deg2rad(1);

% Differential correction settings
DC_TOL      = 1e-4;
DC_DISPLAY  = 'off';   % 'off' | 'iter' | 'final'
DC_MAX_ITER = 300;
DC_MAX_FEVAL = 8000;

% Figure settings
SAVE_FIGS   = true;
FIG_VISIBLE = 'on';

% ── END USER KNOBS ───────────────────────────────────────────────────────────

rs3_setup();

% ── cfg (only propag + diffcorr fields needed — no grid sweep, no atlas) ─────
cfg = rs3_cfg_defaults();
cfg.propag.relTol        = 1e-8;
cfg.propag.absTol        = 1e-8;
cfg.propag.Tmax          = pi;
cfg.grid.dx              = GRID_DX;
cfg.grid.dy              = GRID_DY;
cfg.grid.dtheta          = GRID_DTHETA;
cfg.io.save_figs         = SAVE_FIGS;
cfg.io.fig_visible       = FIG_VISIBLE;
cfg.diffcorr.tol_patch     = DC_TOL;
cfg.diffcorr.tol_converged = DC_TOL;
cfg.diffcorr.display       = DC_DISPLAY;
cfg.diffcorr.MaxIterations = DC_MAX_ITER;
cfg.diffcorr.MaxFunEvals   = DC_MAX_FEVAL;
cfg.diffcorr.N_po_dt       = 0.003;
cfg.diffcorr.N_po_min      = 1001;

grid3   = struct('dx', GRID_DX, 'dy', GRID_DY, 'dtheta', GRID_DTHETA);
TU_days = cfg.units.TU_days;

% ── Scan for *_data.mat files one level deep ──────────────────────────────────
if ~exist(EXPLAINER_DIR, 'dir')
    error('[betweenness_dc] Explainer output folder not found:\n  %s\n', EXPLAINER_DIR);
end

sub_entries = dir(EXPLAINER_DIR);
sub_entries = sub_entries([sub_entries.isdir] & ~startsWith({sub_entries.name}, '.'));

mat_files = struct('folder', {}, 'name', {});
for k = 1:numel(sub_entries)
    hits = dir(fullfile(sub_entries(k).folder, sub_entries(k).name, '*_data.mat'));
    for h = 1:numel(hits)
        mat_files(end+1) = struct('folder', hits(h).folder, 'name', hits(h).name); %#ok<AGROW>
    end
end

if isempty(mat_files)
    error('[betweenness_dc] No *_data.mat files found under:\n  %s\n', EXPLAINER_DIR);
end
fprintf('[betweenness_dc] Found %d data mat(s).\n\n', numel(mat_files));

% ── Colour palette (matches run_betweenness_explainer) ───────────────────────
PALETTE = [
    0.15 0.42 0.80;   % blue   — origin
    0.18 0.62 0.30;   % green  — bridge
    0.82 0.25 0.12;   % red    — destination
];

% ── Process each mat ──────────────────────────────────────────────────────────
for fi = 1:numel(mat_files)
    mat_path = fullfile(mat_files(fi).folder, mat_files(fi).name);
    ex_dir   = mat_files(fi).folder;

    fprintf('============================================================\n');
    fprintf('[%d/%d] %s\n', fi, numel(mat_files), mat_files(fi).name);
    fprintf('============================================================\n');

    D = load(mat_path);

    famA  = D.famA;
    famBr = D.famBr;
    famB  = D.famB;
    fprintf('  %s  →  %s  →  %s\n\n', famA, famBr, famB);

    % ── Build minimal family structs (no cache) ───────────────────────────────
    SA  = local_build_family_struct(famA,  D.PO_origin, grid3);
    SBr = local_build_family_struct(famBr, D.PO_bridge, grid3);
    SB  = local_build_family_struct(famB,  D.PO_dest,   grid3);

    % ── Differential correction — Leg 1 ──────────────────────────────────────
    fprintf('  DC Leg 1 (%s → %s) ...\n', famA, famBr);
    try
        Tc1 = rs4_diffcorr(D.T1, SA, SBr, cfg);
        fprintf('  Leg1 DC: %.1f m/s  (converged=%d, exitflag=%d)\n', ...
            Tc1.DV_total_mps, Tc1.converged, Tc1.exitflag);
    catch ME
        fprintf('  Leg1 DC FAILED: %s — skipping example.\n', ME.message);
        continue
    end

    % ── Differential correction — Leg 2 ──────────────────────────────────────
    fprintf('  DC Leg 2 (%s → %s) ...\n', famBr, famB);
    try
        Tc2 = rs4_diffcorr(D.T2, SBr, SB, cfg);
        fprintf('  Leg2 DC: %.1f m/s  (converged=%d, exitflag=%d)\n', ...
            Tc2.DV_total_mps, Tc2.converged, Tc2.exitflag);
    catch ME
        fprintf('  Leg2 DC FAILED: %s — skipping example.\n', ME.message);
        continue
    end

    % ── DV comparison ────────────────────────────────────────────────────────
    dv_l1_raw = D.T1.DV_total_true_mps;
    dv_l2_raw = D.T2.DV_total_true_mps;
    dv_l1_dc  = Tc1.DV_total_mps;
    dv_l2_dc  = Tc2.DV_total_mps;
    dv_via_raw = dv_l1_raw + dv_l2_raw;
    dv_via_dc  = dv_l1_dc  + dv_l2_dc;

    fprintf('\n  %-8s  %8s  %8s  %9s  m/s\n', '', 'Leg 1', 'Leg 2', 'Via total');
    fprintf('  %-8s  %8.1f  %8.1f  %9.1f\n', 'Raw',   dv_l1_raw, dv_l2_raw, dv_via_raw);
    fprintf('  %-8s  %8.1f  %8.1f  %9.1f\n', 'DC',    dv_l1_dc,  dv_l2_dc,  dv_via_dc);
    fprintf('  %-8s  %8.1f  %8.1f  %9.1f\n', 'Delta', dv_l1_dc-dv_l1_raw, dv_l2_dc-dv_l2_raw, dv_via_dc-dv_via_raw);
    fprintf('  Direct (proxy): %.1f m/s  |  Savings DC: %.1f m/s (%.0f%%)\n\n', ...
        D.dv_direct_mps, D.dv_direct_mps-dv_via_dc, 100*(D.dv_direct_mps-dv_via_dc)/D.dv_direct_mps);

    % ── Recompute coast arc from corrected bridge endpoints ───────────────────
    [coast_arc_dc, coast_time_dc] = local_coast_on_po(SBr, ...
        [Tc1.x_B(1); Tc1.y_B(1)], ...
        [Tc2.XA(1,1); Tc2.XA(1,2)]);

    % ── Save corrected mat ────────────────────────────────────────────────────
    [~, base_name]   = fileparts(mat_files(fi).name);
    ex_tag_dc        = [strrep(base_name, '_data', ''), '_dc'];
    dv_leg1_raw_mps  = dv_l1_raw;
    dv_leg2_raw_mps  = dv_l2_raw;
    dv_leg1_dc_mps   = dv_l1_dc;
    dv_leg2_dc_mps   = dv_l2_dc;
    dv_via_dc_mps    = dv_via_dc;
    dv_direct_mps    = D.dv_direct_mps;

    save(fullfile(ex_dir, [ex_tag_dc '.mat']), ...
        'famA', 'famBr', 'famB', ...
        'Tc1', 'Tc2', ...
        'coast_arc_dc', 'coast_time_dc', ...
        'dv_leg1_raw_mps', 'dv_leg2_raw_mps', ...
        'dv_leg1_dc_mps',  'dv_leg2_dc_mps',  'dv_via_dc_mps', 'dv_direct_mps', ...
        '-v7.3');
    fprintf('  Saved: %s.mat\n', ex_tag_dc);

    % ── Corrected figure ──────────────────────────────────────────────────────
    c_A  = PALETTE(1,:);
    c_Br = PALETTE(2,:);
    c_B  = PALETTE(3,:);

    CJbg = min([SA.CJ, SBr.CJ, SB.CJ]);
    mu   = SA.mu;

    fig = figure('Color','w', ...
        'Name',    sprintf('DC  %s  via  %s  to  %s', famA, famBr, famB), ...
        'Visible', cfg.io.fig_visible, ...
        'Units',   'pixels', 'Position', [80 80 1050 820]);
    ax = axes('Parent', fig);
    rs3_core_plot_cislunar_background(CJbg, mu, ax);
    set(ax.Children, 'HandleVisibility', 'off');
    hold(ax, 'on');
    axis(ax, 'equal');
    grid(ax, 'on');

    % Periodic orbits (dashed + arrows)
    local_plot_po_with_arrows(ax, SA,  c_A,  sprintf('Origin: %s',  famA));
    local_plot_po_with_arrows(ax, SBr, c_Br, sprintf('Bridge: %s',  famBr));
    local_plot_po_with_arrows(ax, SB,  c_B,  sprintf('Dest: %s',    famB));

    % Transfer A — corrected arcs
    plot(ax, Tc1.XA(:,1), Tc1.XA(:,2), '-', 'Color', c_A, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Transfer A DC (%.0f+%.0f d)', Tc1.tof_A_days, Tc1.tof_B_days));
    plot(ax, Tc1.x_B(end:-1:1), Tc1.y_B(end:-1:1), '-', 'Color', c_A, 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');

    % Coast arc on bridge
    if ~isempty(coast_arc_dc) && size(coast_arc_dc, 1) > 1
        plot(ax, coast_arc_dc(:,1), coast_arc_dc(:,2), '-', 'Color', c_Br, ...
            'LineWidth', 3.5, 'DisplayName', sprintf('Bridge coast (%.0f d)', coast_time_dc*TU_days));
    end

    % Transfer B — corrected arcs
    plot(ax, Tc2.XA(:,1), Tc2.XA(:,2), '-', 'Color', c_B, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Transfer B DC (%.0f+%.0f d)', Tc2.tof_A_days, Tc2.tof_B_days));
    plot(ax, Tc2.x_B(end:-1:1), Tc2.y_B(end:-1:1), '-', 'Color', c_B, 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');

    % Patch markers — exact points from Tc (no i_star needed)
    plot(ax, Tc1.xp, Tc1.yp, 'p', 'Color','k', 'MarkerFaceColor',[0.95 0.85 0.05], ...
        'MarkerSize',14, 'LineWidth',1.2, 'DisplayName','\DeltaV patch');
    plot(ax, Tc2.xp, Tc2.yp, 'p', 'Color','k', 'MarkerFaceColor',[0.95 0.85 0.05], ...
        'MarkerSize',14, 'LineWidth',1.2, 'HandleVisibility','off');

    % Event markers
    plot(ax, Tc1.XA(1,1), Tc1.XA(1,2), 'o', 'Color','k', ...
        'MarkerFaceColor', c_A,  'MarkerSize',11, 'LineWidth',1.5, 'DisplayName','Depart A');
    plot(ax, Tc1.x_B(1),  Tc1.y_B(1),  'd', 'Color','k', ...
        'MarkerFaceColor', 'w',  'MarkerSize',11, 'LineWidth',2.0, 'DisplayName','Bridge arrival');
    plot(ax, Tc2.XA(1,1), Tc2.XA(1,2), 'd', 'Color','k', ...
        'MarkerFaceColor', c_Br, 'MarkerSize',11, 'LineWidth',1.5, 'DisplayName','Bridge departure');
    plot(ax, Tc2.x_B(1),  Tc2.y_B(1),  'o', 'Color','k', ...
        'MarkerFaceColor', c_B,  'MarkerSize',11, 'LineWidth',1.5, 'DisplayName','Arrive B');

    savings_dc = D.dv_direct_mps - dv_via_dc;
    title(ax, { ...
        sprintf('%s  \\rightarrow  %s  \\rightarrow  %s', famA, famBr, famB), ...
        sprintf('Via DC: %.0f + %.0f = \\bf%.0f m/s\\rm   |   Direct (proxy): %.0f m/s   |   Savings: \\bf%.0f m/s (%.0f%%)', ...
            dv_l1_dc, dv_l2_dc, dv_via_dc, D.dv_direct_mps, savings_dc, 100*savings_dc/D.dv_direct_mps)}, ...
        'Interpreter','tex', 'FontSize', 10);
    xlabel(ax, 'x [nd]');
    ylabel(ax, 'y [nd]');
    legend(ax, 'Location','best', 'FontSize', 8);

    rs3_io_save_figure(fig, ex_dir, ex_tag_dc, cfg);
    fprintf('  Figure saved.\n');
    fprintf('[%d/%d] Done.\n\n', fi, numel(mat_files));
end

fprintf('[betweenness_dc] All examples complete.\n');

% =============================================================================
% LOCAL HELPERS
% =============================================================================

function S = local_build_family_struct(fam_name, Xpo_mat, grid3)
%LOCAL_BUILD_FAMILY_STRUCT  Minimal struct for rs4_diffcorr — no cache needed.
% mu, CJ, Tf_PO, X0 from rs3_core_family_ic (hardcoded corrected constants).
% Xpo and t_dense reconstructed from the PO array stored in the data mat.
% pp_xy / pp_th are built automatically inside rs4_diffcorr's local_ensure_xpo
% once Xpo and t_dense are present.
    [mu, CJ, Tf_PO, X0] = rs3_core_family_ic(fam_name);
    N_po      = size(Xpo_mat, 1);
    S.name    = fam_name;
    S.mu      = mu;
    S.CJ      = CJ;
    S.Tf_PO   = Tf_PO;
    S.X0      = X0;
    S.Xpo     = Xpo_mat;
    S.t_dense = linspace(0, Tf_PO, N_po)';
    S.grid3   = grid3;
end

% ─────────────────────────────────────────────────────────────────────────────

function [coast_arc, coast_time] = local_coast_on_po(SBr, p_arrive, p_depart)
%LOCAL_COAST_ON_PO  Segment of bridge PO from arrival to departure (physical).
% Identical to the version in run_betweenness_explainer.
    coast_arc  = zeros(0, 2);
    coast_time = 0;
    if ~isfield(SBr,'Xpo') || isempty(SBr.Xpo), return; end
    xy    = SBr.Xpo(:, 1:2);
    N_po  = size(xy, 1);
    dt_po = SBr.Tf_PO / max(1, N_po - 1);

    [~, i_arr] = min(hypot(xy(:,1)-p_arrive(1), xy(:,2)-p_arrive(2)));
    [~, i_dep] = min(hypot(xy(:,1)-p_depart(1), xy(:,2)-p_depart(2)));

    if i_dep >= i_arr
        coast_arc = xy(i_arr:i_dep, :);
        n_steps   = i_dep - i_arr;
    else
        coast_arc = [xy(i_arr:end, :); xy(1:i_dep, :)];
        n_steps   = (N_po - i_arr) + i_dep;
    end
    coast_time = n_steps * dt_po;
end

% ─────────────────────────────────────────────────────────────────────────────

function local_plot_po_with_arrows(ax, S, rgb, dispName)
%LOCAL_PLOT_PO_WITH_ARROWS  Dashed PO + 4 equally-spaced direction arrows.
    xy = S.Xpo(:, 1:2);
    if norm(xy(end,:) - xy(1,:)) > 1e-6
        xy = [xy; xy(1,:)];
    end
    n = size(xy, 1);
    plot(ax, xy(:,1), xy(:,2), '--', 'Color', rgb, 'LineWidth', 1.4, ...
        'DisplayName', dispName);
    N_arr   = 4;
    arr_len = 0.025;
    for k = 1:N_arr
        i1  = max(1, round(k * n / N_arr));
        i2  = min(n, i1 + max(1, round(n/30)));
        dx  = xy(i2,1) - xy(i1,1);
        dy  = xy(i2,2) - xy(i1,2);
        nrm = hypot(dx, dy);
        if nrm < 1e-8, continue; end
        quiver(ax, xy(i1,1), xy(i1,2), arr_len*dx/nrm, arr_len*dy/nrm, 0, ...
            'Color', rgb, 'LineWidth', 1.0, 'MaxHeadSize', 3, ...
            'HandleVisibility', 'off');
    end
end
