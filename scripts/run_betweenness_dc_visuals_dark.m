%% RUN_BETWEENNESS_DC_VISUALS_DARK
%
% Dark-theme, presentation-style redraw of the per-example corrected relay
% figure produced by run_betweenness_dc, kept in a separate script so the
% original (light) figures are untouched.
%
% Does NOT re-run differential correction: reads the already-computed
% *_dc.mat files (Tc1, Tc2, coast arc, DV values) produced by
% run_betweenness_dc and just re-integrates each family's periodic orbit
% for the background plot — same lightweight approach as
% run_betweenness_video.
%
% Adds a "relay path" node-diagram inset (top-left) showing
% origin -> bridge -> destination, matching the paper's presentation style.
%
% Output per example: <tag>_dc_dark.png / .fig  (alongside the light *_dc.*)

clear; clc;

% ── USER KNOBS ────────────────────────────────────────────────────────────────
EXPLAINER_DIR = fullfile(pwd, 'betweenness_explainer');

GRID_DX     = 0.001;
GRID_DY     = 0.001;
GRID_DTHETA = deg2rad(1);

FIG_VISIBLE = 'off';
% ── END USER KNOBS ────────────────────────────────────────────────────────────

setup();

ode_opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
grid3    = struct('dx', GRID_DX, 'dy', GRID_DY, 'dtheta', GRID_DTHETA);
cfg      = atlas_cfg_defaults();
TU_days  = cfg.units.TU_days;

% Bright palette (matches run_betweenness_video_dark for a consistent look
% across the DC figure and the video)
c_A  = [0.30 0.55 1.00];   % origin — bright blue
c_Br = [0.35 0.85 0.45];   % bridge — bright green
c_B  = [0.95 0.35 0.20];   % dest   — bright red-orange

BG      = [0.0  0.0  0.0 ];
TXTC    = [0.92 0.94 0.99];
LEGBG   = [0.09 0.10 0.14];
LEGEDGE = [0.40 0.44 0.52];
EDGEC   = [0.90 0.92 0.97];   % light marker-edge for contrast on black

% ── Scan for *_dc.mat files ───────────────────────────────────────────────────
if ~exist(EXPLAINER_DIR, 'dir')
    error('[betweenness_dc_dark] Explainer folder not found:\n  %s', EXPLAINER_DIR);
end
sub_entries = dir(EXPLAINER_DIR);
sub_entries = sub_entries([sub_entries.isdir] & ~startsWith({sub_entries.name}, '.'));

mat_files = struct('folder', {}, 'name', {});
for k = 1:numel(sub_entries)
    hits = dir(fullfile(sub_entries(k).folder, sub_entries(k).name, '*_dc.mat'));
    for h = 1:numel(hits)
        mat_files(end+1) = struct('folder', hits(h).folder, 'name', hits(h).name); %#ok<AGROW>
    end
end
if isempty(mat_files)
    error('[betweenness_dc_dark] No *_dc.mat files found under:\n  %s', EXPLAINER_DIR);
end
fprintf('[betweenness_dc_dark] Found %d dc mat(s).\n\n', numel(mat_files));

% ── Process each example ──────────────────────────────────────────────────────
for fi = 1:numel(mat_files)
    mat_path = fullfile(mat_files(fi).folder, mat_files(fi).name);
    ex_dir   = mat_files(fi).folder;
    [~, base_name] = fileparts(mat_files(fi).name);

    fprintf('============================================================\n');
    fprintf('[%d/%d] %s\n', fi, numel(mat_files), mat_files(fi).name);
    fprintf('============================================================\n');

    D = load(mat_path);
    famA = D.famA; famBr = D.famBr; famB = D.famB;
    Tc1 = D.Tc1; Tc2 = D.Tc2;

    % ── Minimal family structs + re-integrated POs (dc.mat stores no dense PO) ─
    SA  = local_build_fam(famA,  grid3);
    SBr = local_build_fam(famBr, grid3);
    SB  = local_build_fam(famB,  grid3);
    mu  = SA.mu;

    SA.Xpo  = local_integrate_po(SA,  ode_opts);
    SBr.Xpo = local_integrate_po(SBr, ode_opts);
    SB.Xpo  = local_integrate_po(SB,  ode_opts);

    CJbg = min([SA.CJ, SBr.CJ, SB.CJ]);

    fig = figure('Color', BG, ...
        'Name',    sprintf('DC (dark)  %s  via  %s  to  %s', famA, famBr, famB), ...
        'Visible', FIG_VISIBLE, 'InvertHardcopy', 'off', ...
        'Units',   'pixels', 'Position', [80 80 1050 820]);
    ax = axes('Parent', fig);
    cr3bp_plot_background_dark(CJbg, mu, ax);
    set(ax.Children, 'HandleVisibility', 'off');
    hold(ax, 'on');
    axis(ax, 'equal');
    grid(ax, 'on');

    % Periodic orbits (dashed + arrows)
    local_plot_po_with_arrows(ax, SA,  c_A,  sprintf('Origin: %s',  famA));
    local_plot_po_with_arrows(ax, SBr, c_Br, sprintf('Bridge: %s',  famBr));
    local_plot_po_with_arrows(ax, SB,  c_B,  sprintf('Dest: %s',    famB));

    % Transfer A — corrected arcs
    plot(ax, Tc1.XA(:,1), Tc1.XA(:,2), '-', 'Color', c_A, 'LineWidth', 2.2, ...
        'DisplayName', sprintf('Transfer A DC (%.0f+%.0f d)', Tc1.tof_A_days, Tc1.tof_B_days));
    plot(ax, Tc1.x_B(end:-1:1), Tc1.y_B(end:-1:1), '-', 'Color', c_A, 'LineWidth', 2.2, ...
        'HandleVisibility', 'off');

    % Coast arc on bridge
    if ~isempty(D.coast_arc_dc) && size(D.coast_arc_dc, 1) > 1
        plot(ax, D.coast_arc_dc(:,1), D.coast_arc_dc(:,2), '-', 'Color', c_Br, ...
            'LineWidth', 4.0, 'DisplayName', sprintf('Bridge coast (%.0f d)', D.coast_time_dc*TU_days));
    end

    % Transfer B — corrected arcs
    plot(ax, Tc2.XA(:,1), Tc2.XA(:,2), '-', 'Color', c_B, 'LineWidth', 2.2, ...
        'DisplayName', sprintf('Transfer B DC (%.0f+%.0f d)', Tc2.tof_A_days, Tc2.tof_B_days));
    plot(ax, Tc2.x_B(end:-1:1), Tc2.y_B(end:-1:1), '-', 'Color', c_B, 'LineWidth', 2.2, ...
        'HandleVisibility', 'off');

    % Patch markers
    plot(ax, Tc1.xp, Tc1.yp, 'p', 'Color', EDGEC, 'MarkerFaceColor',[1.00 0.90 0.10], ...
        'MarkerSize',14, 'LineWidth',1.2, 'DisplayName','\DeltaV patch');
    plot(ax, Tc2.xp, Tc2.yp, 'p', 'Color', EDGEC, 'MarkerFaceColor',[1.00 0.90 0.10], ...
        'MarkerSize',14, 'LineWidth',1.2, 'HandleVisibility','off');

    % Event markers
    plot(ax, Tc1.XA(1,1), Tc1.XA(1,2), 'o', 'Color', EDGEC, ...
        'MarkerFaceColor', c_A,  'MarkerSize',11, 'LineWidth',1.5, 'DisplayName','Depart A');
    plot(ax, Tc1.x_B(1),  Tc1.y_B(1),  'd', 'Color', EDGEC, ...
        'MarkerFaceColor', [1 1 1], 'MarkerSize',11, 'LineWidth',2.0, 'DisplayName','Bridge arrival');
    plot(ax, Tc2.XA(1,1), Tc2.XA(1,2), 'd', 'Color', EDGEC, ...
        'MarkerFaceColor', c_Br, 'MarkerSize',11, 'LineWidth',1.5, 'DisplayName','Bridge departure');
    plot(ax, Tc2.x_B(1),  Tc2.y_B(1),  'o', 'Color', EDGEC, ...
        'MarkerFaceColor', c_B,  'MarkerSize',11, 'LineWidth',1.5, 'DisplayName','Arrive B');

    dv_via_dc  = D.dv_via_dc_mps;
    savings_dc = D.dv_direct_mps - dv_via_dc;
    title(ax, { ...
        sprintf('%s  \\rightarrow  %s  \\rightarrow  %s', famA, famBr, famB), ...
        sprintf('Via DC: %.0f + %.0f = \\bf%.0f m/s\\rm   |   Direct (proxy): %.0f m/s   |   Savings: \\bf%.0f m/s (%.0f%%)', ...
            D.dv_leg1_dc_mps, D.dv_leg2_dc_mps, dv_via_dc, D.dv_direct_mps, savings_dc, 100*savings_dc/D.dv_direct_mps)}, ...
        'Interpreter','tex', 'FontSize', 10, 'Color', TXTC);
    xlabel(ax, 'x [nd]', 'Color', TXTC);
    ylabel(ax, 'y [nd]', 'Color', TXTC);
    lg = legend(ax, 'Location','best', 'FontSize', 8);
    set(lg, 'TextColor', TXTC, 'Color', LEGBG, 'EdgeColor', LEGEDGE);

    % Relay-path node-diagram inset (top-left)
    ax_node = axes('Parent', fig, 'Position', [0.06 0.78 0.32 0.18]);
    local_draw_relay_inset(ax_node, famA, famBr, famB, c_A, c_Br, c_B);

    ex_tag_dark = [strrep(base_name, '_dc', ''), '_dc_dark'];
    local_save_dark(fig, ex_dir, ex_tag_dark, BG);
    if strcmpi(FIG_VISIBLE, 'off'), close(fig); end
    fprintf('  Figure saved: %s\n', ex_tag_dark);
    fprintf('[%d/%d] Done.\n\n', fi, numel(mat_files));
end

fprintf('[betweenness_dc_dark] All examples complete.\n');

% =============================================================================
% LOCAL HELPERS
% =============================================================================

function S = local_build_fam(fam_name, grid3)
    [mu, CJ, Tf_PO, X0] = cr3bp_family_ic(fam_name);
    S.name  = fam_name;
    S.mu    = mu;
    S.CJ    = CJ;
    S.Tf_PO = Tf_PO;
    S.X0    = X0;
    S.grid3 = grid3;
end

% ─────────────────────────────────────────────────────────────────────────────

function Xpo = local_integrate_po(S, ode_opts)
    sol = ode113(@(t,X) cr3bp_reduced_ode(t, X, S.CJ, S.mu, false), ...
                 [0, S.Tf_PO], S.X0, ode_opts);
    td  = linspace(0, min(S.Tf_PO, sol.x(end)), 2001)';
    Xpo = deval(sol, td)';
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

% ─────────────────────────────────────────────────────────────────────────────

function [l1, l2] = local_two_line_name(fam_name)
%LOCAL_TWO_LINE_NAME  Two-line abbreviation for the relay-path inset labels.
    map = { ...
        'Lyapunov L1',            'Lyap.',  'L1';     ...
        'Lyapunov L2',            'Lyap.',  'L2';     ...
        'Cycler 21',              'Cycler', '2:1';    ...
        'Cycler 11a',             'Cycler', '(1,1)a'; ...
        'Cycler 11b',             'Cycler', '(1,1)b'; ...
        'Cycler 32',              'Cycler', '3:2';    ...
        'Resonant 2to1 Stable',   'R2:1',   'S';      ...
        'Resonant 2to1 Unstable', 'R2:1',   'U';      ...
        'Resonant 3to1 Stable',   'R3:1',   'S';      ...
        'Resonant 3to1 Unstable', 'R3:1',   'U';      ...
        'Resonant 5to2 Stable',   'R5:2',   'S';      ...
        'Resonant 5to2 Unstable', 'R5:2',   'U';      ...
        'Distant Prograde Orbit', 'DPO',    '';       ...
    };
    idx = find(strcmp(map(:,1), fam_name), 1);
    if ~isempty(idx)
        l1 = map{idx,2}; l2 = map{idx,3};
    else
        l1 = fam_name; l2 = '';
    end
end

% ─────────────────────────────────────────────────────────────────────────────

function local_draw_relay_inset(ax, famA, famBr, famB, c_A, c_Br, c_B)
%LOCAL_DRAW_RELAY_INSET  Static "relay path" node diagram (origin->bridge->dest).
    hold(ax, 'on');
    set(ax, 'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off', ...
        'Color', 'none', 'XTick', [], 'YTick', [], ...
        'DataAspectRatio', [1 1 1], 'DataAspectRatioMode', 'manual');

    % Panel background
    fill(ax, [0.02 0.98 0.98 0.02], [0.04 0.04 0.96 0.96], [0.07 0.08 0.11], ...
        'EdgeColor', [0.40 0.44 0.52], 'LineWidth', 1.0);
    text(ax, 0.5, 0.87, 'relay path', 'Color', [0.92 0.94 0.99], ...
        'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    cx = [0.18, 0.50, 0.82];
    cy = [0.36, 0.36, 0.36];
    r  = 0.16;
    th = linspace(0, 2*pi, 80);
    colors = {c_A, c_Br, c_B};
    [n1A,  n2A]  = local_two_line_name(famA);
    [n1Br, n2Br] = local_two_line_name(famBr);
    [n1B,  n2B]  = local_two_line_name(famB);
    line1 = {n1A, n1Br, n1B};
    line2 = {n2A, n2Br, n2B};

    gap = r + 0.02;
    quiver(ax, cx(1)+gap, cy(1), cx(2)-cx(1)-2*gap, 0, 0, ...
        'Color', [0.85 0.87 0.92], 'LineWidth', 1.6, 'MaxHeadSize', 0.8, 'AutoScale', 'off');
    quiver(ax, cx(2)+gap, cy(2), cx(3)-cx(2)-2*gap, 0, 0, ...
        'Color', [0.85 0.87 0.92], 'LineWidth', 1.6, 'MaxHeadSize', 0.8, 'AutoScale', 'off');

    for k = 1:3
        fill(ax, cx(k) + r*cos(th), cy(k) + r*sin(th), colors{k}, ...
            'EdgeColor', [0.85 0.87 0.92], 'LineWidth', 1.3);
        if isempty(line2{k})
            text(ax, cx(k), cy(k), line1{k}, 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', 'FontSize', 7.5, 'FontWeight','bold', ...
                'Color','w', 'Interpreter','none');
        else
            text(ax, cx(k), cy(k)+0.05, line1{k}, 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', 'FontSize', 7.5, 'FontWeight','bold', ...
                'Color','w', 'Interpreter','none');
            text(ax, cx(k), cy(k)-0.05, line2{k}, 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', 'FontSize', 7.5, 'FontWeight','bold', ...
                'Color','w', 'Interpreter','none');
        end
    end
end

% ─────────────────────────────────────────────────────────────────────────────

function local_save_dark(fig, outdir, baseName, BG)
%LOCAL_SAVE_DARK  Save .png + .fig preserving the dark background.
    set(fig, 'Color', BG, 'InvertHardcopy', 'off');
    pngPath = fullfile(outdir, [baseName '.png']);
    try
        exportgraphics(fig, pngPath, 'Resolution', 220, 'BackgroundColor', 'current');
    catch ME
        warning('[betweenness_dc_dark] Failed to save PNG (%s): %s', pngPath, ME.message);
    end
    figPath = fullfile(outdir, [baseName '.fig']);
    try
        savefig(fig, figPath);
    catch ME
        warning('[betweenness_dc_dark] Failed to save FIG (%s): %s', figPath, ME.message);
    end
end
