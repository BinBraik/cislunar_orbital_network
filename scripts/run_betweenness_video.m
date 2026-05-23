%% RUN_BETWEENNESS_VIDEO
%
% Creates side-by-side MP4 videos (rotating frame + ECI) for each
% *_dc.mat file produced by run_betweenness_dc.
%
% Layout:
%   Left  — synodic (rotating) frame  with cislunar background
%   Right — Earth-centered inertial (ECI) frame, dark background
%   Inset — node diagram (origin→bridge→dest) with leg highlighting
%
% Animation uses true physical time scale (frame allocation proportional
% to ND time, not arc length), so the dot moves faster near close
% approaches and slower in the outer regions.
%
% ECI transform:
%   x_ECI = (x_rot + mu)*cos(t) - y_rot*sin(t)
%   y_ECI = (x_rot + mu)*sin(t) + y_rot*cos(t)
% where t is ND time and mu is the E-M mass ratio.

clear; clc;

% ── USER KNOBS ────────────────────────────────────────────────────────────────
EXPLAINER_DIR  = fullfile(rs3_repo_root(), 'rs3_betweenness_explainer');

GRID_DX     = 0.001;
GRID_DY     = 0.001;
GRID_DTHETA = deg2rad(1);

FPS            = 30;
VIDEO_DURATION = 20;       % seconds per example video
TRAIL_SEC      = 1.5;      % visible fading trail duration [video seconds]
N_TRAIL_SEG    = 16;       % fade segments (more = smoother fade)
VIDEO_QUALITY  = 95;       % VideoWriter quality (0-100)
% ── END USER KNOBS ────────────────────────────────────────────────────────────

rs3_setup();

ode_opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
grid3    = struct('dx', GRID_DX, 'dy', GRID_DY, 'dtheta', GRID_DTHETA);

cfg = rs3_cfg_defaults();
TU_days = cfg.units.TU_days;

% Colour palette (matches betweenness explainer)
c_A  = [0.15 0.42 0.80];
c_Br = [0.18 0.62 0.30];
c_B  = [0.82 0.25 0.12];

BG_ROT = [1.00 1.00 1.00];          % white  — rotating frame background
BG_ECI = [0.05 0.05 0.12];          % dark   — ECI frame background

% ── Scan for *_dc.mat files ───────────────────────────────────────────────────
if ~exist(EXPLAINER_DIR, 'dir')
    error('[video] Explainer folder not found:\n  %s', EXPLAINER_DIR);
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
    error('[video] No *_dc.mat files found under:\n  %s', EXPLAINER_DIR);
end
fprintf('[video] Found %d dc mat(s).\n\n', numel(mat_files));

% ── Process each example ──────────────────────────────────────────────────────
for fi = 1:numel(mat_files)
    mat_path = fullfile(mat_files(fi).folder, mat_files(fi).name);
    ex_dir   = mat_files(fi).folder;
    [~, base_name] = fileparts(mat_files(fi).name);

    fprintf('=================================================================\n');
    fprintf('[%d/%d]  %s\n', fi, numel(mat_files), base_name);
    fprintf('=================================================================\n');

    D = load(mat_path);
    Tc1 = D.Tc1;  Tc2 = D.Tc2;
    coast_arc_dc  = D.coast_arc_dc;
    coast_time_dc = D.coast_time_dc;

    % ── Build minimal family structs ─────────────────────────────────────────
    SA  = local_build_fam(D.famA,  grid3);
    SBr = local_build_fam(D.famBr, grid3);
    SB  = local_build_fam(D.famB,  grid3);
    mu  = SA.mu;

    % ── Re-integrate dense arcs ───────────────────────────────────────────────
    N_frames = FPS * VIDEO_DURATION;
    N_dense  = max(3000, 20 * N_frames);
    fprintf('  Integrating arcs (N_dense=%d) ...\n', N_dense);

    % Phase 1: Leg-1 FRS  (origin → patch-1)
    sol1   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,SA.CJ,SA.mu,false), ...
                    [0, Tc1.t_A], Tc1.IC_A, ode_opts);
    t1_ev  = linspace(0, Tc1.t_A, N_dense)';
    D1     = deval(sol1, t1_ev)';            % N_dense×3 [x y th]

    % Phase 2: Leg-1 BRS reversed  (patch-1 → bridge arrival)
    sol2   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,SBr.CJ,SBr.mu,false), ...
                    [0, Tc1.t_B], Tc1.IC_B_frs, ode_opts);
    D2_frs = flipud(deval(sol2, linspace(0, Tc1.t_B, N_dense))');  % reversed → patch→bridge
    D2_x   =  D2_frs(:,1);
    D2_y   = -D2_frs(:,2);  % R-transform

    % Phase 3: Coast on bridge PO
    N_coast = max(2, size(coast_arc_dc, 1));
    D3_x = coast_arc_dc(:,1);
    D3_y = coast_arc_dc(:,2);

    % Phase 4: Leg-2 FRS  (bridge departure → patch-2)
    sol4   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,SBr.CJ,SBr.mu,false), ...
                    [0, Tc2.t_A], Tc2.IC_A, ode_opts);
    t4_ev  = linspace(0, Tc2.t_A, N_dense)';
    D4     = deval(sol4, t4_ev)';

    % Phase 5: Leg-2 BRS reversed  (patch-2 → dest arrival)
    sol5   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,SB.CJ,SB.mu,false), ...
                    [0, Tc2.t_B], Tc2.IC_B_frs, ode_opts);
    D5_frs = flipud(deval(sol5, linspace(0, Tc2.t_B, N_dense))');
    D5_x   =  D5_frs(:,1);
    D5_y   = -D5_frs(:,2);

    % ── Absolute time axis ────────────────────────────────────────────────────
    t_off2 = Tc1.t_A;
    t_off3 = t_off2 + Tc1.t_B;
    t_off4 = t_off3 + coast_time_dc;
    t_off5 = t_off4 + Tc2.t_A;

    t1_abs = t1_ev;
    t2_abs = t_off2 + linspace(0, Tc1.t_B,     N_dense)';
    t3_abs = t_off3 + linspace(0, coast_time_dc, N_coast)';
    t4_abs = t_off4 + t4_ev;
    t5_abs = t_off5 + linspace(0, Tc2.t_B,     N_dense)';

    % ── Assemble full path ────────────────────────────────────────────────────
    path_x = [D1(:,1); D2_x; D3_x; D4(:,1); D5_x];
    path_y = [D1(:,2); D2_y; D3_y; D4(:,2); D5_y];
    path_t = [t1_abs;  t2_abs; t3_abs; t4_abs; t5_abs];

    % Phase colour at each point  (1=leg1/c_A, 2=coast/c_Br, 3=leg2/c_B)
    n12 = 2*N_dense;  n3 = N_coast;  n45 = 2*N_dense;
    phase_col = [repmat(c_A,  n12, 1); repmat(c_Br, n3, 1); repmat(c_B, n45, 1)];
    phase_id  = [ones(n12,1);           2*ones(n3,1);          3*ones(n45,1)];

    % Strict monotone time (floating-point guard)
    n_total = numel(path_t);
    path_t  = path_t + (0:n_total-1)' * 1e-13 * max(path_t(end), 1);

    % ── ECI transform ─────────────────────────────────────────────────────────
    path_x_eci = (path_x + mu) .* cos(path_t) - path_y .* sin(path_t);
    path_y_eci = (path_x + mu) .* sin(path_t) + path_y .* cos(path_t);

    % ── Integrate POs for background plots ───────────────────────────────────
    fprintf('  Integrating POs ...\n');
    XpoA  = local_integrate_po(SA,  ode_opts);
    XpoBr = local_integrate_po(SBr, ode_opts);
    XpoB  = local_integrate_po(SB,  ode_opts);

    % ── Axis limits ───────────────────────────────────────────────────────────
    pad = 0.12;
    [rot_xlim, rot_ylim] = local_sq([min(path_x)-pad, max(path_x)+pad], ...
                                    [min(path_y)-pad, max(path_y)+pad]);

    eci_all_x = [path_x_eci; 0];   % include Earth
    eci_all_y = [path_y_eci; 0];
    [eci_xlim, eci_ylim] = local_sq([min(eci_all_x)-pad, max(eci_all_x)+pad], ...
                                    [min(eci_all_y)-pad, max(eci_all_y)+pad]);

    % ── Frame index table (true time scale) ───────────────────────────────────
    frame_t   = linspace(path_t(1), path_t(end), N_frames);
    frame_idx = round(interp1(path_t, (1:n_total)', frame_t, 'linear', 'extrap'));
    frame_idx = max(1, min(n_total, frame_idx));

    trail_pts = round((TRAIL_SEC / VIDEO_DURATION) * n_total);

    % ── Build figure ──────────────────────────────────────────────────────────
    fig = figure('Color', 'w', 'Visible', 'off', ...
                 'Units', 'pixels', 'Position', [0 0 1600 800]);

    ax_rot  = axes('Parent', fig, 'Position', [0.03  0.05  0.44  0.92]);
    ax_eci  = axes('Parent', fig, 'Position', [0.535 0.05  0.44  0.92]);
    ax_node = axes('Parent', fig, 'Position', [0.03  0.70  0.16  0.26]);

    % ── Static rotating frame ─────────────────────────────────────────────────
    CJbg = min([SA.CJ, SBr.CJ, SB.CJ]);
    rs3_core_plot_cislunar_background(CJbg, mu, ax_rot);
    set(ax_rot.Children, 'HandleVisibility', 'off');
    hold(ax_rot, 'on');
    set(ax_rot, 'XLim', rot_xlim, 'YLim', rot_ylim, 'DataAspectRatio', [1 1 1]);
    set(ax_rot, 'XTick', [], 'YTick', [], 'Box', 'off');

    % Earth and Moon in rotating frame (fixed)
    plot(ax_rot, -mu,   0, 'o', 'MarkerFaceColor', [0.20 0.47 0.90], ...
         'MarkerEdgeColor', 'k', 'MarkerSize', 11, 'HandleVisibility', 'off');
    plot(ax_rot, 1-mu,  0, 'o', 'MarkerFaceColor', [0.68 0.68 0.68], ...
         'MarkerEdgeColor', 'k', 'MarkerSize', 7,  'HandleVisibility', 'off');

    % POs (dashed)
    local_plot_po(ax_rot, XpoA,  c_A  * 0.7 + 0.3);
    local_plot_po(ax_rot, XpoBr, c_Br * 0.7 + 0.3);
    local_plot_po(ax_rot, XpoB,  c_B  * 0.7 + 0.3);

    % Updateable trail + dot
    h_hist_rot  = plot(ax_rot, NaN, NaN, '-',  'Color', [0.80 0.80 0.80], 'LineWidth', 0.6);
    h_seg_rot   = gobjects(N_TRAIL_SEG, 1);
    for s = 1:N_TRAIL_SEG
        a = s / N_TRAIL_SEG;
        h_seg_rot(s) = plot(ax_rot, NaN, NaN, '-', ...
            'Color', a*c_A + (1-a)*BG_ROT, 'LineWidth', 0.8 + a*2.2);
    end
    h_dot_rot = plot(ax_rot, NaN, NaN, 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

    % ── Static ECI frame ─────────────────────────────────────────────────────
    set(ax_eci, 'Color', BG_ECI, 'XColor', [0.35 0.35 0.45], 'YColor', [0.35 0.35 0.45]);
    hold(ax_eci, 'on');
    set(ax_eci, 'XLim', eci_xlim, 'YLim', eci_ylim, 'DataAspectRatio', [1 1 1]);
    set(ax_eci, 'XTick', [], 'YTick', [], 'Box', 'on');

    theta_c = linspace(0, 2*pi, 300);
    plot(ax_eci, cos(theta_c), sin(theta_c), '-', ...
         'Color', [0.25 0.28 0.38], 'LineWidth', 0.9);   % Moon orbit

    plot(ax_eci, 0, 0, 'o', ...                           % Earth
         'MarkerFaceColor', [0.20 0.45 0.90], ...
         'MarkerEdgeColor', [0.55 0.75 1.00], ...
         'MarkerSize', 14, 'LineWidth', 1.5);

    h_moon_eci = plot(ax_eci, 1, 0, 'o', ...             % Moon (moves)
         'MarkerFaceColor', [0.60 0.62 0.68], ...
         'MarkerEdgeColor', [0.80 0.82 0.88], ...
         'MarkerSize', 8,  'LineWidth', 1.2);

    h_hist_eci  = plot(ax_eci, NaN, NaN, '-', 'Color', [0.18 0.20 0.30], 'LineWidth', 0.6);
    h_seg_eci   = gobjects(N_TRAIL_SEG, 1);
    for s = 1:N_TRAIL_SEG
        a = s / N_TRAIL_SEG;
        h_seg_eci(s) = plot(ax_eci, NaN, NaN, '-', ...
            'Color', a*c_A + (1-a)*BG_ECI, 'LineWidth', 0.8 + a*2.2);
    end
    h_dot_eci = plot(ax_eci, NaN, NaN, 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.85 0.85 0.90], 'LineWidth', 1.5);

    % ── Node diagram ─────────────────────────────────────────────────────────
    [h_nc, h_na] = local_draw_node_diagram(ax_node, D.famA, D.famBr, D.famB, c_A, c_Br, c_B);
    local_update_node_diagram(h_nc, h_na, 1, c_A, c_Br, c_B);  % init highlight

    % ── Video writer ──────────────────────────────────────────────────────────
    vid_path = fullfile(ex_dir, [base_name '_video.mp4']);
    vid = VideoWriter(vid_path, 'MPEG-4');
    vid.FrameRate = FPS;
    vid.Quality   = VIDEO_QUALITY;
    open(vid);

    fprintf('  Writing %d frames to:\n    %s\n', N_frames, vid_path);

    prev_phase = 0;

    for f = 1:N_frames
        idx = frame_idx(f);
        ph  = phase_id(idx);
        col = phase_col(idx, :);
        t_now = path_t(idx);

        % Trail window
        i_end   = idx;
        i_start = max(1, idx - trail_pts);
        seg_n   = max(1, floor((i_end - i_start) / N_TRAIL_SEG));

        % History (everything before the trail window)
        if i_start > 1
            set(h_hist_rot, 'XData', path_x(1:i_start),     'YData', path_y(1:i_start));
            set(h_hist_eci, 'XData', path_x_eci(1:i_start), 'YData', path_y_eci(1:i_start));
        end

        % Fading trail segments
        for s = 1:N_TRAIL_SEG
            is = i_start + (s-1)*seg_n;
            ie = min(i_start + s*seg_n, i_end);
            if is >= ie
                set(h_seg_rot(s), 'XData', NaN, 'YData', NaN);
                set(h_seg_eci(s), 'XData', NaN, 'YData', NaN);
                continue
            end
            a       = s / N_TRAIL_SEG;
            col_rot = a * col + (1-a) * BG_ROT;
            col_eci = a * col + (1-a) * BG_ECI;
            set(h_seg_rot(s), 'XData', path_x(is:ie),     'YData', path_y(is:ie),     'Color', col_rot);
            set(h_seg_eci(s), 'XData', path_x_eci(is:ie), 'YData', path_y_eci(is:ie), 'Color', col_eci);
        end

        % Spacecraft dot
        set(h_dot_rot, 'XData', path_x(idx),     'YData', path_y(idx));
        set(h_dot_eci, 'XData', path_x_eci(idx), 'YData', path_y_eci(idx));

        % Moon in ECI
        set(h_moon_eci, 'XData', cos(t_now), 'YData', sin(t_now));

        % Node highlight (only update on phase change)
        if ph ~= prev_phase
            local_update_node_diagram(h_nc, h_na, ph, c_A, c_Br, c_B);
            prev_phase = ph;
        end

        drawnow;
        writeVideo(vid, getframe(fig));

        if mod(f, round(N_frames/10)) == 0
            fprintf('    %3d%%\n', round(100*f/N_frames));
        end
    end

    close(vid);
    close(fig);
    fprintf('  Done: %s\n\n', vid_path);
end

fprintf('[video] All examples complete.\n');

% =============================================================================
% LOCAL HELPERS
% =============================================================================

function S = local_build_fam(fam_name, grid3)
    [mu, CJ, Tf_PO, X0] = rs3_core_family_ic(fam_name);
    S.name  = fam_name;
    S.mu    = mu;
    S.CJ    = CJ;
    S.Tf_PO = Tf_PO;
    S.X0    = X0;
    S.grid3 = grid3;
end

% ─────────────────────────────────────────────────────────────────────────────

function Xpo = local_integrate_po(S, ode_opts)
    sol = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t, X, S.CJ, S.mu, false), ...
                 [0, S.Tf_PO], S.X0, ode_opts);
    td  = linspace(0, min(S.Tf_PO, sol.x(end)), 1001)';
    Xpo = deval(sol, td)';
end

% ─────────────────────────────────────────────────────────────────────────────

function local_plot_po(ax, Xpo, rgb)
    xy = [Xpo(:,1:2); Xpo(1,1:2)];
    plot(ax, xy(:,1), xy(:,2), '--', 'Color', rgb, 'LineWidth', 1.1, ...
         'HandleVisibility', 'off');
end

% ─────────────────────────────────────────────────────────────────────────────

function [xl, yl] = local_sq(xl_in, yl_in)
% Return square (equal-aspect) x and y limits centred on the input ranges.
    cx   = mean(xl_in);  cy = mean(yl_in);
    half = max(diff(xl_in), diff(yl_in)) / 2;
    xl   = [cx - half, cx + half];
    yl   = [cy - half, cy + half];
end

% ─────────────────────────────────────────────────────────────────────────────

function s = local_short_name(fam_name)
    map = { ...
        'Lyapunov L1',            'Lyap. L1'; ...
        'Lyapunov L2',            'Lyap. L2'; ...
        'Cycler 21',              'Cycler 21'; ...
        'Cycler 11a',             sprintf('Cycler\n(1,1)a'); ...
        'Cycler 11b',             sprintf('Cycler\n(1,1)b'); ...
        'Cycler 32',              'Cycler 32'; ...
        'Resonant 2to1 Stable',   'R2:1 S'; ...
        'Resonant 2to1 Unstable', 'R2:1 U'; ...
        'Resonant 3to1 Stable',   'R3:1 S'; ...
        'Resonant 3to1 Unstable', 'R3:1 U'; ...
        'Resonant 5to2 Stable',   'R5:2 S'; ...
        'Resonant 5to2 Unstable', 'R5:2 U'; ...
        'Distant Prograde Orbit', 'DPO'; ...
    };
    idx = find(strcmp(map(:,1), fam_name), 1);
    if ~isempty(idx), s = map{idx,2}; else, s = fam_name; end
end

% ─────────────────────────────────────────────────────────────────────────────

function [h_nc, h_na] = local_draw_node_diagram(ax, famA, famBr, famB, c_A, c_Br, c_B)
% Draws the 3-node relay diagram.  Returns handles for dynamic highlighting.
    set(ax, 'XLim', [0 1], 'YLim', [0.25 0.75], 'Visible', 'off', ...
        'Color', 'none', 'XTick', [], 'YTick', []);
    hold(ax, 'on');

    cx = [0.15, 0.50, 0.85];
    cy = [0.50, 0.50, 0.50];
    r  = 0.13;
    th = linspace(0, 2*pi, 80);
    colors = {c_A, c_Br, c_B};
    names  = {local_short_name(famA), local_short_name(famBr), local_short_name(famB)};

    h_nc = gobjects(3, 1);
    for k = 1:3
        h_nc(k) = fill(ax, cx(k) + r*cos(th), cy(k) + r*sin(th), colors{k}, ...
            'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5);
        text(ax, cx(k), cy(k), names{k}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 7, 'FontWeight', 'bold', 'Color', 'w', 'Interpreter', 'none');
    end

    % Arrows (quiver with AutoScale off)
    gap = r + 0.01;
    h_na = gobjects(2, 1);
    h_na(1) = quiver(ax, cx(1)+gap, cy(1), cx(2)-cx(1)-2*gap, 0, 0, ...
        'Color', [0.2 0.2 0.2], 'LineWidth', 2.0, 'MaxHeadSize', 0.8, 'AutoScale', 'off');
    h_na(2) = quiver(ax, cx(2)+gap, cy(2), cx(3)-cx(2)-2*gap, 0, 0, ...
        'Color', [0.2 0.2 0.2], 'LineWidth', 2.0, 'MaxHeadSize', 0.8, 'AutoScale', 'off');

    % Thin border
    plot(ax, [0 1 1 0 0], [0.25 0.25 0.75 0.75 0.25], '-', ...
         'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
end

% ─────────────────────────────────────────────────────────────────────────────

function local_update_node_diagram(h_nc, h_na, phase, c_A, c_Br, c_B)
    dim = @(c) 0.35*c + 0.65*[1 1 1];
    COL_DIM_ARROW = [0.75 0.75 0.75];
    COL_ACT_ARROW = [0.15 0.15 0.15];

    switch phase
        case 1   % Leg 1 active: origin → bridge
            set(h_nc(1), 'FaceColor', c_A,       'LineWidth', 2.5);
            set(h_nc(2), 'FaceColor', c_Br,       'LineWidth', 2.5);
            set(h_nc(3), 'FaceColor', dim(c_B),   'LineWidth', 1.0);
            set(h_na(1), 'Color', COL_ACT_ARROW,  'LineWidth', 2.5);
            set(h_na(2), 'Color', COL_DIM_ARROW,  'LineWidth', 1.0);
        case 2   % Coast: on bridge
            set(h_nc(1), 'FaceColor', dim(c_A),   'LineWidth', 1.0);
            set(h_nc(2), 'FaceColor', c_Br,        'LineWidth', 3.0);
            set(h_nc(3), 'FaceColor', dim(c_B),   'LineWidth', 1.0);
            set(h_na(1), 'Color', COL_DIM_ARROW,  'LineWidth', 1.0);
            set(h_na(2), 'Color', COL_DIM_ARROW,  'LineWidth', 1.0);
        case 3   % Leg 2 active: bridge → dest
            set(h_nc(1), 'FaceColor', dim(c_A),   'LineWidth', 1.0);
            set(h_nc(2), 'FaceColor', c_Br,        'LineWidth', 2.5);
            set(h_nc(3), 'FaceColor', c_B,         'LineWidth', 2.5);
            set(h_na(1), 'Color', COL_DIM_ARROW,  'LineWidth', 1.0);
            set(h_na(2), 'Color', COL_ACT_ARROW,  'LineWidth', 2.5);
    end
end
