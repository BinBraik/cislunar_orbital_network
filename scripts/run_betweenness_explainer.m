%% RUN_BETWEENNESS_EXPLAINER
%
% Demonstrates betweenness centrality by showing that routing transfers
% through a high-betweenness "bridge" orbit (default: Cycler 11a) is
% cheaper than flying direct between origin and destination.
%
% Workflow:
%   1. Load the DV sweep matrix and pick a baseline (DV_cap, Tmax) snapshot.
%   2. Compute savings = DV_direct(A,B) - [DV(A,bridge) + DV(bridge,B)]
%      for every undirected pair that doesn't involve the bridge.
%   3. Select the top N_EXAMPLES pairs (ranked by savings).
%   4. For each selected pair:
%       a. Load atlases for origin, bridge, destination.
%       b. Extract Leg-1 arc  (origin  → bridge) via overlap_voxel_traj_extract.
%       c. Extract Leg-2 arc  (bridge  → dest  ) via overlap_voxel_traj_extract.
%       d. Compute coast segment on bridge PO between the two seeds.
%       e. Plot all elements (POs dashed + direction arrows, arcs in
%          distinct colours, patch markers, DV annotation).
%       f. [optional] Export an animated GIF of a spacecraft following
%          the full multi-leg path with a growing trail.
%
% Coordinate note:
%   All arcs are converted to physical (x, y) space for plotting:
%     FRS arcs   XA(:,1:2)          — already physical.
%     BRS arcs   reversed(x_B, y_B) — y_B already physical after R-transform;
%                                      reversed so arc runs patch→PO in forward time.
%   Periodic orbits are plotted from S.Xpo without modification.
%   A small geometric gap at each patch point is intentional; it
%   represents the DV_patch manoeuvre.
%
% Colour convention:
%   c_A  (blue)  — full Transfer A: origin departure arc + bridge arrival arc
%   c_Br (green) — bridge periodic orbit + coast segment on bridge
%   c_B  (red)   — full Transfer B: bridge departure arc + destination arrival arc

clear; clc;

% ── USER KNOBS ───────────────────────────────────────────────────────────────

% Input data (produced by run_overlap_dv_tmax_sweep / run_network_centrality_sweep)
SWEEP_MAT = fullfile(repo_root(), 'atlas_sweep_results', ...
                     'sweep_DVmatrix_results.mat');
NET_MAT   = fullfile(repo_root(), 'atlas_network_results', ...
                     'network_results.mat');

% Bridge (high-betweenness) orbit family
BRIDGE_FAMILY = 'Cycler 11a';

% Number of (origin, dest) examples to generate
N_EXAMPLES = 3;

% Baseline snapshot indices — [] uses the maximum-budget cell (last di, dj)
BASELINE_DI = [];
BASELINE_DJ = [];

% Animation
ANIM_ENABLE      = true;   % generate a GIF per example
ANIM_N_FRAMES    = 150;    % total frames in the GIF
ANIM_FRAME_DELAY = 0.04;   % seconds per frame

% Output directory
OUT_DIR = fullfile(repo_root(), 'betweenness_explainer');

% ── END OF USER KNOBS ────────────────────────────────────────────────────────

setup();
addpath(fullfile(repo_root(), 'src', 'network'));

if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

cfg = atlas_cfg_defaults();

cfg.grid.dx               = 0.001;
cfg.grid.dy               = 0.001;
cfg.grid.dtheta           = deg2rad(1);
cfg.seed.ds_seed          = 0.01;
cfg.propag.Tmax           = pi;
cfg.fan.DV_cap_nd         = 0.2;
cfg.fan.dtheta_fan        = deg2rad(0.5);
cfg.propag.absTol         = 1e-8;
cfg.propag.relTol         = 1e-8;
cfg.propag.v2tol          = 1e-8;
cfg.log.step_len_factor   = 0.75;
cfg.log.maxstep_factor    = 2;

cfg.io.save_figs        = true;
cfg.io.fig_visible      = 'on';
cfg.cache.enable        = true;
cfg.cache.rebuild       = false;
cfg.propag.Tmax         = pi;
cfg.fan.DV_cap_nd       = 0.2;
% Suppress intermediate bound plots (we only need B.imin)
cfg.plot.overlap.bounds_lb    = false;
cfg.plot.overlap.bounds_ub    = false;
cfg.plot.overlap.bounds_proxy = false;

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;

%% ── Load sweep data ─────────────────────────────────────────────────────────
fprintf('Loading sweep data:\n  %s\n', SWEEP_MAT);
if ~isfile(SWEEP_MAT)
    error('[betweenness_explainer] Sweep .mat not found:\n  %s\n', SWEEP_MAT);
end
S_sw = load(SWEEP_MAT, 'DV_cap_list', 'Tmax_list', 'DVmatrix_sweep', 'families');

DV_cap_list    = S_sw.DV_cap_list(:);
Tmax_list      = S_sw.Tmax_list(:);
DVmatrix_sweep = S_sw.DVmatrix_sweep;
families       = S_sw.families(:);

nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);
N     = numel(families);

%% ── Resolve baseline snapshot indices ───────────────────────────────────────
if isfile(NET_MAT) && (isempty(BASELINE_DI) || isempty(BASELINE_DJ))
    tmp = load(NET_MAT, 'BASELINE_DI', 'BASELINE_DJ');
    if isempty(BASELINE_DI) && isfield(tmp, 'BASELINE_DI')
        BASELINE_DI = tmp.BASELINE_DI;
    end
    if isempty(BASELINE_DJ) && isfield(tmp, 'BASELINE_DJ')
        BASELINE_DJ = tmp.BASELINE_DJ;
    end
end
if isempty(BASELINE_DI), BASELINE_DI = nDV;   end
if isempty(BASELINE_DJ), BASELINE_DJ = nTmax; end
BASELINE_DI = max(1, min(BASELINE_DI, nDV));
BASELINE_DJ = max(1, min(BASELINE_DJ, nTmax));

fprintf('Baseline snapshot: di=%d (DV_cap=%.3f nd), dj=%d (Tmax=%.3f nd)\n', ...
    BASELINE_DI, DV_cap_list(BASELINE_DI), BASELINE_DJ, Tmax_list(BASELINE_DJ));

%% ── Symmetrised DV matrix at baseline ───────────────────────────────────────
DV_raw = DVmatrix_sweep{BASELINE_DI, BASELINE_DJ};   % N×N, in m/s
DV_sym = min(DV_raw, DV_raw');                        % take cheaper direction
DV_sym(isnan(DV_sym)) = Inf;
DV_sym(1:N+1:end) = 0;                               % zero diagonal

%% ── Find bridge family index ─────────────────────────────────────────────────
bridge_idx = find(strcmpi(families, BRIDGE_FAMILY), 1);
if isempty(bridge_idx)
    error('[betweenness_explainer] Bridge family "%s" not found.\nAvailable: %s\n', ...
        BRIDGE_FAMILY, strjoin(families, ', '));
end
fprintf('Bridge family: "%s"  (index %d)\n', BRIDGE_FAMILY, bridge_idx);

%% ── Compute savings for all candidate pairs ──────────────────────────────────
savings_mat = NaN(N, N);
for iA = 1:N
    for iB = 1:N
        if iA == iB || iA == bridge_idx || iB == bridge_idx, continue; end
        dv_d  = DV_sym(iA, iB);
        dv_l1 = DV_sym(iA, bridge_idx);
        dv_l2 = DV_sym(bridge_idx, iB);
        if ~isinf(dv_d) && ~isinf(dv_l1) && ~isinf(dv_l2)
            savings_mat(iA, iB) = dv_d - (dv_l1 + dv_l2);
        end
    end
end

% Sort flat indices by savings descending, remove NaN/Inf
flat_savings = savings_mat(:);
[~, sorted_idx] = sort(flat_savings, 'descend', 'MissingPlacement', 'last');
sorted_idx = sorted_idx(isfinite(flat_savings(sorted_idx)));

% Pick top N_EXAMPLES (avoid selecting both A→B and B→A)
selected = zeros(0, 2);
for k = 1:numel(sorted_idx)
    if size(selected, 1) >= N_EXAMPLES, break; end
    [iA, iB] = ind2sub([N, N], sorted_idx(k));
    if ~any(selected(:,1)==iB & selected(:,2)==iA)
        selected(end+1, :) = [iA, iB];  %#ok<AGROW>
    end
end

if isempty(selected)
    error('[betweenness_explainer] No valid (origin, dest) pairs found via "%s".', BRIDGE_FAMILY);
end

fprintf('\nTop %d betweenness examples via "%s":\n', size(selected,1), BRIDGE_FAMILY);
for k = 1:size(selected, 1)
    iA = selected(k,1); iB = selected(k,2);
    fprintf('  %d: %-30s → %-30s | savings = %+.0f m/s  (direct~%.0f, via=%.0f+%.0f)\n', ...
        k, families{iA}, families{iB}, savings_mat(iA,iB), ...
        DV_sym(iA,iB), DV_sym(iA,bridge_idx), DV_sym(bridge_idx,iB));
end

%% ── Build shared grid ────────────────────────────────────────────────────────
grid3 = atlas_grid_make(cfg);

%% ── Per-example colour palette ───────────────────────────────────────────────
% [origin, bridge, dest] colours
PALETTE = [
    0.15 0.42 0.80;   % blue   — origin
    0.18 0.62 0.30;   % green  — bridge
    0.82 0.25 0.12;   % red    — destination
];

%% ── Process each example ────────────────────────────────────────────────────
for ex = 1:size(selected, 1)
    iA    = selected(ex, 1);
    iB    = selected(ex, 2);
    famA  = families{iA};
    famBr = BRIDGE_FAMILY;
    famB  = families{iB};

    fprintf('\n============================================================\n');
    fprintf('[ex %d]  %s  →  %s  →  %s\n', ex, famA, famBr, famB);
    fprintf('============================================================\n');

    ex_tag = sprintf('ex%d_%s_VIA_%s_TO_%s', ex, ...
        regexprep(famA,  '[^A-Za-z0-9]','_'), ...
        regexprep(famBr, '[^A-Za-z0-9]','_'), ...
        regexprep(famB,  '[^A-Za-z0-9]','_'));
    ex_dir = fullfile(OUT_DIR, ex_tag);
    if ~exist(ex_dir, 'dir'), mkdir(ex_dir); end

    % ── Load & densify atlases ───────────────────────────────────────────────
    fprintf('[ex %d] Loading atlases ...\n', ex);
    [SA,  ~] = atlas_prepare_or_load(famA,  cfg, grid3);
    [SBr, ~] = atlas_prepare_or_load(famBr, cfg, grid3);
    [SB,  ~] = atlas_prepare_or_load(famB,  cfg, grid3);

    relTol = cfg.propag.relTol;
    absTol = cfg.propag.absTol;
    SA  = local_ensure_xpo(SA,  relTol, absTol, 1001);
    SBr = local_ensure_xpo(SBr, relTol, absTol, 1001);
    SB  = local_ensure_xpo(SB,  relTol, absTol, 1001);

    % ── Leg 1: origin → bridge ───────────────────────────────────────────────
    fprintf('[ex %d] Leg 1: %s → %s\n', ex, famA, famBr);
    O1 = overlap_pair(SA, SBr, cfg);
    if isempty(O1) || ~isfield(O1,'ids') || isempty(O1.ids)
        fprintf('[ex %d]  WARNING: no overlap for Leg 1 — skipping example.\n', ex);
        continue
    end
    V1 = overlap_extract_voxel_info(SA, SBr, O1, cfg);
    B1 = overlap_visualize_bounds(V1, SA, SBr, O1, cfg, ex_dir, [ex_tag '_leg1']);
    T1 = overlap_voxel_traj_extract(SA, SBr, V1, B1, cfg);

    % ── Leg 2: bridge → dest ─────────────────────────────────────────────────
    fprintf('[ex %d] Leg 2: %s → %s\n', ex, famBr, famB);
    O2 = overlap_pair(SBr, SB, cfg);
    if isempty(O2) || ~isfield(O2,'ids') || isempty(O2.ids)
        fprintf('[ex %d]  WARNING: no overlap for Leg 2 — skipping example.\n', ex);
        continue
    end
    V2 = overlap_extract_voxel_info(SBr, SB, O2, cfg);
    B2 = overlap_visualize_bounds(V2, SBr, SB, O2, cfg, ex_dir, [ex_tag '_leg2']);
    T2 = overlap_voxel_traj_extract(SBr, SB, V2, B2, cfg);

    % ── DV summary ───────────────────────────────────────────────────────────
    dv_l1      = T1.DV_total_true_mps;
    dv_l2      = T2.DV_total_true_mps;
    dv_via     = dv_l1 + dv_l2;
    dv_direct  = DV_sym(iA, iB);   % proxy from sweep (upper bound for direct)
    savings    = dv_direct - dv_via;

    fprintf('[ex %d] Leg1=%.1f  +  Leg2=%.1f  =  %.1f m/s\n', ex, dv_l1, dv_l2, dv_via);
    fprintf('[ex %d] Direct proxy=%.1f m/s  |  Savings=%.1f m/s (%.0f%%)\n', ...
        ex, dv_direct, savings, 100*savings/dv_direct);

    % ── Coast arc on bridge PO ───────────────────────────────────────────────
    % Leg-1 bridge arrival seed (physical):  T1.seed_B_frs(1:2)
    % Leg-2 bridge departure seed (physical): T2.seed_A(1:2)
    % Bridge arrival in physical (x,y):  R-transform of FRS IC = (x_B(1), y_B(1))
    % Bridge departure in physical (x,y): FRS IC of Leg-2 = T2.XA(1,1:2)
    [coast_arc, coast_time] = local_coast_on_po(SBr, ...
        [T1.x_B(1); T1.y_B(1)], ...
        [T2.XA(1,1); T2.XA(1,2)]);

    % ── Static figure ────────────────────────────────────────────────────────
    c_A  = PALETTE(1,:);
    c_Br = PALETTE(2,:);
    c_B  = PALETTE(3,:);

    CJbg = min([SA.CJ, SBr.CJ, SB.CJ]);
    mu   = SA.mu;

    fig = figure('Color','w', ...
        'Name',    sprintf('Betweenness Example %d: %s via %s', ex, famA, famBr), ...
        'Visible', cfg.io.fig_visible, ...
        'Units',  'pixels', 'Position', [80 80 1050 820]);
    ax = axes('Parent', fig);
    cr3bp_plot_background(CJbg, mu, ax);
    set(ax.Children, 'HandleVisibility', 'off');
    hold(ax, 'on');
    axis(ax, 'equal');
    grid(ax, 'on');

    % ── Periodic orbits (dashed) with direction arrows ────────────────────────
    local_plot_po_with_arrows(ax, SA,  c_A,  sprintf('Origin: %s',  famA));
    local_plot_po_with_arrows(ax, SBr, c_Br, sprintf('Bridge: %s',  famBr));
    local_plot_po_with_arrows(ax, SB,  c_B,  sprintf('Dest: %s',    famB));

    % ── Transfer A: both arcs in c_A ─────────────────────────────────────────
    % FRS: origin → patch-1  (physical, forward in time)
    plot(ax, T1.XA(:,1), T1.XA(:,2), '-', 'Color', c_A, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Transfer A (%.0f + %.0f d)', T1.tof_A_days, T1.tof_B_days));
    % BRS: patch-1 → bridge  (reversed; y_B already physical after R-transform)
    plot(ax, T1.x_B(end:-1:1), T1.y_B(end:-1:1), '-', 'Color', c_A, 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');

    % ── Coast arc on bridge ───────────────────────────────────────────────────
    if ~isempty(coast_arc) && size(coast_arc,1) > 1
        plot(ax, coast_arc(:,1), coast_arc(:,2), '-', 'Color', c_Br, ...
            'LineWidth', 3.5, 'DisplayName', sprintf('Bridge coast (%.0f d)', coast_time*TU_days));
    end

    % ── Transfer B: both arcs in c_B ─────────────────────────────────────────
    % FRS: bridge → patch-2  (physical, forward in time)
    plot(ax, T2.XA(:,1), T2.XA(:,2), '-', 'Color', c_B, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Transfer B (%.0f + %.0f d)', T2.tof_A_days, T2.tof_B_days));
    % BRS: patch-2 → dest  (reversed; y_B already physical after R-transform)
    plot(ax, T2.x_B(end:-1:1), T2.y_B(end:-1:1), '-', 'Color', c_B, 'LineWidth', 2.0, ...
        'HandleVisibility', 'off');

    % ── ΔV patch markers (yellow stars) ──────────────────────────────────────
    xp1 = T1.XA(T1.i_star, 1);  yp1 = T1.XA(T1.i_star, 2);
    xp2 = T2.XA(T2.i_star, 1);  yp2 = T2.XA(T2.i_star, 2);
    plot(ax, xp1, yp1, 'p', 'Color','k', 'MarkerFaceColor',[0.95 0.85 0.05], ...
        'MarkerSize', 14, 'LineWidth', 1.2, 'DisplayName', '\DeltaV patch');
    plot(ax, xp2, yp2, 'p', 'Color','k', 'MarkerFaceColor',[0.95 0.85 0.05], ...
        'MarkerSize', 14, 'LineWidth', 1.2, 'HandleVisibility','off');

    % ── Event markers: depart A, bridge arrival/departure, arrive B ───────────
    % Departure from A
    plot(ax, T1.XA(1,1),  T1.XA(1,2),  'o', 'Color','k', ...
        'MarkerFaceColor', c_A, 'MarkerSize', 11, 'LineWidth', 1.5, ...
        'DisplayName', 'Depart A');
    % Bridge arrival (end of Transfer A inbound arc)
    plot(ax, T1.x_B(1),   T1.y_B(1),   'd', 'Color','k', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 11, 'LineWidth', 2.0, ...
        'DisplayName', 'Bridge arrival');
    % Bridge departure (start of Transfer B outbound arc)
    plot(ax, T2.XA(1,1),  T2.XA(1,2),  'd', 'Color','k', ...
        'MarkerFaceColor', c_Br, 'MarkerSize', 11, 'LineWidth', 1.5, ...
        'DisplayName', 'Bridge departure');
    % Arrival at B
    plot(ax, T2.x_B(1),   T2.y_B(1),   'o', 'Color','k', ...
        'MarkerFaceColor', c_B, 'MarkerSize', 11, 'LineWidth', 1.5, ...
        'DisplayName', 'Arrive B');

    % ── DV annotation ─────────────────────────────────────────────────────────
    title(ax, { ...
        sprintf('%s  \\rightarrow  %s  \\rightarrow  %s', famA, famBr, famB), ...
        sprintf('Via bridge: %.0f + %.0f = \\bf%.0f m/s\\rm   |   Direct (proxy): %.0f m/s   |   Savings: \\bf%.0f m/s (%.0f%%)', ...
            dv_l1, dv_l2, dv_via, dv_direct, savings, 100*savings/dv_direct)}, ...
        'Interpreter','tex', 'FontSize', 10);

    xlabel(ax, 'x [nd]');
    ylabel(ax, 'y [nd]');
    legend(ax, 'Location','best', 'FontSize', 8);

    io_save_figure(fig, ex_dir, ['betweenness_' ex_tag], cfg);

    % ── Optional GIF animation ────────────────────────────────────────────────
    if ANIM_ENABLE
        local_make_gif(T1, T2, coast_arc, coast_time, SA, SBr, SB, c_A, c_Br, c_B, ...
            CJbg, mu, famA, famBr, famB, dv_l1, dv_l2, dv_via, dv_direct, ...
            ANIM_N_FRAMES, ANIM_FRAME_DELAY, ex_dir, ex_tag, cfg);
    end

    fprintf('[ex %d] Done. Output: %s\n', ex, ex_dir);
end

fprintf('\n[betweenness_explainer] All examples complete. Output root: %s\n', OUT_DIR);

% =============================================================================
% LOCAL HELPERS
% =============================================================================

function local_plot_po_with_arrows(ax, S, rgb, dispName)
%LOCAL_PLOT_PO_WITH_ARROWS  Dashed PO + 4 equally-spaced direction arrows.
if isfield(S,'Xpo') && ~isempty(S.Xpo)
    xy = S.Xpo(:, 1:2);
elseif isfield(S,'PO_xy') && ~isempty(S.PO_xy)
    xy = S.PO_xy;
else
    return
end
if norm(xy(end,:) - xy(1,:)) > 1e-6
    xy = [xy; xy(1,:)];
end
n = size(xy, 1);

plot(ax, xy(:,1), xy(:,2), '--', 'Color', rgb, 'LineWidth', 1.4, ...
    'DisplayName', dispName);

% Direction arrows (quiver)  — 4 equally-spaced along the orbit
N_arr = 4;
arr_len = 0.025;   % arrow length [nd]
for k = 1:N_arr
    i1 = max(1, round(k * n / N_arr));
    i2 = min(n, i1 + max(1, round(n/30)));
    dx = xy(i2,1) - xy(i1,1);
    dy = xy(i2,2) - xy(i1,2);
    nrm = hypot(dx, dy);
    if nrm < 1e-8, continue; end
    quiver(ax, xy(i1,1), xy(i1,2), arr_len*dx/nrm, arr_len*dy/nrm, 0, ...
        'Color', rgb, 'LineWidth', 1.0, 'MaxHeadSize', 3, ...
        'HandleVisibility','off');
end
end

% ─────────────────────────────────────────────────────────────────────────────

function [coast_arc, coast_time] = local_coast_on_po(SBr, p_arrive, p_depart)
%LOCAL_COAST_ON_PO  Segment of bridge PO from arrival to departure (physical).
% p_arrive [2×1]: physical (x,y) of bridge arrival point  (end of Transfer A)
% p_depart [2×1]: physical (x,y) of bridge departure point (start of Transfer B)
%
% Xpo is stored uniformly in time over [0, Tf_PO].  The coast time is
% computed from the fraction of the orbit traversed.
%
% Returns:
%   coast_arc  [n×2] (x,y) sampled from Xpo at the Xpo resolution
%   coast_time  time duration of the coast in non-dimensional units

coast_arc  = zeros(0, 2);
coast_time = 0;
if ~isfield(SBr,'Xpo') || isempty(SBr.Xpo), return; end

xy  = SBr.Xpo(:, 1:2);
N_po = size(xy, 1);
dt_po = SBr.Tf_PO / max(1, N_po - 1);   % time step between Xpo samples

[~, i_arr] = min(hypot(xy(:,1)-p_arrive(1), xy(:,2)-p_arrive(2)));
[~, i_dep] = min(hypot(xy(:,1)-p_depart(1), xy(:,2)-p_depart(2)));

if i_dep >= i_arr
    coast_arc = xy(i_arr:i_dep, :);
    n_steps   = i_dep - i_arr;
else
    % Wrap: continue to end of array, then from start to i_dep
    coast_arc = [xy(i_arr:end, :); xy(1:i_dep, :)];
    n_steps   = (N_po - i_arr) + i_dep;
end
coast_time = n_steps * dt_po;
end

% ─────────────────────────────────────────────────────────────────────────────

function local_make_gif(T1, T2, coast_arc, coast_time, SA, SBr, SB, c_A, c_Br, c_B, ...
    CJbg, mu, famA, famBr, famB, dv_l1, dv_l2, dv_via, dv_direct, ...
    N_frames, frame_delay, out_dir, ex_tag, cfg)
%LOCAL_MAKE_GIF  Animated GIF: spacecraft travels the full multi-leg path.
%
% Path phases (physical coords, forward in time):
%   1. Transfer A out : T1.XA               (origin  → patch-1)   colour: c_A
%   2. Transfer A in  : T1.BRS reversed     (patch-1 → bridge)    colour: c_A
%   3. Coast          : coast_arc           (bridge arrival → departure) colour: c_Br
%   4. Transfer B out : T2.XA               (bridge  → patch-2)   colour: c_B
%   5. Transfer B in  : T2.BRS reversed     (patch-2 → dest)      colour: c_B
%
% Frames are allocated proportional to each phase's physical time duration
% so the animated dot moves at a uniform rate across all segments.

TU_days  = local_cfg_get(cfg, 'units.TU_days', 1.0);
relTol   = local_cfg_get(cfg, 'propag.relTol', 1e-9);
absTol   = local_cfg_get(cfg, 'propag.absTol', 1e-9);
ode_opts = odeset('RelTol', relTol, 'AbsTol', absTol);

% ── Step 1: Dense integration via deval — exact trajectories, no shortcuts ────
% N_dense >> N_frames so subsequent arc-length resampling is smooth.
N_dense = max(2000, 20 * N_frames);

sol1 = ode113(@(t,X) cr3bp_reduced_ode(t,X,SA.CJ, SA.mu, false), ...
               [0, T1.t_A], T1.IC_A, ode_opts);
D1   = deval(sol1, linspace(0, T1.t_A, N_dense))';   % N_dense×3, physical

sol2 = ode113(@(t,X) cr3bp_reduced_ode(t,X,SBr.CJ,SBr.mu,false), ...
               [0, T1.t_B], T1.IC_B_frs, ode_opts);
D2r  = deval(sol2, linspace(0, T1.t_B, N_dense))';   % FRS frame
D2   = [D2r(end:-1:1,1), -D2r(end:-1:1,2)];          % reversed + R-transform → physical

sol4 = ode113(@(t,X) cr3bp_reduced_ode(t,X,SBr.CJ,SBr.mu,false), ...
               [0, T2.t_A], T2.IC_A, ode_opts);
D4   = deval(sol4, linspace(0, T2.t_A, N_dense))';

sol5 = ode113(@(t,X) cr3bp_reduced_ode(t,X,SB.CJ, SB.mu, false), ...
               [0, T2.t_B], T2.IC_B_frs, ode_opts);
D5r  = deval(sol5, linspace(0, T2.t_B, N_dense))';
D5   = [D5r(end:-1:1,1), -D5r(end:-1:1,2)];

% Coast: Xpo is already dense and uniform in time
D3 = coast_arc;   % may be empty

% ── Step 2: Allocate frames proportional to arc length ────────────────────────
% This gives uniform *visual* speed regardless of physical speed variations.
arc_len = @(xy) sum(hypot(diff(xy(:,1)), diff(xy(:,2))));
L3 = 0;
if size(D3,1) > 1, L3 = arc_len(D3); end
L = max([arc_len(D1(:,1:2)), arc_len(D2), L3, arc_len(D4(:,1:2)), arc_len(D5)], 1e-12);
L_total = sum(L);

n_ph = max(2, round(N_frames * L / L_total));
n_ph(end) = max(2, N_frames - sum(n_ph(1:end-1)));

% ── Step 3: Resample each phase to uniform arc-length steps ──────────────────
[ph1_x, ph1_y] = local_arc_resample(D1(:,1),   D1(:,2),   n_ph(1));
[ph2_x, ph2_y] = local_arc_resample(D2(:,1),   D2(:,2),   n_ph(2));
if size(D3,1) > 1
    [ph3_x, ph3_y] = local_arc_resample(D3(:,1), D3(:,2), n_ph(3));
else
    ph3_x = zeros(0,1); ph3_y = zeros(0,1);
end
[ph4_x, ph4_y] = local_arc_resample(D4(:,1),   D4(:,2),   n_ph(4));
[ph5_x, ph5_y] = local_arc_resample(D5(:,1),   D5(:,2),   n_ph(5));

path_x = [ph1_x; ph2_x; ph3_x; ph4_x; ph5_x];
path_y = [ph1_y; ph2_y; ph3_y; ph4_y; ph5_y];

% Phase colour at each resampled point
n_c_actual = numel(ph3_x);
phase_col = [repmat(c_A,  n_ph(1),    1); ...
             repmat(c_A,  n_ph(2),    1); ...
             repmat(c_Br, n_c_actual, 1); ...
             repmat(c_B,  n_ph(4),    1); ...
             repmat(c_B,  n_ph(5),    1)];

n_total = numel(path_x);
if n_total < 2
    warning('[betweenness_gif] Path too short — skipping animation.');
    return
end

% ── Create figure ─────────────────────────────────────────────────────────────
fig_anim = figure('Color','w', 'Visible','off', ...
                  'Units','pixels', 'Position',[100 100 950 760]);
ax = axes('Parent', fig_anim);
cr3bp_plot_background(CJbg, mu, ax);
set(ax.Children, 'HandleVisibility','off');
hold(ax,'on');
axis(ax,'equal');
grid(ax,'on');

local_plot_po_with_arrows(ax, SA,  c_A,  '');
local_plot_po_with_arrows(ax, SBr, c_Br, '');
local_plot_po_with_arrows(ax, SB,  c_B,  '');

% ── Static event markers (plotted once) ──────────────────────────────────────
plot(ax, T1.XA(1,1),  T1.XA(1,2),  'o', 'Color','k', ...
    'MarkerFaceColor', c_A,  'MarkerSize',10, 'LineWidth',1.5, 'HandleVisibility','off');
plot(ax, T1.XA(T1.i_star,1), T1.XA(T1.i_star,2), 'p', 'Color','k', ...
    'MarkerFaceColor',[0.95 0.85 0.05], 'MarkerSize',13, 'LineWidth',1.2, 'HandleVisibility','off');
plot(ax, T1.x_B(1),   T1.y_B(1),   'd', 'Color','k', ...
    'MarkerFaceColor', 'w',   'MarkerSize',10, 'LineWidth',2.0, 'HandleVisibility','off');
plot(ax, T2.XA(1,1),  T2.XA(1,2),  'd', 'Color','k', ...
    'MarkerFaceColor', c_Br, 'MarkerSize',10, 'LineWidth',1.5, 'HandleVisibility','off');
plot(ax, T2.XA(T2.i_star,1), T2.XA(T2.i_star,2), 'p', 'Color','k', ...
    'MarkerFaceColor',[0.95 0.85 0.05], 'MarkerSize',13, 'LineWidth',1.2, 'HandleVisibility','off');
plot(ax, T2.x_B(1),   T2.y_B(1),   'o', 'Color','k', ...
    'MarkerFaceColor', c_B,  'MarkerSize',10, 'LineWidth',1.5, 'HandleVisibility','off');

t_coast_days = coast_time * TU_days;
title(ax, { ...
    sprintf('%s  \\rightarrow  %s  \\rightarrow  %s', famA, famBr, famB), ...
    sprintf('Transfer A: %.0f m/s  |  Coast: %.0f d  |  Transfer B: %.0f m/s  |  Total: %.0f m/s  |  Direct: %.0f m/s', ...
        dv_l1, t_coast_days, dv_l2, dv_via, dv_direct)}, ...
    'Interpreter','tex', 'FontSize', 8);
xlabel(ax,'x [nd]'); ylabel(ax,'y [nd]');

% Dot + growing trail (updated each frame)
h_trail = plot(ax, NaN, NaN, '-', 'Color', [0.3 0.3 0.3], ...
    'LineWidth', 1.4, 'HandleVisibility','off');
h_dot   = plot(ax, NaN, NaN, 'o', 'Color','k', ...
    'MarkerFaceColor', c_A, 'MarkerSize', 10, 'HandleVisibility','off');

gif_path = fullfile(out_dir, ['betweenness_' ex_tag '_animation.gif']);

for fi = 1:n_total
    set(h_trail, 'XData', path_x(1:fi), 'YData', path_y(1:fi));
    set(h_dot,   'XData', path_x(fi),   'YData', path_y(fi), ...
                 'MarkerFaceColor', phase_col(fi,:));
    drawnow;
    frame = getframe(fig_anim);
    [im_idx, cmap] = rgb2ind(frame2im(frame), 256);
    if fi == 1
        imwrite(im_idx, cmap, gif_path, 'gif', 'LoopCount', Inf, 'DelayTime', frame_delay);
    else
        imwrite(im_idx, cmap, gif_path, 'gif', 'WriteMode','append', 'DelayTime', frame_delay);
    end
end

close(fig_anim);
fprintf('[anim] GIF saved: %s\n', gif_path);
end

% ─────────────────────────────────────────────────────────────────────────────

function [xr, yr] = local_arc_resample(x, y, n_out)
%LOCAL_ARC_RESAMPLE  Resample (x,y) curve to n_out points uniform in arc length.
% Gives visually uniform dot speed regardless of physical velocity variations.
x = x(:); y = y(:);
ds = hypot(diff(x), diff(y));
s  = [0; cumsum(ds)];
% Clamp any floating-point duplicates to keep interp1 happy
s  = s + (0:numel(s)-1)' * 1e-14 * max(s(end), 1);
s_q = linspace(0, s(end), n_out);
xr  = interp1(s, x, s_q, 'pchip')';
yr  = interp1(s, y, s_q, 'pchip')';
end

% ─────────────────────────────────────────────────────────────────────────────

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

% ─────────────────────────────────────────────────────────────────────────────

function S = local_ensure_xpo(S, relTol, absTol, N_po)
%LOCAL_ENSURE_XPO  Re-integrate PO if dense Xpo was stripped from cache.
if isfield(S,'Xpo') && ~isempty(S.Xpo) && ...
   isfield(S,'t_dense') && ~isempty(S.t_dense)
    return
end
fprintf('[betweenness_explainer] Re-integrating PO for "%s" ...\n', S.name);
opts  = odeset('RelTol', relTol, 'AbsTol', absTol);
sol   = ode113(@(t,X) cr3bp_reduced_ode(t,X,S.CJ,S.mu,true), ...
               [0, S.Tf_PO], S.X0, opts);
td    = linspace(0, S.Tf_PO, N_po)';
S.t_dense = td;
S.Xpo     = deval(sol, td)';
end
