%% =========================================================================
%  run_network_centrality_sweep.m
%
%  Network centrality analysis over the DV-cap / Tmax parameter sweep.
%
%  For each snapshot (di, dj) in the sweep grid the script:
%    1. Loads and sanitises the DV/TOF transfer matrix
%    2. Applies the feasibility filter
%    3. Runs Floyd-Warshall for all-pairs shortest-path distances
%    4. Computes the Largest Connected Component (LCC)
%    5-9. Computes budgeted reachability, reach, gateway, hub, and
%         articulation-point metrics (see ALGORITHM in task description)
%    10. Resolves ties and records winners
%
%  Outputs (all under OUT_DIR):
%    snapshot_summary.csv, node_metrics.csv, network_results.mat
%    winner_map_reach.{png,fig}, winner_map_gateway.{png,fig}
%    winner_map_hub.{png,fig}, lcc_map.{png,fig}, articulation_map.{png,fig}
%    [optional] animation_Tmax<dj>.gif
%
%  Dependencies: rs3_cfg_defaults, src/network/*.m  (added to path below)
% =========================================================================

%% ========================================================================
%  USER KNOBS — edit this block; do not change code below the divider
% =========================================================================

% Path to the sweep .mat produced by run_rs4_dv_tmax_sweep.m
SWEEP_MAT = fullfile(rs3_repo_root(), 'rs3_sweep_results', ...
                     'sweep_DVmatrix_results.mat');

% Output directory (created automatically if absent)
OUT_DIR = fullfile(rs3_repo_root(), 'rs3_network_results');

% Physical budget multiplier  (2 = departure manoeuvre + arrival manoeuvre)
BUDGET_FACTOR = 2;

% Reach guard: only nodes with R_budget >= RMIN_GUARD qualify as reach
% winner.  Set to 0 to disable (all nodes always qualify).
RMIN_GUARD = 0.75;

% Tie tolerances
TIE_TOL_REL         = 1e-4;    % relative tolerance
TIE_TOL_ABS_SCALE   = 1e-2;   % absolute floor = TIE_TOL_ABS_SCALE * VU_mps

% GIF animation settings
GIF_ENABLE      = false;   % set true to generate animated GIF
GIF_TMAX_IDX    = 10;      % fixed dj (Tmax column index) for the GIF
GIF_FRAME_DELAY = 0.5;     % seconds per frame

% =========================================================================
%  END OF USER KNOBS
% =========================================================================

%% Setup ------------------------------------------------------------------
% Add repo src/ and src/network/ to path
rs3_setup();
addpath(fullfile(rs3_repo_root(), 'src', 'network'));

if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

% Load unit constants — never hardcode these
cfg     = rs3_cfg_defaults();
VU_mps  = cfg.units.VU_mps;    % ~1023.2 m/s
TU_days = cfg.units.TU_days;   % ~4.348 days

TIE_TOL_ABS = TIE_TOL_ABS_SCALE * VU_mps;   % ~10 m/s floor

%% Load sweep data --------------------------------------------------------
fprintf('Loading sweep data:\n  %s\n', SWEEP_MAT);
if ~exist(SWEEP_MAT, 'file')
    error('run_network_centrality_sweep: sweep .mat not found:\n  %s\n', SWEEP_MAT);
end
S = load(SWEEP_MAT, 'DV_cap_list', 'Tmax_list', 'DVmatrix_sweep', ...
         'TOFmatrix_sweep', 'families');

DV_cap_list     = S.DV_cap_list(:);   % [nDV×1]
Tmax_list       = S.Tmax_list(:);     % [nTmax×1]
DVmatrix_sweep  = S.DVmatrix_sweep;   % {nDV×nTmax}
TOFmatrix_sweep = S.TOFmatrix_sweep;  % {nDV×nTmax}
full_names      = S.families(:);      % {N×1} full family names

nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);
N     = 13;                            % fixed number of families

if numel(full_names) ~= N
    warning('Expected %d families in sweep .mat; found %d.', N, numel(full_names));
    N = numel(full_names);
end

% Physical axis vectors (used as figure axes and stored in outputs)
DV_vec   = BUDGET_FACTOR * DV_cap_list * VU_mps;    % [nDV×1]  m/s
Tmax_vec = BUDGET_FACTOR * Tmax_list   * TU_days;   % [nTmax×1] days

short_names = net_family_short_names();

fprintf('  Grid: %d DV-cap values × %d Tmax values = %d snapshots\n', ...
        nDV, nTmax, nDV*nTmax);

%% Preallocate result arrays ----------------------------------------------

% Snapshot-level maps [nDV×nTmax]
skip_map       = true(nDV, nTmax);
edges_kept_map = zeros(nDV, nTmax);
lcc_size_map   = zeros(nDV, nTmax);
lcc_full_map   = zeros(nDV, nTmax);
ap_count_map   = zeros(nDV, nTmax);

% Per-node metrics [N × nDV × nTmax]
R_budget_all = NaN(N, nDV, nTmax);
reach_all    = NaN(N, nDV, nTmax);
gateway_all  = NaN(N, nDV, nTmax);
hub_all      = NaN(N, nDV, nTmax);
is_ap_all    = false(N, nDV, nTmax);

% Winner maps [nDV×nTmax]:  k = family index, 0 = Tie, NaN = skip
reach_winner_idx   = NaN(nDV, nTmax);
gateway_winner_idx = NaN(nDV, nTmax);
hub_winner_idx     = NaN(nDV, nTmax);

reach_tie_sz       = zeros(nDV, nTmax);
gateway_tie_sz     = zeros(nDV, nTmax);
hub_tie_sz         = zeros(nDV, nTmax);

% Winner name strings (short, semicolon-joined if tie)
reach_winner_names   = repmat({''}, nDV, nTmax);
gateway_winner_names = repmat({''}, nDV, nTmax);
hub_winner_names     = repmat({''}, nDV, nTmax);

%% Main sweep loop --------------------------------------------------------
fprintf('Processing snapshots...\n');
t_start = tic;

for di = 1:nDV
    for dj = 1:nTmax

        % ── Steps 1-2: build graph ────────────────────────────────────────
        [A, W, D_sym, ~, edges_kept, DVcap_true, Tmax_true, skip] = ...
            net_build_graph( ...
                DVmatrix_sweep{di, dj}, TOFmatrix_sweep{di, dj}, ...
                DV_cap_list(di), Tmax_list(dj), ...
                VU_mps, TU_days, BUDGET_FACTOR);

        if skip, continue; end

        skip_map(di, dj)       = false;
        edges_kept_map(di, dj) = edges_kept;

        % ── Step 3: all-pairs shortest paths ─────────────────────────────
        dist = net_floyd_warshall(W);

        % ── Step 4: LCC ──────────────────────────────────────────────────
        [lcc_sz, lcc_full, ~] = net_lcc(A);
        lcc_size_map(di, dj)  = lcc_sz;
        lcc_full_map(di, dj)  = lcc_full;

        % ── Steps 5-9: centrality metrics ────────────────────────────────
        metrics = net_centrality(A, W, D_sym, dist, DVcap_true, VU_mps, RMIN_GUARD);

        R_budget_all(:, di, dj) = metrics.R_budget;
        reach_all(:,    di, dj) = metrics.reach;
        gateway_all(:,  di, dj) = metrics.gateway;
        hub_all(:,      di, dj) = metrics.hub;
        is_ap_all(:,    di, dj) = metrics.is_articulation;
        ap_count_map(di, dj)    = sum(metrics.is_articulation);

        % ── Step 10: tie handling per metric ─────────────────────────────
        [reach_winner_idx(di,dj), reach_tie_sz(di,dj), reach_winner_names{di,dj}] = ...
            local_get_winner(metrics.reach, short_names, N, ...
                             TIE_TOL_REL, TIE_TOL_ABS, ...
                             RMIN_GUARD, metrics.R_budget);

        [gateway_winner_idx(di,dj), gateway_tie_sz(di,dj), gateway_winner_names{di,dj}] = ...
            local_get_winner(metrics.gateway, short_names, N, ...
                             TIE_TOL_REL, TIE_TOL_ABS, 0, []);

        [hub_winner_idx(di,dj), hub_tie_sz(di,dj), hub_winner_names{di,dj}] = ...
            local_get_winner(metrics.hub, short_names, N, ...
                             TIE_TOL_REL, TIE_TOL_ABS, 0, []);
    end

    if mod(di, max(1, floor(nDV/5))) == 0
        fprintf('  di=%d/%d done (%.1f s elapsed)\n', di, nDV, toc(t_start));
    end
end

fprintf('Sweep complete in %.1f s.\n', toc(t_start));

%% Save snapshot_summary.csv ----------------------------------------------
csv1_path = fullfile(OUT_DIR, 'snapshot_summary.csv');
fprintf('Writing %s ...\n', csv1_path);

fid = fopen(csv1_path, 'w');
fprintf(fid, ['di,dj,DVcap_nd,Tmax_nd,DVcap_true_mps,Tmax_true_days,' ...
              'edges_kept,lcc_size,lcc_full,' ...
              'reach_winner,reach_tie_size,' ...
              'gateway_winner,gateway_tie_size,' ...
              'hub_winner,hub_tie_size,' ...
              'n_articulation_points\n']);

for di = 1:nDV
    for dj = 1:nTmax
        if skip_map(di, dj)
            % Write a skip row with NaN for all computed fields
            fprintf(fid, '%d,%d,%.6g,%.6g,%.4f,%.4f,%d,%d,%d,%s,%d,%s,%d,%s,%d,%d\n', ...
                di, dj, ...
                DV_cap_list(di), Tmax_list(dj), ...
                DV_vec(di), Tmax_vec(dj), ...
                0, 0, 0, ...
                'skip', 0, ...
                'skip', 0, ...
                'skip', 0, ...
                0);
            continue
        end

        fprintf(fid, '%d,%d,%.6g,%.6g,%.4f,%.4f,%d,%d,%d,%s,%d,%s,%d,%s,%d,%d\n', ...
            di, dj, ...
            DV_cap_list(di), Tmax_list(dj), ...
            DV_vec(di), Tmax_vec(dj), ...
            edges_kept_map(di,dj), lcc_size_map(di,dj), lcc_full_map(di,dj), ...
            local_winner_str(reach_winner_names{di,dj}, reach_winner_idx(di,dj)), ...
            reach_tie_sz(di,dj), ...
            local_winner_str(gateway_winner_names{di,dj}, gateway_winner_idx(di,dj)), ...
            gateway_tie_sz(di,dj), ...
            local_winner_str(hub_winner_names{di,dj}, hub_winner_idx(di,dj)), ...
            hub_tie_sz(di,dj), ...
            ap_count_map(di,dj));
    end
end

fclose(fid);

%% Save node_metrics.csv --------------------------------------------------
csv2_path = fullfile(OUT_DIR, 'node_metrics.csv');
fprintf('Writing %s ...\n', csv2_path);

fid2 = fopen(csv2_path, 'w');
fprintf(fid2, ['di,dj,DVcap_true_mps,Tmax_true_days,' ...
               'node_idx,node_name,' ...
               'R_budget,reach,gateway,hub,is_articulation\n']);

for di = 1:nDV
    for dj = 1:nTmax
        if skip_map(di, dj), continue; end

        for k = 1:N
            fprintf(fid2, '%d,%d,%.4f,%.4f,%d,"%s",%.6g,%.6g,%.6g,%.6g,%d\n', ...
                di, dj, ...
                DV_vec(di), Tmax_vec(dj), ...
                k, full_names{k}, ...
                R_budget_all(k, di, dj), ...
                reach_all(k,    di, dj), ...
                gateway_all(k,  di, dj), ...
                hub_all(k,      di, dj), ...
                int8(is_ap_all(k, di, dj)));
        end
    end
end

fclose(fid2);

%% Save .mat bundle -------------------------------------------------------
mat_path = fullfile(OUT_DIR, 'network_results.mat');
fprintf('Writing %s ...\n', mat_path);

results = struct();
results.DV_cap_list      = DV_cap_list;
results.Tmax_list        = Tmax_list;
results.DV_vec           = DV_vec;
results.Tmax_vec         = Tmax_vec;
results.BUDGET_FACTOR    = BUDGET_FACTOR;
results.RMIN_GUARD       = RMIN_GUARD;
results.TIE_TOL_REL      = TIE_TOL_REL;
results.TIE_TOL_ABS      = TIE_TOL_ABS;
results.short_names      = short_names;
results.full_names       = full_names;
results.skip_map         = skip_map;
results.edges_kept_map   = edges_kept_map;
results.lcc_size_map     = lcc_size_map;
results.lcc_full_map     = lcc_full_map;
results.ap_count_map     = ap_count_map;
results.R_budget_all     = R_budget_all;
results.reach_all        = reach_all;
results.gateway_all      = gateway_all;
results.hub_all          = hub_all;
results.is_ap_all        = is_ap_all;
results.reach_winner_idx   = reach_winner_idx;
results.gateway_winner_idx = gateway_winner_idx;
results.hub_winner_idx     = hub_winner_idx;
results.reach_tie_sz       = reach_tie_sz;
results.gateway_tie_sz     = gateway_tie_sz;
results.hub_tie_sz         = hub_tie_sz;

save(mat_path, '-struct', 'results');

%% Figures ----------------------------------------------------------------
fprintf('Creating figures...\n');

% Assemble data structs for the plotting function
reach_map_data.winner_idx = reach_winner_idx;
reach_map_data.tie_size   = reach_tie_sz;

gw_map_data.winner_idx = gateway_winner_idx;
gw_map_data.tie_size   = gateway_tie_sz;

hub_map_data.winner_idx = hub_winner_idx;
hub_map_data.tie_size   = hub_tie_sz;

lcc_map_data.lcc_size = lcc_size_map;

ap_map_data.ap_count = ap_count_map;
ap_map_data.skip_map = skip_map;

net_plot_winner_map('winner_map_reach',    reach_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, short_names, N, OUT_DIR);

net_plot_winner_map('winner_map_gateway',  gw_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, short_names, N, OUT_DIR);

net_plot_winner_map('winner_map_hub',      hub_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, short_names, N, OUT_DIR);

net_plot_winner_map('lcc_map',             lcc_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, short_names, N, OUT_DIR);

net_plot_winner_map('articulation_map',    ap_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, short_names, N, OUT_DIR);

%% GIF animation (optional) -----------------------------------------------
if GIF_ENABLE
    fprintf('Generating GIF animation (dj=%d)...\n', GIF_TMAX_IDX);
    net_make_gif(DVmatrix_sweep, TOFmatrix_sweep, DV_cap_list, Tmax_list, ...
        GIF_TMAX_IDX, VU_mps, TU_days, BUDGET_FACTOR, ...
        short_names, OUT_DIR, GIF_FRAME_DELAY);
end

%% Print file summary -----------------------------------------------------
fprintf('\n=== Output files created in: %s ===\n', OUT_DIR);
created = { ...
    'snapshot_summary.csv', ...
    'node_metrics.csv',     ...
    'network_results.mat',  ...
    'winner_map_reach.png', 'winner_map_reach.fig',    ...
    'winner_map_gateway.png','winner_map_gateway.fig', ...
    'winner_map_hub.png',    'winner_map_hub.fig',     ...
    'lcc_map.png',           'lcc_map.fig',             ...
    'articulation_map.png',  'articulation_map.fig',   ...
    };
if GIF_ENABLE
    created{end+1} = sprintf('animation_Tmax%d.gif', GIF_TMAX_IDX);
end
for fi = 1:numel(created)
    p = fullfile(OUT_DIR, created{fi});
    if exist(p, 'file')
        fprintf('  ✓  %s\n', created{fi});
    else
        fprintf('  ✗  %s  (NOT FOUND — check for errors above)\n', created{fi});
    end
end
fprintf('Done.\n');

% =========================================================================
%  Local functions (MATLAB R2016b+)
% =========================================================================

function [winner_idx, tie_sz, name_str] = local_get_winner( ...
        metric, short_names, N, tol_rel, tol_abs, rmin_guard, R_budget)
%LOCAL_GET_WINNER  Step 10: apply tie tolerances and return winner info.
%
%   Returns
%     winner_idx  k (1..N) if unique winner, 0 if tie, NaN if no valid data
%     tie_sz      number of co-winners
%     name_str    short name, or semicolon-joined short names if tie

% Determine candidate set (Rmin_guard applies to reach only)
if rmin_guard > 0 && ~isempty(R_budget)
    candidates = find(R_budget(:) >= rmin_guard);
    if isempty(candidates)
        candidates = (1:N)';   % fallback: all nodes
    end
else
    candidates = (1:N)';
end

metric_cand = metric(candidates);

% Handle degenerate cases
if isempty(metric_cand) || all(~isfinite(metric_cand))
    winner_idx = NaN;
    tie_sz     = 0;
    name_str   = '';
    return
end

best_val = max(metric_cand(isfinite(metric_cand)));
tol      = max(tol_rel * abs(best_val), tol_abs);

co_idx   = candidates(metric_cand >= best_val - tol);   % 1-based node indices
tie_sz   = numel(co_idx);

if tie_sz == 0
    winner_idx = NaN;
    name_str   = '';
elseif tie_sz == 1
    winner_idx = co_idx(1);
    name_str   = short_names{winner_idx};
else
    winner_idx = 0;   % Tie marker
    name_str   = strjoin(short_names(co_idx)', ';');
end

end

function s = local_winner_str(name_str, winner_idx)
%LOCAL_WINNER_STR  Format the winner field for CSV output.
if isnan(winner_idx)
    s = 'skip';
elseif isempty(name_str)
    s = 'skip';
else
    s = name_str;   % already the short name or semicolon-joined tie list
end
end
