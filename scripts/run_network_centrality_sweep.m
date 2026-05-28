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
%    5-9. Computes budgeted reachability, harmonic closeness, betweenness,
%         strength, and articulation-point metrics
%    10. Resolves ties and records winners
%
%  Outputs (all under OUT_DIR):
%    snapshot_summary.csv, node_metrics.csv, network_results.mat
%    edges_count_map.{pdf,png,svg,fig}
%    lcc_map.{pdf,png,svg,fig}
%    articulation_map.{pdf,png,svg,fig}
%    winner_map_strength.{pdf,png,svg,fig}
%    winner_map_harmonic_closeness.{pdf,png,svg,fig}
%    winner_map_betweenness.{pdf,png,svg,fig}
%    strength_contour.{pdf,png,svg,fig}
%    harmonic_closeness_contour.{pdf,png,svg,fig}
%    budget_feasible_pairs_map.{pdf,png,svg,fig}
%    baseline_strength.{pdf,png,svg,fig}
%    baseline_harmonic_closeness.{pdf,png,svg,fig}
%    baseline_betweenness.{pdf,png,svg,fig}
%    [optional] animation_Tmax<dj>.gif
%
%  Dependencies: atlas_cfg_defaults, src/network/*.m  (added to path below)
% =========================================================================

%% ========================================================================
%  USER KNOBS — edit this block; do not change code below the divider
% =========================================================================

% Path to the sweep .mat produced by run_overlap_dv_tmax_sweep.m
SWEEP_MAT = fullfile(repo_root(), 'atlas_sweep_results', ...
                     'sweep_DVmatrix_results.mat');

% Output directory (created automatically if absent)
OUT_DIR = fullfile(repo_root(), 'atlas_network_results');

% Physical budget multiplier  (2 = departure manoeuvre + arrival manoeuvre)
BUDGET_FACTOR = 2;

% Tie tolerances
TIE_TOL_REL       = 1e-4;
TIE_TOL_ABS_SCALE = 0;

% Baseline snapshot for per-family metric bar charts.
% Set to [] to use the maximum-budget cell (last di, last dj).
BASELINE_DI = [];   % DV-cap index  ([] = max budget)
BASELINE_DJ = [];   % Tmax index    ([] = max budget)

% GIF animation settings
GIF_ENABLE      = false;
GIF_TMAX_IDX    = 10;
GIF_FRAME_DELAY = 0.5;

% =========================================================================
%  END OF USER KNOBS
% =========================================================================

%% Setup ------------------------------------------------------------------
setup();
addpath(fullfile(repo_root(), 'src', 'network'));

if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

cfg     = atlas_cfg_defaults();
VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;

TIE_TOL_ABS = TIE_TOL_ABS_SCALE * VU_mps;

%% Load sweep data --------------------------------------------------------
fprintf('Loading sweep data:\n  %s\n', SWEEP_MAT);
if ~exist(SWEEP_MAT, 'file')
    error('run_network_centrality_sweep: sweep .mat not found:\n  %s\n', SWEEP_MAT);
end
S = load(SWEEP_MAT, 'DV_cap_list', 'Tmax_list', 'DVmatrix_sweep', ...
         'TOFmatrix_sweep', 'families');

DV_cap_list     = S.DV_cap_list(:);
Tmax_list       = S.Tmax_list(:);
DVmatrix_sweep  = S.DVmatrix_sweep;
TOFmatrix_sweep = S.TOFmatrix_sweep;
full_names      = S.families(:);

nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);
N     = 13;

if numel(full_names) ~= N
    warning('Expected %d families in sweep .mat; found %d.', N, numel(full_names));
    N = numel(full_names);
end

DV_vec   = BUDGET_FACTOR * DV_cap_list * VU_mps;
Tmax_vec = BUDGET_FACTOR * Tmax_list   * TU_days;

% Max undirected pairs for N families
MAX_PAIRS = N * (N - 1) / 2;   % = 78

short_names = net_family_short_names();

fprintf('  Grid: %d DV-cap values × %d Tmax values = %d snapshots\n', ...
        nDV, nTmax, nDV*nTmax);

% Resolve baseline snapshot indices
if isempty(BASELINE_DI), BASELINE_DI = nDV;   end
if isempty(BASELINE_DJ), BASELINE_DJ = nTmax; end
BASELINE_DI = max(1, min(BASELINE_DI, nDV));
BASELINE_DJ = max(1, min(BASELINE_DJ, nTmax));

%% Preallocate result arrays ----------------------------------------------

skip_map       = true(nDV, nTmax);
edges_kept_map = zeros(nDV, nTmax);
lcc_size_map   = zeros(nDV, nTmax);
lcc_full_map   = zeros(nDV, nTmax);
ap_count_map   = zeros(nDV, nTmax);

R_budget_all = NaN(N, nDV, nTmax);
hc_all       = NaN(N, nDV, nTmax);
bw_all       = NaN(N, nDV, nTmax);
str_all      = NaN(N, nDV, nTmax);
is_ap_all    = false(N, nDV, nTmax);

hc_winner_idx  = NaN(nDV, nTmax);
bw_winner_idx  = NaN(nDV, nTmax);
str_winner_idx = NaN(nDV, nTmax);

hc_tie_sz  = zeros(nDV, nTmax);
bw_tie_sz  = zeros(nDV, nTmax);
str_tie_sz = zeros(nDV, nTmax);

hc_winner_names  = repmat({''}, nDV, nTmax);
bw_winner_names  = repmat({''}, nDV, nTmax);
str_winner_names = repmat({''}, nDV, nTmax);

%% Main sweep loop --------------------------------------------------------
fprintf('Processing snapshots...\n');
t_start = tic;

for di = 1:nDV
    for dj = 1:nTmax

        [A, W, D_sym, ~, edges_kept, DVcap_true, ~, skip] = ...
            net_build_graph( ...
                DVmatrix_sweep{di, dj}, TOFmatrix_sweep{di, dj}, ...
                DV_cap_list(di), Tmax_list(dj), ...
                VU_mps, TU_days, BUDGET_FACTOR);

        if skip, continue; end

        skip_map(di, dj)       = false;
        edges_kept_map(di, dj) = edges_kept;

        dist = net_floyd_warshall(W);

        [lcc_sz, lcc_full, ~] = net_lcc(A);
        lcc_size_map(di, dj)  = lcc_sz;
        lcc_full_map(di, dj)  = lcc_full;

        metrics = net_centrality(A, W, D_sym, dist, DVcap_true);

        R_budget_all(:, di, dj) = metrics.R_budget;
        hc_all(:,      di, dj)  = metrics.harmonic_closeness;
        bw_all(:,      di, dj)  = metrics.betweenness;
        str_all(:,     di, dj)  = metrics.strength;
        is_ap_all(:,   di, dj)  = metrics.is_articulation;
        ap_count_map(di, dj)    = sum(metrics.is_articulation);

        [hc_winner_idx(di,dj),  hc_tie_sz(di,dj),  hc_winner_names{di,dj}] = ...
            local_get_winner(metrics.harmonic_closeness, short_names, N, ...
                             TIE_TOL_REL, TIE_TOL_ABS);

        [bw_winner_idx(di,dj),  bw_tie_sz(di,dj),  bw_winner_names{di,dj}] = ...
            local_get_winner(metrics.betweenness, short_names, N, ...
                             TIE_TOL_REL, TIE_TOL_ABS);

        [str_winner_idx(di,dj), str_tie_sz(di,dj), str_winner_names{di,dj}] = ...
            local_get_winner(metrics.strength, short_names, N, ...
                             TIE_TOL_REL, TIE_TOL_ABS);
    end

    if mod(di, max(1, floor(nDV/5))) == 0
        fprintf('  di=%d/%d done (%.1f s elapsed)\n', di, nDV, toc(t_start));
    end
end

fprintf('Sweep complete in %.1f s.\n', toc(t_start));

%% Derived maps -----------------------------------------------------------
% Undirected direct-edge count (upper triangle only, max = MAX_PAIRS = 78)
direct_pairs_map = edges_kept_map / 2;

% Undirected budget-feasible-path count (max = MAX_PAIRS = 78)
% R_budget(k) = fraction of (N-1) destinations reachable from k
% → total reachable ordered pairs = sum(R_budget) * (N-1)
% → total reachable unordered pairs = sum(R_budget) * (N-1) / 2
budget_pairs_map = squeeze(sum(R_budget_all, 1, 'omitnan')) * (N - 1) / 2;
% [nDV×nTmax]; NaN where skip
budget_pairs_map(skip_map) = NaN;

%% Save snapshot_summary.csv ----------------------------------------------
csv1_path = fullfile(OUT_DIR, 'snapshot_summary.csv');
fprintf('Writing %s ...\n', csv1_path);

fid = fopen(csv1_path, 'w');
fprintf(fid, ['di,dj,DVcap_nd,Tmax_nd,DVcap_true_mps,Tmax_true_days,' ...
              'edges_kept,direct_pairs,lcc_size,lcc_full,' ...
              'budget_pairs,' ...
              'hc_winner,hc_tie_size,' ...
              'bw_winner,bw_tie_size,' ...
              'str_winner,str_tie_size,' ...
              'n_articulation_points\n']);

for di = 1:nDV
    for dj = 1:nTmax
        if skip_map(di, dj)
            fprintf(fid, '%d,%d,%.6g,%.6g,%.4f,%.4f,%d,%d,%d,%d,%d,%s,%d,%s,%d,%s,%d,%d\n', ...
                di, dj, DV_cap_list(di), Tmax_list(dj), DV_vec(di), Tmax_vec(dj), ...
                0, 0, 0, 0, 0, 'skip', 0, 'skip', 0, 'skip', 0, 0);
            continue
        end
        fprintf(fid, '%d,%d,%.6g,%.6g,%.4f,%.4f,%d,%d,%d,%d,%.1f,%s,%d,%s,%d,%s,%d,%d\n', ...
            di, dj, DV_cap_list(di), Tmax_list(dj), DV_vec(di), Tmax_vec(dj), ...
            edges_kept_map(di,dj), direct_pairs_map(di,dj), ...
            lcc_size_map(di,dj), lcc_full_map(di,dj), ...
            budget_pairs_map(di,dj), ...
            local_winner_str(hc_winner_names{di,dj},  hc_winner_idx(di,dj)),  hc_tie_sz(di,dj), ...
            local_winner_str(bw_winner_names{di,dj},  bw_winner_idx(di,dj)),  bw_tie_sz(di,dj), ...
            local_winner_str(str_winner_names{di,dj}, str_winner_idx(di,dj)), str_tie_sz(di,dj), ...
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
               'R_budget,harmonic_closeness,betweenness,strength,is_articulation\n']);

for di = 1:nDV
    for dj = 1:nTmax
        if skip_map(di, dj), continue; end
        for k = 1:N
            fprintf(fid2, '%d,%d,%.4f,%.4f,%d,"%s",%.6g,%.6g,%.6g,%.6g,%d\n', ...
                di, dj, DV_vec(di), Tmax_vec(dj), k, full_names{k}, ...
                R_budget_all(k, di, dj), hc_all(k, di, dj), ...
                bw_all(k, di, dj), str_all(k, di, dj), ...
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
results.TIE_TOL_REL      = TIE_TOL_REL;
results.TIE_TOL_ABS      = TIE_TOL_ABS;
results.short_names      = short_names;
results.full_names       = full_names;
results.skip_map         = skip_map;
results.edges_kept_map   = edges_kept_map;
results.direct_pairs_map = direct_pairs_map;
results.budget_pairs_map = budget_pairs_map;
results.lcc_size_map     = lcc_size_map;
results.lcc_full_map     = lcc_full_map;
results.ap_count_map     = ap_count_map;
results.R_budget_all     = R_budget_all;
results.harmonic_closeness_all = hc_all;
results.betweenness_all  = bw_all;
results.strength_all     = str_all;
results.is_ap_all        = is_ap_all;
results.hc_winner_idx    = hc_winner_idx;
results.bw_winner_idx    = bw_winner_idx;
results.str_winner_idx   = str_winner_idx;
results.hc_tie_sz        = hc_tie_sz;
results.bw_tie_sz        = bw_tie_sz;
results.str_tie_sz       = str_tie_sz;
results.hc_winner_names  = hc_winner_names;
results.bw_winner_names  = bw_winner_names;
results.str_winner_names = str_winner_names;
results.BASELINE_DI      = BASELINE_DI;
results.BASELINE_DJ      = BASELINE_DJ;

save(mat_path, '-struct', 'results');

%% Figures ----------------------------------------------------------------
fprintf('Creating figures...\n');

% ── Fig 1: Direct family-pair count ──────────────────────────────────────────
edge_map_data.direct_pairs = direct_pairs_map;
edge_map_data.skip_map     = skip_map;
net_plot_winner_map('edges_count_map', edge_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Fig 2: LCC size ──────────────────────────────────────────────────────────
lcc_map_data.lcc_size = lcc_size_map;
net_plot_winner_map('lcc_map', lcc_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Fig 3: Articulation-point count ──────────────────────────────────────────
ap_map_data.ap_count = ap_count_map;
ap_map_data.skip_map = skip_map;
net_plot_winner_map('articulation_map', ap_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Fig 4: Strength winner map ────────────────────────────────────────────────
str_map_data.winner_idx   = str_winner_idx;
str_map_data.tie_size     = str_tie_sz;
str_map_data.winner_names = str_winner_names;
net_plot_winner_map('winner_map_strength', str_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Fig 5: Harmonic closeness winner map ─────────────────────────────────────
hc_map_data.winner_idx   = hc_winner_idx;
hc_map_data.tie_size     = hc_tie_sz;
hc_map_data.winner_names = hc_winner_names;
net_plot_winner_map('winner_map_harmonic_closeness', hc_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Fig 6: Betweenness winner map ────────────────────────────────────────────
bw_map_data.winner_idx   = bw_winner_idx;
bw_map_data.tie_size     = bw_tie_sz;
bw_map_data.winner_names = bw_winner_names;
net_plot_winner_map('winner_map_betweenness', bw_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Fig 7: Strength continuous heatmap ───────────────────────────────────────
str_ctr_data.values   = str_all;
str_ctr_data.skip_map = skip_map;
net_plot_winner_map('strength_contour', str_ctr_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Fig 8: Harmonic closeness continuous heatmap ─────────────────────────────
hc_ctr_data.values   = hc_all;
hc_ctr_data.skip_map = skip_map;
net_plot_winner_map('harmonic_closeness_contour', hc_ctr_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Fig 9: Budget-feasible transfer pairs ────────────────────────────────────
bp_map_data.budget_pairs = budget_pairs_map;
bp_map_data.skip_map     = skip_map;
net_plot_winner_map('budget_feasible_pairs_map', bp_map_data, ...
    DV_vec, Tmax_vec, lcc_full_map, budget_pairs_map, short_names, N, OUT_DIR);

% ── Figs 10-12: Baseline per-family bar charts ───────────────────────────────
baseline_data.strength            = str_all(:, BASELINE_DI, BASELINE_DJ);
baseline_data.harmonic_closeness  = hc_all(:,  BASELINE_DI, BASELINE_DJ);
baseline_data.betweenness         = bw_all(:,  BASELINE_DI, BASELINE_DJ);
baseline_data.is_articulation     = is_ap_all(:, BASELINE_DI, BASELINE_DJ);
baseline_data.DVcap_mps           = DV_vec(BASELINE_DI);
baseline_data.Tmax_days           = Tmax_vec(BASELINE_DJ);

net_plot_baseline(baseline_data, short_names, N, OUT_DIR);

%% GIF animation (optional) -----------------------------------------------
if GIF_ENABLE
    fprintf('Generating GIF animation (dj=%d)...\n', GIF_TMAX_IDX);
    net_make_gif(DVmatrix_sweep, TOFmatrix_sweep, DV_cap_list, Tmax_list, ...
        GIF_TMAX_IDX, VU_mps, TU_days, BUDGET_FACTOR, ...
        short_names, OUT_DIR, GIF_FRAME_DELAY);
end

%% Print file summary -----------------------------------------------------
fprintf('\n=== Output files in: %s ===\n', OUT_DIR);

sweep_figs = { ...
    'edges_count_map',              ...
    'lcc_map',                      ...
    'articulation_map',             ...
    'winner_map_strength',          ...
    'winner_map_harmonic_closeness',...
    'winner_map_betweenness',       ...
    'strength_contour',             ...
    'harmonic_closeness_contour',   ...
    'budget_feasible_pairs_map',    ...
    };
baseline_figs = { ...
    'baseline_strength',            ...
    'baseline_harmonic_closeness',  ...
    'baseline_betweenness',         ...
    };
data_files = {'snapshot_summary.csv', 'node_metrics.csv', 'network_results.mat'};
if GIF_ENABLE
    data_files{end+1} = sprintf('animation_Tmax%d.gif', GIF_TMAX_IDX);
end

all_figs = [sweep_figs, baseline_figs];
for fi = 1:numel(all_figs)
    for ext = {'.pdf', '.png', '.svg', '.fig'}
        p = fullfile(OUT_DIR, [all_figs{fi} ext{1}]);
        if exist(p, 'file'), fprintf('  OK  %s%s\n', all_figs{fi}, ext{1});
        else,                fprintf('  --  %s%s  (missing)\n', all_figs{fi}, ext{1}); end
    end
end
for fi = 1:numel(data_files)
    p = fullfile(OUT_DIR, data_files{fi});
    if exist(p, 'file'), fprintf('  OK  %s\n', data_files{fi});
    else,                fprintf('  --  %s  (missing)\n', data_files{fi}); end
end
fprintf('Done.\n');

% =========================================================================
%  Local functions
% =========================================================================

function [winner_idx, tie_sz, name_str] = local_get_winner( ...
        metric, short_names, N, tol_rel, tol_abs)

metric_cand = metric(:);
if isempty(metric_cand) || all(~isfinite(metric_cand))
    winner_idx = NaN;  tie_sz = 0;  name_str = '';  return
end

best_val = max(metric_cand(isfinite(metric_cand)));
tol      = max(tol_rel * abs(best_val), tol_abs);
co_idx   = find(metric_cand >= best_val - tol);
tie_sz   = numel(co_idx);

if tie_sz == 0
    winner_idx = NaN;  name_str = '';
elseif tie_sz == 1
    winner_idx = co_idx(1);
    name_str   = short_names{winner_idx};
else
    winner_idx = 0;
    name_str   = strjoin(short_names(co_idx)', ';');
end
end

function s = local_winner_str(name_str, winner_idx)
if isnan(winner_idx) || isempty(name_str)
    s = 'skip';
else
    s = name_str;
end
end
