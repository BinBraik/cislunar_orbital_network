%% RUN_OVERLAP_MINTOF_MAXBUDGET_PREVIEW
% Single-snapshot preview of the min-TOF network, at the MAXIMUM budget
% (DV_cap_nd = 0.2, Tmax = pi — i.e. exactly the base/full atlas config,
% no subset derivation needed). Meant to be reviewed BEFORE committing to
% the full (DV_cap, Tmax) sweep in run_overlap_dv_tmax_sweep.m.
%
% For every family pair, scores overlap voxels two independent ways
% (src/overlap_proxy_pair.m), each using per-voxel MIN TOF (not mean —
% see src/overlap_proxy_footprint.m):
%   - min-DV-proxy winner  → today's existing network (baseline, for comparison)
%   - min-TOF-proxy winner → the new min-TOF network
%
% Both networks are then built with the same feasibility/centrality pipeline
% as run_network_centrality_sweep.m (net_build_graph / net_floyd_warshall /
% net_centrality), for this one snapshot only.
%
% PRE-FLIGHT: this script REQUIRES a cached atlas for every family already
% present under cfg.cache.dir. It will NOT silently rebuild — if any atlas
% is missing it errors out immediately (rebuilding one atlas can take
% hours; that should be an explicit, separate step).
%
% Outputs written to OUT_DIR:
%   node_metrics_maxbudget.csv   — per-family metrics, both networks side by side
%   network_compare_maxbudget.mat
%   mindv/baseline_*.{pdf,png,svg,fig}   — existing min-DV network bar charts
%   mintof/baseline_*.{pdf,png,svg,fig}  — new min-TOF network bar charts

clear; clc;

% ── repo paths ────────────────────────────────────────────────────────────────
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'src', 'network'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% ══════════════════════════════════════════════════════════════════════════════
%  USER KNOBS
% ══════════════════════════════════════════════════════════════════════════════

OUT_DIR       = fullfile(repoRoot, 'atlas_network_results_mintof_preview');
BUDGET_FACTOR = 2;     % departure + arrival manoeuvre, matches run_network_centrality_sweep.m

% Families — must match the base cached atlases.
families = { ...
    'Lyapunov L1', ...
    'Lyapunov L2', ...
    'Cycler 21', ...
    'Cycler 11a', ...
    'Cycler 11b', ...
    'Cycler 32', ...
    'Resonant 2to1 Stable', ...
    'Resonant 2to1 Unstable', ...
    'Resonant 3to1 Stable', ...
    'Resonant 3to1 Unstable', ...
    'Resonant 5to2 Stable', ...
    'Resonant 5to2 Unstable', ...
    'Distant Prograde Orbit' ...
};
N = numel(families);

% ══════════════════════════════════════════════════════════════════════════════
%  BASE (= MAXIMUM BUDGET) CONFIGURATION  — must match your cached atlases exactly
% ══════════════════════════════════════════════════════════════════════════════
cfg = atlas_cfg_defaults();

cfg.families.list      = families;
cfg.families.test_only = false;

cfg.grid.dx               = 0.001;
cfg.grid.dy               = 0.001;
cfg.grid.dtheta           = deg2rad(1);
cfg.seed.ds_seed          = 0.01;
cfg.propag.Tmax           = pi;      % maximum budget
cfg.fan.DV_cap_nd         = 0.2;     % maximum budget
cfg.fan.dtheta_fan        = deg2rad(0.5);
cfg.propag.absTol         = 1e-8;
cfg.propag.relTol         = 1e-8;
cfg.propag.v2tol          = 1e-8;
cfg.log.step_len_factor   = 0.75;
cfg.log.maxstep_factor    = 2;

cfg.cache.enable  = true;
cfg.cache.dir     = fullfile(repoRoot, 'atlas_cache');
cfg.cache.rebuild = false;

cfg.io.save_figs   = false;
cfg.io.save_fig    = false;
cfg.io.fig_visible = 'off';

cfg.plot.overlap.overlap_xy   = false;
cfg.plot.overlap.overlap_xyz  = false;
cfg.plot.overlap.combo_xy     = false;
cfg.plot.overlap.combo_xyz    = false;
cfg.plot.overlap.bounds_lb    = false;
cfg.plot.overlap.bounds_ub    = false;
cfg.plot.overlap.bounds_proxy = false;

if exist('atlas_cfg_validate', 'file') == 2
    atlas_cfg_validate(cfg);
end

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
DV_cap_nd = cfg.fan.DV_cap_nd;
Tmax_nd   = cfg.propag.Tmax;

if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

% ══════════════════════════════════════════════════════════════════════════════
%  PRE-FLIGHT: verify the atlas cache exists for every family
% ══════════════════════════════════════════════════════════════════════════════
[cache_ok, cache_report] = atlas_check_cache_exists(families, cfg); %#ok<ASGLU>
if ~cache_ok
    error(['run_overlap_mintof_maxbudget_preview: one or more cached atlases ' ...
           'are missing (see list above). This preview script will not build ' ...
           'atlases from scratch — run the base atlas build for those ' ...
           'families first (e.g. via run_atlas_one_family.m), then re-run.']);
end

% ══════════════════════════════════════════════════════════════════════════════
%  LOAD ATLASES  (cache hit expected for all — checked above)
% ══════════════════════════════════════════════════════════════════════════════
fprintf('[preview] Loading atlases (Tmax=%.4f, DV_cap_nd=%.4f)...\n', Tmax_nd, DV_cap_nd);
grid3 = atlas_grid_make(cfg);
Sall  = cell(N, 1);
for i = 1:N
    fprintf('[preview]   family %d/%d: %s\n', i, N, families{i});
    [Sall{i}, cacheInfo] = atlas_prepare_or_load(families{i}, cfg, grid3);
    if ~cacheInfo.hit
        warning('[preview] Family "%s" was NOT a cache hit — it was just rebuilt from scratch.', ...
            families{i});
    end
end

% ══════════════════════════════════════════════════════════════════════════════
%  BUILD PER-FAMILY VOXEL FOOTPRINTS  (min DV, min TOF per voxel)
% ══════════════════════════════════════════════════════════════════════════════
fprintf('[preview] Building voxel footprints...\n');
Fall = cell(N, 1);
for i = 1:N
    Fall{i} = overlap_proxy_footprint(Sall{i}, grid3, VU_mps, TU_days);
end
clear Sall

% ══════════════════════════════════════════════════════════════════════════════
%  PAIR LOOP — both min-DV-proxy and min-TOF-proxy winners
% ══════════════════════════════════════════════════════════════════════════════
nPairs = N * (N - 1) / 2;
pairI = zeros(nPairs, 1);  pairJ = zeros(nPairs, 1);
p = 0;
for ii = 1:N
    for jj = ii+1:N
        p = p + 1;
        pairI(p) = ii;  pairJ(p) = jj;
    end
end

pair_minDV      = nan(nPairs, 1);
pair_DVlb       = nan(nPairs, 1);
pair_DVpatch    = nan(nPairs, 1);
pair_TOF        = nan(nPairs, 1);
pair_voxelId    = nan(nPairs, 1);
pair_minTOF     = nan(nPairs, 1);
pair_DVatMinTOF = nan(nPairs, 1);
pair_voxelIdTOF = nan(nPairs, 1);

fprintf('[preview] Running %d pairwise overlaps...\n', nPairs);
for p = 1:nPairs
    [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), pair_TOF(p), pair_voxelId(p), ...
        pair_minTOF(p), pair_DVatMinTOF(p), pair_voxelIdTOF(p)] = ...
        overlap_proxy_pair(Fall{pairI(p)}, Fall{pairJ(p)}, grid3, cfg, VU_mps);
end
clear Fall

minDVproxyMat  = nan(N, N);
TOFatMinDVmat  = nan(N, N);
minTOFproxyMat = nan(N, N);
DVatMinTOFmat  = nan(N, N);
for p = 1:nPairs
    i = pairI(p);  j = pairJ(p);
    minDVproxyMat(i, j)  = pair_minDV(p);   minDVproxyMat(j, i)  = pair_minDV(p);
    TOFatMinDVmat(i, j)  = pair_TOF(p);     TOFatMinDVmat(j, i)  = pair_TOF(p);
    minTOFproxyMat(i, j) = pair_minTOF(p);  minTOFproxyMat(j, i) = pair_minTOF(p);
    DVatMinTOFmat(i, j)  = pair_DVatMinTOF(p); DVatMinTOFmat(j, i) = pair_DVatMinTOF(p);
end

% ══════════════════════════════════════════════════════════════════════════════
%  BUILD BOTH NETWORKS  (same feasibility + centrality pipeline as the sweep)
% ══════════════════════════════════════════════════════════════════════════════
short_names = net_family_short_names();

fprintf('[preview] Building min-DV network (baseline)...\n');
[A_dv, W_dv, D_sym_dv, T_sym_dv, edges_dv, DVcap_true, Tmax_true, skip_dv] = ...
    net_build_graph(minDVproxyMat, TOFatMinDVmat, DV_cap_nd, Tmax_nd, ...
                     VU_mps, TU_days, BUDGET_FACTOR);
if skip_dv, error('[preview] min-DV snapshot is empty/all-NaN — check atlas data.'); end
dist_dv    = net_floyd_warshall(W_dv);
metrics_dv = net_centrality(A_dv, W_dv, D_sym_dv, dist_dv, DVcap_true);
[lcc_sz_dv, lcc_full_dv, ~] = net_lcc(A_dv);

fprintf('[preview] Building min-TOF network (new)...\n');
[A_tof, ~, D_sym_tof, T_sym_tof, edges_tof, ~, ~, skip_tof] = ...
    net_build_graph(DVatMinTOFmat, minTOFproxyMat, DV_cap_nd, Tmax_nd, ...
                     VU_mps, TU_days, BUDGET_FACTOR);
if skip_tof, error('[preview] min-TOF snapshot is empty/all-NaN — check atlas data.'); end

% TOF-weighted adjacency: same feasibility mask (A_tof), edge cost = TOF (days)
W_tof = Inf(N);
feasMask = logical(A_tof);
W_tof(feasMask)      = T_sym_tof(feasMask);
W_tof(1:N+1:end)     = 0;

dist_tof    = net_floyd_warshall(W_tof);
metrics_tof = net_centrality(A_tof, W_tof, T_sym_tof, dist_tof, Tmax_true);
[lcc_sz_tof, lcc_full_tof, ~] = net_lcc(A_tof);

% ══════════════════════════════════════════════════════════════════════════════
%  CONSOLE SUMMARY
% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n=== Snapshot: DV_cap = %.1f m/s, Tmax = %.1f days ===\n', DVcap_true, Tmax_true);

local_print_summary('MIN-DV network  (baseline)', short_names, N, ...
    edges_dv, lcc_sz_dv, lcc_full_dv, metrics_dv);
local_print_summary('MIN-TOF network (new)', short_names, N, ...
    edges_tof, lcc_sz_tof, lcc_full_tof, metrics_tof);

% ══════════════════════════════════════════════════════════════════════════════
%  SAVE node_metrics_maxbudget.csv  (both networks, side by side)
% ══════════════════════════════════════════════════════════════════════════════
csv_path = fullfile(OUT_DIR, 'node_metrics_maxbudget.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, ['node_idx,node_name,' ...
    'dv_R_budget,dv_harmonic_closeness,dv_betweenness,dv_strength,dv_is_articulation,' ...
    'tof_R_budget,tof_harmonic_closeness,tof_betweenness,tof_strength,tof_is_articulation\n']);
for k = 1:N
    fprintf(fid, '%d,"%s",%.6g,%.6g,%.6g,%.6g,%d,%.6g,%.6g,%.6g,%.6g,%d\n', ...
        k, families{k}, ...
        metrics_dv.R_budget(k), metrics_dv.harmonic_closeness(k), ...
        metrics_dv.betweenness(k), metrics_dv.strength(k), int8(metrics_dv.is_articulation(k)), ...
        metrics_tof.R_budget(k), metrics_tof.harmonic_closeness(k), ...
        metrics_tof.betweenness(k), metrics_tof.strength(k), int8(metrics_tof.is_articulation(k)));
end
fclose(fid);
fprintf('\n[preview] Wrote %s\n', csv_path);

% ══════════════════════════════════════════════════════════════════════════════
%  SAVE .mat bundle
% ══════════════════════════════════════════════════════════════════════════════
mat_path = fullfile(OUT_DIR, 'network_compare_maxbudget.mat');
save(mat_path, 'families', 'short_names', 'DV_cap_nd', 'Tmax_nd', ...
    'DVcap_true', 'Tmax_true', 'BUDGET_FACTOR', ...
    'minDVproxyMat', 'TOFatMinDVmat', 'minTOFproxyMat', 'DVatMinTOFmat', ...
    'A_dv', 'W_dv', 'dist_dv', 'metrics_dv', 'edges_dv', 'lcc_sz_dv', 'lcc_full_dv', ...
    'A_tof', 'W_tof', 'dist_tof', 'metrics_tof', 'edges_tof', 'lcc_sz_tof', 'lcc_full_tof', ...
    '-v7.3');
fprintf('[preview] Wrote %s\n', mat_path);

% ══════════════════════════════════════════════════════════════════════════════
%  FIGURES  (per-family bar charts, one set per network)
% ══════════════════════════════════════════════════════════════════════════════
dv_dir  = fullfile(OUT_DIR, 'mindv');
tof_dir = fullfile(OUT_DIR, 'mintof');

bd_dv.strength           = metrics_dv.strength;
bd_dv.harmonic_closeness = metrics_dv.harmonic_closeness;
bd_dv.betweenness        = metrics_dv.betweenness;
bd_dv.is_articulation    = metrics_dv.is_articulation;
bd_dv.DVcap_mps          = DVcap_true;
bd_dv.Tmax_days          = Tmax_true;
net_plot_baseline(bd_dv, short_names, N, dv_dir);

bd_tof = bd_dv;
bd_tof.strength           = metrics_tof.strength;
bd_tof.harmonic_closeness = metrics_tof.harmonic_closeness;
bd_tof.betweenness        = metrics_tof.betweenness;
bd_tof.is_articulation    = metrics_tof.is_articulation;
net_plot_baseline(bd_tof, short_names, N, tof_dir);

fprintf('[preview] Wrote bar charts to:\n  %s\n  %s\n', dv_dir, tof_dir);
fprintf('\n[preview] Done. Review the above before running the full sweep:\n');
fprintf('  scripts/run_overlap_dv_tmax_sweep.m\n');
fprintf('  scripts/run_network_centrality_sweep.m  (min-DV, existing)\n');

% ══════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════════

function local_print_summary(label, short_names, N, edges_kept, lcc_sz, lcc_full, metrics)
fprintf('\n--- %s ---\n', label);
fprintf('  Feasible edges (undirected) : %d / %d\n', edges_kept/2, N*(N-1)/2);
fprintf('  LCC size                    : %d / %d  (full=%d)\n', lcc_sz, N, lcc_full);
fprintf('  Articulation points         : %d\n', sum(metrics.is_articulation));

[~, ord_s]  = sort(metrics.strength, 'descend');
[~, ord_hc] = sort(metrics.harmonic_closeness, 'descend');
[~, ord_bw] = sort(metrics.betweenness, 'descend');

top = @(ord) strjoin(short_names(ord(1:min(3,N)))', ', ');
fprintf('  Top-3 strength              : %s\n', top(ord_s));
fprintf('  Top-3 harmonic closeness    : %s\n', top(ord_hc));
fprintf('  Top-3 betweenness           : %s\n', top(ord_bw));
end
