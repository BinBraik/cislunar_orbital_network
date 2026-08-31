%% RUN_RS4_DV_TMAX_SWEEP
% Sweep DVmatrix/TOFmatrix (min-DV-proxy) AND TOFproxyMatrix/DVatMinTOFMatrix
% (min-TOF-proxy) over all (DV_cap_nd, Tmax) combinations.
%
% Each family pair's overlap voxels are scored two independent ways
% (see src/overlap_proxy_pair.m):
%   - min-DV-proxy winner:  argmin(dv_proxy)   → DVmatrix_sweep / TOFmatrix_sweep
%   - min-TOF-proxy winner: argmin(tof_proxy)  → TOFproxyMatrix_sweep / DVatMinTOFMatrix_sweep
% where tof_proxy = per-voxel MIN TOF of family A + MIN TOF of family B
% (src/overlap_proxy_footprint.m), matching how dv_proxy already uses MIN DV.
%
% Strategy:
%   0. Pre-flight: verify a cached atlas exists for every family (fails fast
%      instead of silently rebuilding for hours).
%   1. Load base atlases once  (Tmax = 3*pi/2, DV_cap = 0.3).
%   2. For each (DV_cap, Tmax) cell:
%        a. Check atlas_results for a finished run whose config matches
%           → extract minDVproxyMat + TOFmatrix directly, skip recompute.
%           NOTE: legacy atlas_results runs predate the min-TOF-proxy
%           selection, so reused cells get NaN TOFproxyMatrix/DVatMinTOFMatrix.
%        b. Otherwise: derive in-memory subset atlases and run the pairwise
%           overlap loop  (no figures, no per-pair .mat files).
%   3. Checkpoint after EVERY completed cell  (safe to kill + requeue on HPC).
%   4. Write final .mat  +  three Excel workbooks when all cells are done.
%
% Parallelism:
%   Serial over (DV_cap, Tmax) combinations  (predictable memory footprint).
%   Optional parfor over the N*(N-1)/2 pairs within each combination.
%   Set N_WORKERS = 0 to force fully serial (safest for memory-limited nodes).
%
% Outputs written to OUTPUT_DIR:
%   sweep_DVmatrix_results.mat  — authoritative data store (DV-proxy AND TOF-proxy)
%   sweep_DVmatrix.xlsx         — one sheet per combination, N×N DVproxy matrix
%   sweep_TOFmatrix.xlsx        — one sheet per combination, N×N mean-TOF matrix
%   sweep_winners.xlsx          — one sheet per combination, pair-winner table
%
% See also: scripts/run_overlap_mintof_maxbudget_preview.m — single-snapshot
% (max-budget) preview of the min-TOF network, meant to be inspected BEFORE
% committing to a full sweep run.

clear; clc;

% ── repo paths ────────────────────────────────────────────────────────────────
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% ══════════════════════════════════════════════════════════════════════════════
%  USER KNOBS
% ══════════════════════════════════════════════════════════════════════════════

% Parallel workers for the PAIR loop inside each combination.
%   0 → fully serial  (use on memory-limited nodes)
%   N → parfor with N workers
N_WORKERS = 4;


% Output controls
% NOTE: With a fine sweep grid, the per-combination Excel-sheet approach can create
% hundreds of sheets and become slow/heavy. The .mat file is the authoritative output.
WRITE_EXCEL    = false;   % true → write Excel workbooks (one sheet per combination)
WRITE_FLAT_CSV = true;    % true → write a single flat CSV from winners_sweep
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

% HPC paths
CHECKPOINT_FILE = fullfile(repoRoot, 'atlas_sweep_checkpoint.mat');
OUTPUT_DIR      = fullfile(repoRoot, 'atlas_sweep_results');
SWEEP_CACHE_DIR = fullfile(repoRoot, 'atlas_cache_sweep');   % derived caches (if ever saved)

% ── Sweep grid ────────────────────────────────────────────────────────────────
% Upper bounds MUST NOT exceed the base atlas's own generation bounds
% (cfg.fan.DV_cap_nd / cfg.propag.Tmax below) — atlas_derive_subset can only
% shrink the candidate row set, never grow it beyond what the atlas has.
DV_cap_list = linspace(0.025, 0.300,   20)';   % 20 values  (0.025 → 0.300)
Tmax_list   = linspace(pi/4,  3*pi/2,  20)';   % 20 values  (pi/4  → 3*pi/2)

% Short labels for Excel sheet names  (≤31 chars; keep them compact)
% Encode Tmax as a multiple of pi, rounded to 3 decimals: e.g., Tp0p250 = 0.250*pi
Tmax_labels = arrayfun(@(t) strrep(sprintf('Tp%0.3f', t/pi), '.', 'p'), ...
                       Tmax_list, 'UniformOutput', false);

% ══════════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION  (must match your full cached atlases exactly)
% ══════════════════════════════════════════════════════════════════════════════
cfg = atlas_cfg_defaults();

cfg.families.list      = families;
cfg.families.test_only = false;

% NOTE: these are the extended-budget atlas config (Rdom=1.2, dx=dy=0.001,
% dtheta=1°, ds_seed=0.01, dtheta_fan=0.5°, DV_a=0.3, Ta=3*pi/2,
% tolerances=1e-8 — Table 3 of the paper except DV_a/Ta raised beyond their
% nominal 0.2/pi) and MUST byte-for-byte match TARGET_CFG in
% run_atlas_cache_rebuild_to_target.m (the script that reconciles/rebuilds
% atlas_cache/ to this config) — they feed atlas_grid_make and the cache
% fingerprint. If you change one, change it in both places.
cfg.grid.dx               = 0.001;
cfg.grid.dy               = 0.001;
cfg.grid.dtheta           = deg2rad(1);
cfg.seed.ds_seed          = 0.01;
cfg.propag.Tmax           = 3*pi/2;
cfg.fan.DV_cap_nd         = 0.3;
cfg.fan.dtheta_fan        = deg2rad(0.5);
cfg.propag.absTol         = 1e-8;
cfg.propag.relTol         = 1e-8;
cfg.propag.v2tol          = 1e-8;
cfg.log.step_len_factor   = 0.75;
cfg.log.maxstep_factor    = 2;

cfg.cache.enable      = true;
cfg.cache.dir         = fullfile(repoRoot, 'atlas_cache');
cfg.cache.rebuild     = false;
cfg.cache.version_tag = 'atlas_v1_keep_masked';   % must match run_atlas_cache_rebuild_to_target.m's TARGET_CFG

% Suppress all figure/file output from sub-functions
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

% ── Derived constants ─────────────────────────────────────────────────────────
VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
N      = numel(families);
nDV    = numel(DV_cap_list);
nTmax  = numel(Tmax_list);
nPairs = N * (N - 1) / 2;

% Pre-compute flat pair-index vectors (needed for parfor slicing)
pairI = zeros(nPairs, 1);
pairJ = zeros(nPairs, 1);
p = 0;
for ii = 1:N
    for jj = ii+1:N
        p = p + 1;
        pairI(p) = ii;
        pairJ(p) = jj;
    end
end

% Ensure output directories exist
if ~exist(OUTPUT_DIR, 'dir'),      mkdir(OUTPUT_DIR);      end
if ~exist(SWEEP_CACHE_DIR, 'dir'), mkdir(SWEEP_CACHE_DIR); end

fprintf('[sweep] Sweep grid: %d DV_cap values × %d Tmax values = %d combinations.\n', ...
    nDV, nTmax, nDV * nTmax);
fprintf('[sweep] Checkpoint file: %s\n', CHECKPOINT_FILE);

% ══════════════════════════════════════════════════════════════════════════════
%  LOAD / INITIALISE CHECKPOINT
% ══════════════════════════════════════════════════════════════════════════════
done_mask             = false(nDV, nTmax);
DVmatrix_sweep        = cell(nDV, nTmax);   % min-DV-proxy winner:  DV at that voxel
TOFmatrix_sweep       = cell(nDV, nTmax);   % TOF read off at the min-DV winner voxel
TOFproxyMatrix_sweep  = cell(nDV, nTmax);   % min-TOF-proxy winner: TOF at that voxel
DVatMinTOFMatrix_sweep = cell(nDV, nTmax);  % DV read off at the min-TOF winner voxel
winners_sweep         = cell(nDV, nTmax);
source_sweep          = cell(nDV, nTmax);

if isfile(CHECKPOINT_FILE)
    fprintf('[sweep] Checkpoint found — loading...\n');
    try
        ck = load(CHECKPOINT_FILE, ...
            'done_mask', 'DVmatrix_sweep', 'TOFmatrix_sweep', ...
            'TOFproxyMatrix_sweep', 'DVatMinTOFMatrix_sweep', ...
            'winners_sweep', 'source_sweep');
        if isequal(size(ck.done_mask), [nDV, nTmax])
            done_mask       = ck.done_mask;
            DVmatrix_sweep  = ck.DVmatrix_sweep;
            TOFmatrix_sweep = ck.TOFmatrix_sweep;
            if isfield(ck, 'TOFproxyMatrix_sweep')
                TOFproxyMatrix_sweep   = ck.TOFproxyMatrix_sweep;
                DVatMinTOFMatrix_sweep = ck.DVatMinTOFMatrix_sweep;
            end
            winners_sweep   = ck.winners_sweep;
            source_sweep    = ck.source_sweep;
            fprintf('[sweep] Checkpoint loaded. %d/%d cells already done.\n', ...
                sum(done_mask(:)), nDV * nTmax);
        else
            warning('[sweep] Checkpoint grid size mismatch — starting fresh.');
        end
    catch ME
        warning('[sweep] Could not load checkpoint (%s) — starting fresh.', ME.message);
    end
end

% ══════════════════════════════════════════════════════════════════════════════
%  SCAN atlas_results FOR REUSABLE COMPLETED RUNS
% ══════════════════════════════════════════════════════════════════════════════
fprintf('[sweep] Scanning atlas_results for matching completed runs...\n');
existing = local_scan_results(cfg.io.out_root, families, ...
    DV_cap_list, Tmax_list, cfg);

nFound = 0;
for di = 1:nDV
    for dj = 1:nTmax
        if done_mask(di, dj), continue; end
        if ~isempty(existing{di, dj})
            e = existing{di, dj};
            DVmatrix_sweep{di, dj}  = e.minDVproxyMat;
            TOFmatrix_sweep{di, dj} = e.TOFmat;
            % Legacy atlas_results runs predate the min-TOF-proxy selection
            % (overlap_proxy_pair) and only recorded the min-DV winner, so
            % there is no min-TOF voxel to recover here. Leave NaN — these
            % cells will read as "skip" in the min-TOF network unless
            % recomputed (delete the matching entry from CHECKPOINT_FILE,
            % or clear atlas_results, to force a recompute).
            TOFproxyMatrix_sweep{di, dj}   = nan(size(e.minDVproxyMat));
            DVatMinTOFMatrix_sweep{di, dj} = nan(size(e.minDVproxyMat));
            winners_sweep{di, dj}   = e.T;
            source_sweep{di, dj}    = 'loaded_from_results';
            done_mask(di, dj)       = true;
            nFound = nFound + 1;
            fprintf('[sweep]   loaded (%d,%d): DV=%.3f  Tmax=%s  (min-TOF-proxy NOT available for reused cells)\n', ...
                di, dj, DV_cap_list(di), local_tmax_str(Tmax_list(dj)));
        end
    end
end
fprintf('[sweep] %d cell(s) loaded from existing results. %d remain to compute.\n', ...
    nFound, sum(~done_mask(:)));

local_save_checkpoint(CHECKPOINT_FILE, done_mask, ...
    DVmatrix_sweep, TOFmatrix_sweep, TOFproxyMatrix_sweep, DVatMinTOFMatrix_sweep, ...
    winners_sweep, source_sweep);

% Early exit if everything is already done
if all(done_mask(:))
    fprintf('[sweep] All cells already done — skipping computation.\n');
else
    % ══════════════════════════════════════════════════════════════════════════
    %  PRE-FLIGHT: verify the base atlas cache exists before doing any work
    % ══════════════════════════════════════════════════════════════════════════
    [cache_ok, cache_report] = atlas_check_cache_exists(families, cfg); %#ok<ASGLU>
    if ~cache_ok
        error(['run_overlap_dv_tmax_sweep: one or more base-atlas cache files ' ...
               'are missing (see list above). Either run the atlas build for ' ...
               'those families first, or set cfg.cache.rebuild appropriately.']);
    end

    % ══════════════════════════════════════════════════════════════════════════
    %  LOAD BASE ATLASES  (Tmax = 3*pi/2, DV_cap = 0.3)  — once for the whole sweep
    % ══════════════════════════════════════════════════════════════════════════
    fprintf('[sweep] Loading base atlases (Tmax=3π/2, DV_cap=0.3)...\n');
    grid3_base = atlas_grid_make(cfg);
    Sall_base  = cell(N, 1);
    for i = 1:N
        fprintf('[sweep]   family %d/%d: %s\n', i, N, families{i});
        % Cache-only load (current scheme, then legacy_rs3 fallback) — the
        % pre-flight check above already guarantees a hit under one of the
        % two, so this never falls through to a from-scratch rebuild.
        [Sall_base{i}, scheme] = atlas_load_cached_compat(families{i}, cfg);
        if ~strcmp(scheme, 'current')
            fprintf('[sweep]     (loaded via %s fingerprint scheme)\n', scheme);
        end
    end
    fprintf('[sweep] Base atlases loaded.\n\n');

    % ══════════════════════════════════════════════════════════════════════════
    %  PARALLEL POOL SETUP
    % ══════════════════════════════════════════════════════════════════════════
    if N_WORKERS > 0
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local', N_WORKERS);
            fprintf('[sweep] Started parpool with %d workers.\n', N_WORKERS);
        else
            fprintf('[sweep] Using existing parpool (%d workers).\n', pool.NumWorkers);
        end
    else
        fprintf('[sweep] N_WORKERS = 0 → fully serial pair loop.\n');
    end

    % ══════════════════════════════════════════════════════════════════════════
    %  MAIN SWEEP LOOP  — serial over (DV_cap, Tmax) combinations
    %  NOTE: reversed order (large→small) so workers see peak-size data on
    %  cell 1 and their heap never grows beyond that for the rest of the sweep.
    % ══════════════════════════════════════════════════════════════════════════
    for di = nDV:-1:1
        for dj = nTmax:-1:1

            if done_mask(di, dj)
                fprintf('[sweep] (%d,%d) DV=%.3f  Tmax=%s — already done, skipping.\n', ...
                    di, dj, DV_cap_list(di), local_tmax_str(Tmax_list(dj)));
                continue;
            end

            dv_cap = DV_cap_list(di);
            tmax   = Tmax_list(dj);
            cellNo = (di - 1) * nTmax + dj;

            fprintf('\n[sweep] ══ cell %d/%d  |  DV_cap=%.3f  Tmax=%s ══\n', ...
                cellNo, nDV * nTmax, dv_cap, local_tmax_str(tmax));

            tCell = tic;

            % ── Build cfg for this combination ────────────────────────────────
            cfg_sub                = cfg;
            cfg_sub.propag.Tmax    = tmax;
            cfg_sub.fan.DV_cap_nd  = dv_cap;
            % Note: derived atlases are NOT saved to disk (no atlas_cache_save call).
            % SWEEP_CACHE_DIR is set here for completeness only.
            cfg_sub.cache.dir      = SWEEP_CACHE_DIR;

            % ── Derive in-memory subset atlases ───────────────────────────────
            fprintf('[sweep] Deriving subset atlases...\n');
            Sall_sub = cell(N, 1);
            if N_WORKERS > 0
                parfor i = 1:N
                    Sall_sub{i} = atlas_derive_subset(Sall_base{i}, cfg_sub);
                end
            else
                for i = 1:N
                    Sall_sub{i} = atlas_derive_subset(Sall_base{i}, cfg_sub);
                end
            end

            % ── Build compact per-family voxel footprints ─────────────────────
            % Each footprint is ~5-25 MB (unique voxel IDs + per-voxel dv_min/t_min)
            % vs. full atlas ~0.5-2 GB.  Workers in the pair parfor receive footprints
            % instead of full row structs → drastically reduces worker heap growth.
            fprintf('[sweep] Building voxel footprints...\n');
            Fall_sub = cell(N, 1);
            if N_WORKERS > 0
                parfor i = 1:N
                    Fall_sub{i} = overlap_proxy_footprint( ...
                        Sall_sub{i}, grid3_base, VU_mps, TU_days);
                end
            else
                for i = 1:N
                    Fall_sub{i} = overlap_proxy_footprint( ...
                        Sall_sub{i}, grid3_base, VU_mps, TU_days);
                end
            end
            clear Sall_sub   % row data no longer needed — footprints hold all pair needs

            % ── Pre-extract per-pair footprint pointers (required for parfor) ──
            FA_arr = cell(nPairs, 1);
            FB_arr = cell(nPairs, 1);
            for p = 1:nPairs
                FA_arr{p} = Fall_sub{pairI(p)};
                FB_arr{p} = Fall_sub{pairJ(p)};
            end
            clear Fall_sub   % individual footprints held by FA/FB_arr

            % ── Run pair loop ─────────────────────────────────────────────────
            % Each pair yields TWO independent winners over the same candidate
            % voxels: the min-DV-proxy voxel (unchanged from before) and the
            % min-TOF-proxy voxel (new).  See overlap_proxy_pair.m.
            pair_minDV      = nan(nPairs, 1);
            pair_DVlb       = nan(nPairs, 1);
            pair_DVpatch    = nan(nPairs, 1);
            pair_TOF        = nan(nPairs, 1);   % TOF at the min-DV winner voxel
            pair_voxelId    = nan(nPairs, 1);
            pair_minTOF     = nan(nPairs, 1);   % min-TOF-proxy winner
            pair_DVatMinTOF = nan(nPairs, 1);   % DV at the min-TOF winner voxel
            pair_voxelIdTOF = nan(nPairs, 1);

            if N_WORKERS > 0
                % parfor over pairs — FA/FB_arr sliced, grid3_base/cfg_sub/VU_mps broadcast
                parfor p = 1:nPairs
                    [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
                        pair_TOF(p), pair_voxelId(p), ...
                        pair_minTOF(p), pair_DVatMinTOF(p), pair_voxelIdTOF(p)] = ...
                        overlap_proxy_pair(FA_arr{p}, FB_arr{p}, grid3_base, cfg_sub, VU_mps);
                end
            else
                for p = 1:nPairs
                    fprintf('[sweep]   pair %d/%d: %-30s → %s\n', ...
                        p, nPairs, families{pairI(p)}, families{pairJ(p)});
                    [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
                        pair_TOF(p), pair_voxelId(p), ...
                        pair_minTOF(p), pair_DVatMinTOF(p), pair_voxelIdTOF(p)] = ...
                        overlap_proxy_pair(FA_arr{p}, FB_arr{p}, grid3_base, cfg_sub, VU_mps);
                end
            end

            clear FA_arr FB_arr   % free pair-footprint memory

            % ── Assemble N×N matrices ─────────────────────────────────────────
            minDVproxyMat  = nan(N, N);
            TOFmat         = nan(N, N);
            minTOFproxyMat = nan(N, N);
            DVatMinTOFmat  = nan(N, N);
            for p = 1:nPairs
                i = pairI(p);  j = pairJ(p);
                minDVproxyMat(i, j)  = pair_minDV(p);
                minDVproxyMat(j, i)  = pair_minDV(p);
                TOFmat(i, j)         = pair_TOF(p);
                TOFmat(j, i)         = pair_TOF(p);
                minTOFproxyMat(i, j) = pair_minTOF(p);
                minTOFproxyMat(j, i) = pair_minTOF(p);
                DVatMinTOFmat(i, j)  = pair_DVatMinTOF(p);
                DVatMinTOFmat(j, i)  = pair_DVatMinTOF(p);
            end

            % ── Build winners table (mirrors pair_winners_top1.csv) ───────────
            pair_famA = families(pairI);   % cell column
            pair_famB = families(pairJ);
            T = table(pair_famA(:), pair_famB(:), ...
                pair_minDV, pair_DVlb, pair_DVpatch, pair_TOF, pair_voxelId, ...
                pair_minTOF, pair_DVatMinTOF, pair_voxelIdTOF, ...
                'VariableNames', { ...
                    'FamilyA', 'FamilyB', ...
                    'minDVproxy_mps', 'DVlb_mps', 'DVpatch_ub_mps', ...
                    'TOFatMinDV_days', 'VoxelId_DV', ...
                    'minTOFproxy_days', 'DVatMinTOF_mps', 'VoxelId_TOF'});

            % ── Store in sweep arrays ─────────────────────────────────────────
            DVmatrix_sweep{di, dj}         = minDVproxyMat;
            TOFmatrix_sweep{di, dj}        = TOFmat;
            TOFproxyMatrix_sweep{di, dj}   = minTOFproxyMat;
            DVatMinTOFMatrix_sweep{di, dj} = DVatMinTOFmat;
            winners_sweep{di, dj}   = T;
            source_sweep{di, dj}    = 'computed';
            done_mask(di, dj)       = true;

            % ── Checkpoint immediately after cell completes ───────────────────
            local_save_checkpoint(CHECKPOINT_FILE, done_mask, ...
                DVmatrix_sweep, TOFmatrix_sweep, TOFproxyMatrix_sweep, DVatMinTOFMatrix_sweep, ...
                winners_sweep, source_sweep);

            fprintf('[sweep] Cell (%d,%d) done in %.1f s — checkpoint saved.\n', ...
                di, dj, toc(tCell));

        end  % dj
    end  % di

    clear Sall_base grid3_base   % base atlases no longer needed

end  % if ~all(done_mask)

% ══════════════════════════════════════════════════════════════════════════════
%  SAVE FINAL .MAT
% ══════════════════════════════════════════════════════════════════════════════
outMat = fullfile(OUTPUT_DIR, 'sweep_DVmatrix_results.mat');
save(outMat, ...
    'DV_cap_list', 'Tmax_list', 'Tmax_labels', 'families', ...
    'DVmatrix_sweep', 'TOFmatrix_sweep', ...
    'TOFproxyMatrix_sweep', 'DVatMinTOFMatrix_sweep', ...
    'winners_sweep', 'source_sweep', ...
    '-v7.3');
fprintf('\n[sweep] Final .mat saved:\n  %s\n', outMat);

% ══════════════════════════════════════════════════════════════════════════════
%  WRITE EXCEL WORKBOOKS (OPTIONAL)
% ══════════════════════════════════════════════════════════════════════════════
if WRITE_EXCEL
    fprintf('[sweep] Writing Excel files...\n');

    local_write_matrix_excel(OUTPUT_DIR, 'sweep_DVmatrix.xlsx', ...
        families, DV_cap_list, Tmax_list, Tmax_labels, DVmatrix_sweep);

    local_write_matrix_excel(OUTPUT_DIR, 'sweep_TOFmatrix.xlsx', ...
        families, DV_cap_list, Tmax_list, Tmax_labels, TOFmatrix_sweep);

    local_write_winners_excel(OUTPUT_DIR, 'sweep_winners.xlsx', ...
        DV_cap_list, Tmax_list, Tmax_labels, winners_sweep);
else
    fprintf('[sweep] Skipping Excel output (WRITE_EXCEL=false).\n');
end

% ══════════════════════════════════════════════════════════════════════════════
%  WRITE FLAT SUMMARY (OPTIONAL)
% ══════════════════════════════════════════════════════════════════════════════
if WRITE_FLAT_CSV
    fprintf('[sweep] Writing flat winners CSV...\n');
    Tall   = local_flatten_winners(DV_cap_list, Tmax_list, winners_sweep);
    outCsv = fullfile(OUTPUT_DIR, 'sweep_winners_flat.csv');
    writetable(Tall, outCsv);
    fprintf('[sweep]   Written: %s\n', outCsv);
end

fprintf('\n[sweep] ══════════ SWEEP COMPLETE ══════════\n');
fprintf('  Output dir : %s\n', OUTPUT_DIR);
fprintf('  Final .mat : %s\n', outMat);
if WRITE_EXCEL
    fprintf('  DVmatrix   : sweep_DVmatrix.xlsx\n');
    fprintf('  TOFmatrix  : sweep_TOFmatrix.xlsx\n');
    fprintf('  Winners    : sweep_winners.xlsx\n');
end
if WRITE_FLAT_CSV
    fprintf('  Flat CSV   : sweep_winners_flat.csv\n');
end

% ══════════════════════════════════════════════════════════════════════════════

%  LOCAL FUNCTIONS
%
%  NOTE: the per-family footprint builder and per-pair winner search used to
%  live here as local_compute_footprint / local_run_pair. They now live in
%  src/ as overlap_proxy_footprint.m / overlap_proxy_pair.m (reused by the
%  single-snapshot preview script run_overlap_mintof_maxbudget_preview.m).
% ══════════════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────────────
function local_save_checkpoint(fpath, done_mask, ...
        DVmatrix_sweep, TOFmatrix_sweep, TOFproxyMatrix_sweep, DVatMinTOFMatrix_sweep, ...
        winners_sweep, source_sweep)
%LOCAL_SAVE_CHECKPOINT  Atomically save sweep progress.
tmp = [fpath '.tmp'];
save(tmp, ...
    'done_mask', 'DVmatrix_sweep', 'TOFmatrix_sweep', ...
    'TOFproxyMatrix_sweep', 'DVatMinTOFMatrix_sweep', ...
    'winners_sweep', 'source_sweep', ...
    '-v7.3');
if isfile(fpath), delete(fpath); end
movefile(tmp, fpath);
end

% ─────────────────────────────────────────────────────────────────────────────
function existing = local_scan_results(resultsRoot, families, ...
        DV_cap_list, Tmax_list, cfg_ref)
%LOCAL_SCAN_RESULTS  Walk atlas_results, find batch_summary_workspace.mat files
% whose config matches (grid, tolerances, log) and whose (Tmax, DV_cap) falls
% on a sweep grid point.  Returns a nDV×nTmax cell: empty or data struct.
nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);
existing = cell(nDV, nTmax);

if ~isfolder(resultsRoot)
    fprintf('[sweep]   atlas_results not found — skipping scan.\n');
    return;
end

hits = dir(fullfile(resultsRoot, '**', 'batch_summary_workspace.mat'));
if isempty(hits)
    fprintf('[sweep]   No batch_summary_workspace.mat found in atlas_results.\n');
    return;
end

fprintf('[sweep]   Found %d workspace file(s) to inspect.\n', numel(hits));

for k = 1:numel(hits)
    fpath = fullfile(hits(k).folder, hits(k).name);
    try
        ws = load(fpath, 'cfg', 'minDVproxyMat', 'T', 'families');
    catch
        continue;
    end

    % Must have all required fields
    if ~all(isfield(ws, {'cfg', 'minDVproxyMat', 'T', 'families'}))
        continue;
    end

    % Family list must match exactly (same names, same order)
    if ~isequal(cellstr(ws.families(:)), cellstr(families(:)))
        continue;
    end

    % Non-sweep config fields must match (grid, tolerances, log, seed, fan res)
    if ~local_cfg_match(ws.cfg, cfg_ref)
        continue;
    end

    % Identify which (di, dj) cell this workspace corresponds to
    ws_tmax = ws.cfg.propag.Tmax;
    ws_dv   = ws.cfg.fan.DV_cap_nd;

    for di = 1:nDV
        for dj = 1:nTmax
            if abs(DV_cap_list(di) - ws_dv)   > 1e-9, continue; end
            if abs(Tmax_list(dj)   - ws_tmax) > 1e-9, continue; end
            if ~isempty(existing{di, dj}),             continue; end  % first match wins

            % Reconstruct symmetric TOF matrix from the winners table
            N      = numel(families);
            TOFmat = nan(N, N);
            if istable(ws.T) && height(ws.T) > 0
                famA_col = cellstr(ws.T.FamilyA);
                famB_col = cellstr(ws.T.FamilyB);
                for r = 1:height(ws.T)
                    ia = find(strcmp(families, famA_col{r}), 1);
                    ib = find(strcmp(families, famB_col{r}), 1);
                    if ~isempty(ia) && ~isempty(ib)
                        TOFmat(ia, ib) = ws.T.EstimatedTOF_days(r);
                        TOFmat(ib, ia) = ws.T.EstimatedTOF_days(r);
                    end
                end
            end

            existing{di, dj} = struct( ...
                'minDVproxyMat', ws.minDVproxyMat, ...
                'TOFmat',        TOFmat, ...
                'T',             ws.T);
        end
    end
end
end

% ─────────────────────────────────────────────────────────────────────────────
function ok = local_cfg_match(cfg_a, cfg_ref)
%LOCAL_CFG_MATCH  Return true when all non-sweep config fields agree.
% Intentionally excludes propag.Tmax and fan.DV_cap_nd (the sweep axes).
tol = 1e-9;
ok  = false;

% Each row: {struct_field, sub_field, reference_value}
checks = { ...
    'grid',   'dx',                cfg_ref.grid.dx; ...
    'grid',   'dy',                cfg_ref.grid.dy; ...
    'grid',   'dtheta',            cfg_ref.grid.dtheta; ...
    'seed',   'ds_seed',           cfg_ref.seed.ds_seed; ...
    'fan',    'dtheta_fan',        cfg_ref.fan.dtheta_fan; ...
    'propag', 'absTol',            cfg_ref.propag.absTol; ...
    'propag', 'relTol',            cfg_ref.propag.relTol; ...
    'propag', 'v2tol',             cfg_ref.propag.v2tol; ...
    'log',    'step_len_factor',   cfg_ref.log.step_len_factor; ...
    'log',    'maxstep_factor',    cfg_ref.log.maxstep_factor; ...
};

for i = 1:size(checks, 1)
    f1 = checks{i, 1};
    f2 = checks{i, 2};
    rv = checks{i, 3};
    if ~isfield(cfg_a, f1) || ~isfield(cfg_a.(f1), f2)
        return;
    end
    if abs(cfg_a.(f1).(f2) - rv) > tol
        return;
    end
end
ok = true;
end

% ─────────────────────────────────────────────────────────────────────────────
% ─────────────────────────────────────────────────────────────────────────────
function Tall = local_flatten_winners(DV_cap_list, Tmax_list, winners_sweep)
%LOCAL_FLATTEN_WINNERS  Combine per-cell winners tables into one flat table.
nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);

Ts = cell(nDV * nTmax, 1);
k  = 0;

for di = 1:nDV
    for dj = 1:nTmax
        T = winners_sweep{di, dj};
        if isempty(T), continue; end

        T.DV_cap_nd = repmat(DV_cap_list(di), height(T), 1);
        T.Tmax      = repmat(Tmax_list(dj),   height(T), 1);

        % Put sweep axes first (if movevars exists)
        if exist('movevars', 'file') == 2
            T = movevars(T, {'DV_cap_nd','Tmax'}, 'Before', 1);
        end

        k = k + 1;
        Ts{k} = T;
    end
end

Ts = Ts(1:k);
if isempty(Ts)
    Tall = table();
else
    Tall = vertcat(Ts{:});
end
end

function local_write_matrix_excel(outdir, fname, families, ...
        DV_cap_list, Tmax_list, Tmax_labels, data_sweep)
%LOCAL_WRITE_MATRIX_EXCEL  Write N×N matrices to an Excel workbook.
% One sheet per (DV_cap, Tmax) combination; family names as row/column headers.
fpath = fullfile(outdir, fname);
if isfile(fpath), delete(fpath); end

N     = numel(families);
nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);

headerRow = [{'Family'}, families(:)'];

for di = 1:nDV
    for dj = 1:nTmax
        sheetName = local_sheet_name(DV_cap_list(di), Tmax_labels{dj});
        mat = data_sweep{di, dj};
        if isempty(mat)
            cellData = [headerRow; repmat({'(no data)'}, N, N + 1)];
        else
            cellData = [headerRow; [families(:), num2cell(mat)]];
        end
        writecell(cellData, fpath, 'Sheet', sheetName);
    end
end

fprintf('[sweep]   Written: %s\n', fname);
end

% ─────────────────────────────────────────────────────────────────────────────
function local_write_winners_excel(outdir, fname, ...
        DV_cap_list, Tmax_list, Tmax_labels, winners_sweep)
%LOCAL_WRITE_WINNERS_EXCEL  Write pair-winner tables to an Excel workbook.
% One sheet per (DV_cap, Tmax) combination.
fpath = fullfile(outdir, fname);
if isfile(fpath), delete(fpath); end

nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);

for di = 1:nDV
    for dj = 1:nTmax
        sheetName = local_sheet_name(DV_cap_list(di), Tmax_labels{dj});
        T = winners_sweep{di, dj};
        if isempty(T)
            writecell({'(no data)'}, fpath, 'Sheet', sheetName);
        else
            writetable(T, fpath, 'Sheet', sheetName);
        end
    end
end

fprintf('[sweep]   Written: %s\n', fname);
end

% ─────────────────────────────────────────────────────────────────────────────
function s = local_sheet_name(dv_cap, tmax_label)
%LOCAL_SHEET_NAME  Build an Excel sheet name ≤31 chars.
% Format: "DV0p025_Tpi4"
dv_str = strrep(sprintf('%.3f', dv_cap), '.', 'p');
s = sprintf('DV%s_%s', dv_str, tmax_label);
% Truncate to 31 chars just in case
if numel(s) > 31, s = s(1:31); end
end

% ─────────────────────────────────────────────────────────────────────────────
function s = local_tmax_str(tmax)
%LOCAL_TMAX_STR  Human-readable label for a Tmax value.
known_val = [pi/4,   pi/3,   pi/2,   2*pi/3,  3*pi/4,  pi  ];
known_str = {'pi/4', 'pi/3', 'pi/2', '2*pi/3','3*pi/4','pi' };
[~, idx] = min(abs(tmax - known_val));
if abs(tmax - known_val(idx)) < 1e-9
    s = known_str{idx};
else
    s = sprintf('%.5g', tmax);
end
end
