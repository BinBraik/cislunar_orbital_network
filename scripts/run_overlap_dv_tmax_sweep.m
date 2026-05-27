%% RUN_RS4_DV_TMAX_SWEEP
% Sweep DVmatrix (and TOFmatrix) over all (DV_cap_nd, Tmax) combinations.
%
% Strategy:
%   1. Load base atlases once  (Tmax = pi, DV_cap = 0.2).
%   2. For each (DV_cap, Tmax) cell:
%        a. Check atlas_results for a finished run whose config matches
%           → extract minDVproxyMat + TOFmatrix directly, skip recompute.
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
%   sweep_DVmatrix_results.mat  — authoritative data store
%   sweep_DVmatrix.xlsx         — one sheet per combination, N×N DVproxy matrix
%   sweep_TOFmatrix.xlsx        — one sheet per combination, N×N mean-TOF matrix
%   sweep_winners.xlsx          — one sheet per combination, pair-winner table

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
DV_cap_list = linspace(0.025, 0.200, 20)';   % 20 values  (0.025 → 0.200)
Tmax_list   = linspace(pi/4,  pi,    20)';   % 20 values  (pi/4  → pi)

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

cfg.cache.enable  = true;
cfg.cache.dir     = fullfile(repoRoot, 'atlas_cache');
cfg.cache.rebuild = false;

% Suppress all figure/file output from sub-functions
cfg.io.save_figs   = false;
cfg.io.save_fig    = false;
cfg.io.fig_visible = 'off';

cfg.plot.rs4.overlap_xy   = false;
cfg.plot.rs4.overlap_xyz  = false;
cfg.plot.rs4.combo_xy     = false;
cfg.plot.rs4.combo_xyz    = false;
cfg.plot.rs4.bounds_lb    = false;
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = false;

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
done_mask       = false(nDV, nTmax);
DVmatrix_sweep  = cell(nDV, nTmax);
TOFmatrix_sweep = cell(nDV, nTmax);
winners_sweep   = cell(nDV, nTmax);
source_sweep    = cell(nDV, nTmax);

if isfile(CHECKPOINT_FILE)
    fprintf('[sweep] Checkpoint found — loading...\n');
    try
        ck = load(CHECKPOINT_FILE, ...
            'done_mask', 'DVmatrix_sweep', 'TOFmatrix_sweep', 'winners_sweep', 'source_sweep');
        if isequal(size(ck.done_mask), [nDV, nTmax])
            done_mask       = ck.done_mask;
            DVmatrix_sweep  = ck.DVmatrix_sweep;
            TOFmatrix_sweep = ck.TOFmatrix_sweep;
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
            winners_sweep{di, dj}   = e.T;
            source_sweep{di, dj}    = 'loaded_from_results';
            done_mask(di, dj)       = true;
            nFound = nFound + 1;
            fprintf('[sweep]   loaded (%d,%d): DV=%.3f  Tmax=%s\n', ...
                di, dj, DV_cap_list(di), local_tmax_str(Tmax_list(dj)));
        end
    end
end
fprintf('[sweep] %d cell(s) loaded from existing results. %d remain to compute.\n', ...
    nFound, sum(~done_mask(:)));

local_save_checkpoint(CHECKPOINT_FILE, done_mask, ...
    DVmatrix_sweep, TOFmatrix_sweep, winners_sweep, source_sweep);

% Early exit if everything is already done
if all(done_mask(:))
    fprintf('[sweep] All cells already done — skipping computation.\n');
else
    % ══════════════════════════════════════════════════════════════════════════
    %  LOAD BASE ATLASES  (Tmax = pi, DV_cap = 0.2)  — once for the whole sweep
    % ══════════════════════════════════════════════════════════════════════════
    fprintf('[sweep] Loading base atlases (Tmax=π, DV_cap=0.2)...\n');
    grid3_base = atlas_grid_make(cfg);
    Sall_base  = cell(N, 1);
    for i = 1:N
        fprintf('[sweep]   family %d/%d: %s\n', i, N, families{i});
        [Sall_base{i}, ~] = atlas_prepare_or_load(families{i}, cfg, grid3_base);
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
            % Each footprint is ~5-25 MB (unique voxel IDs + per-voxel dv_min/t_mean)
            % vs. full atlas ~0.5-2 GB.  Workers in the pair parfor receive footprints
            % instead of full row structs → drastically reduces worker heap growth.
            fprintf('[sweep] Building voxel footprints...\n');
            Fall_sub = cell(N, 1);
            if N_WORKERS > 0
                parfor i = 1:N
                    Fall_sub{i} = local_compute_footprint( ...
                        Sall_sub{i}, grid3_base, VU_mps, TU_days);
                end
            else
                for i = 1:N
                    Fall_sub{i} = local_compute_footprint( ...
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
            pair_minDV   = nan(nPairs, 1);
            pair_DVlb    = nan(nPairs, 1);
            pair_DVpatch = nan(nPairs, 1);
            pair_TOF     = nan(nPairs, 1);
            pair_voxelId = nan(nPairs, 1);

            if N_WORKERS > 0
                % parfor over pairs — FA/FB_arr sliced, grid3_base/cfg_sub/VU_mps broadcast
                parfor p = 1:nPairs
                    [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
                        pair_TOF(p), pair_voxelId(p)] = ...
                        local_run_pair(FA_arr{p}, FB_arr{p}, grid3_base, cfg_sub, VU_mps);
                end
            else
                for p = 1:nPairs
                    fprintf('[sweep]   pair %d/%d: %-30s → %s\n', ...
                        p, nPairs, families{pairI(p)}, families{pairJ(p)});
                    [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
                        pair_TOF(p), pair_voxelId(p)] = ...
                        local_run_pair(FA_arr{p}, FB_arr{p}, grid3_base, cfg_sub, VU_mps);
                end
            end

            clear FA_arr FB_arr   % free pair-footprint memory

            % ── Assemble N×N matrices ─────────────────────────────────────────
            minDVproxyMat = nan(N, N);
            TOFmat        = nan(N, N);
            for p = 1:nPairs
                i = pairI(p);  j = pairJ(p);
                minDVproxyMat(i, j) = pair_minDV(p);
                minDVproxyMat(j, i) = pair_minDV(p);
                TOFmat(i, j)        = pair_TOF(p);
                TOFmat(j, i)        = pair_TOF(p);
            end

            % ── Build winners table (mirrors pair_winners_top1.csv) ───────────
            pair_famA = families(pairI);   % cell column
            pair_famB = families(pairJ);
            T = table(pair_famA(:), pair_famB(:), ...
                pair_minDV, pair_DVlb, pair_DVpatch, pair_TOF, pair_voxelId, ...
                'VariableNames', { ...
                    'FamilyA', 'FamilyB', ...
                    'minDVproxy_mps', 'DVlb_mps', 'DVpatch_ub_mps', ...
                    'EstimatedTOF_days', 'VoxelId'});

            % ── Store in sweep arrays ─────────────────────────────────────────
            DVmatrix_sweep{di, dj}  = minDVproxyMat;
            TOFmatrix_sweep{di, dj} = TOFmat;
            winners_sweep{di, dj}   = T;
            source_sweep{di, dj}    = 'computed';
            done_mask(di, dj)       = true;

            % ── Checkpoint immediately after cell completes ───────────────────
            local_save_checkpoint(CHECKPOINT_FILE, done_mask, ...
                DVmatrix_sweep, TOFmatrix_sweep, winners_sweep, source_sweep);

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
    'DVmatrix_sweep', 'TOFmatrix_sweep', 'winners_sweep', 'source_sweep', ...
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
% ══════════════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────────────
function [minDV, dvlb, dvpatch, tof, voxelId] = ...
        local_run_pair(FA, FB, grid3, cfg, VU_mps)
%LOCAL_RUN_PAIR  Compute DV proxy and TOF for one pair using pre-computed footprints.
%
% FA / FB are compact structs produced by local_compute_footprint; they contain
% unique sorted voxel-ID vectors + per-voxel dv_min and t_mean for FRS and BRS.
% No full row data is needed — workers only receive ~5-25 MB instead of ~1-3 GB.
%
% Results are numerically identical to the overlap_pair +
% overlap_extract_voxel_info + inline-proxy pipeline.
minDV = NaN;  dvlb = NaN;  dvpatch = NaN;  tof = NaN;  voxelId = NaN;
try
    % ── 1. Intersect FRS(A) with BRS(B) (both already sorted unique) ──────
    idsO = intersect(FA.uid_frs, FB.uid_brs);
    if isempty(idsO), return; end

    % ── 2. Unpack voxel grid indices ───────────────────────────────────────
    Ny = numel(grid3.y_centers);
    Nx = numel(grid3.x_centers);
    Nt = numel(grid3.th_centers);
    [iy, ix, ~] = ind2sub([Ny, Nx, Nt], idsO);

    % ── 3. Keep + primary-buffer filter (mirrors overlap_pair) ─────────
    bufFrac = 0.05;
    if isfield(cfg,'overlap') && isfield(cfg.overlap,'primary_buffer_frac') ...
            && ~isempty(cfg.overlap.primary_buffer_frac)
        bufFrac = cfg.overlap.primary_buffer_frac;
    end
    if ~(isfield(cfg,'sys') && isfield(cfg.sys,'RE_nd') && isfield(cfg.sys,'RM_nd'))
        error('cfg.sys.RE_nd and cfg.sys.RM_nd required for primary buffer filter.');
    end
    RE = cfg.sys.RE_nd;
    RM = cfg.sys.RM_nd;
    mu = FA.mu;

    % Keep mask — all families share grid3_base so keepA = keepB = grid3.Keep
    if isfield(grid3,'Keep') && ~isempty(grid3.Keep)
        keepXY = logical(grid3.Keep);
        if ~isequal(size(keepXY), [Ny, Nx]), keepXY = keepXY.'; end
        okKeep = keepXY(sub2ind([Ny, Nx], iy, ix));
    else
        okKeep = true(numel(idsO), 1);
    end

    x = grid3.x_centers(ix);
    y = grid3.y_centers(iy);
    okEarth = hypot(x + mu, y)     > (1 + bufFrac) * RE;
    okMoon  = hypot(x - (1-mu), y) > (1 + bufFrac) * RM;
    ok = okKeep(:) & okEarth(:) & okMoon(:);

    idsO = idsO(ok);
    if isempty(idsO), return; end
    ix = ix(ok);  iy = iy(ok);

    % ── 4. Look up pre-computed per-voxel DV and TOF ──────────────────────
    % FA.uid_frs and FB.uid_brs are sorted → ismember is fast
    [~, locA] = ismember(idsO, FA.uid_frs);
    [~, locB] = ismember(idsO, FB.uid_brs);
    dv_min_A = FA.dv_min_frs(locA);
    dv_min_B = FB.dv_min_brs(locB);
    t_mean_A = FA.t_mean_frs(locA);
    t_mean_B = FB.t_mean_brs(locB);

    % ── 5. DV proxy (identical formula to original local_run_pair) ─────────
    x_ok = grid3.x_centers(ix);
    y_ok = grid3.y_centers(iy);
    CJstar = min(FA.CJ, FB.CJ);
    pot = cr3bp_potential(x_ok(:), y_ok(:), mu);
    v_box = sqrt(max(2 * pot.U - CJstar, 0));
    dv_patch_vec = 2 * v_box .* sin(abs(grid3.dtheta) / 2) * VU_mps;
    dv_lb_vec    = dv_min_A(:) + dv_min_B(:);
    dv_proxy     = dv_lb_vec + dv_patch_vec;

    valid = isfinite(dv_proxy);
    if ~any(valid), return; end

    idxValid  = find(valid);
    [~, iLoc] = min(dv_proxy(idxValid));
    iWin      = idxValid(iLoc);

    minDV   = dv_proxy(iWin);
    dvlb    = dv_lb_vec(iWin);
    dvpatch = dv_patch_vec(iWin);
    voxelId = idsO(iWin);
    tof     = t_mean_A(iWin) + t_mean_B(iWin);
catch ME
    warning('[sweep:pair] %s', ME.message);
end
end

% ─────────────────────────────────────────────────────────────────────────────
function F = local_compute_footprint(S, grid3, VU_mps, TU_days)
%LOCAL_COMPUTE_FOOTPRINT  Build compact per-voxel summary for one atlas family.
%
% Computes, for FRS and BRS separately:
%   uid       — sorted unique voxel IDs  (linear index into [Ny,Nx,Nt])
%   dv_min    — min dv_turn (m/s) over all rows in that voxel
%   t_mean    — mean |TOF| (days) over all rows in that voxel
%
% FRS  = direct rows (FRS_upper + FRS_lower).
% BRS  = R(FRS): mirror of FRS_upper + mirror of FRS_lower.
%        Mirror: iy → Ny-iy+1,  it → it_lut(it)  (same as overlap_pair).
%
% DV reuse: since U(x,y)=U(x,-y) in CR3BP, the DV computed at a seed and at
% its y-mirror are identical.  BRS voxels therefore reuse the DV values from
% the corresponding FRS rows — no extra potential evaluation needed.

Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

% ── Theta mirror LUT (identical formula to overlap_pair) ───────────────
thm = wrap_to_pi(pi - grid3.th_centers(:));
lut = discretize(thm, grid3.th_edges);
lut(isnan(lut)) = 0;
it_lut = uint16(lut);

% ── Pre-build delta-angle lookup matrix (vectorised, avoids cell loop) ──────
dlists = S.Step4.delta_lists;
Ns     = numel(dlists);
max_h  = max(1, max(cellfun(@numel, dlists)));
delta_mat = zeros(Ns, max_h);
for s = 1:Ns
    v = double(dlists{s});
    delta_mat(s, 1:numel(v)) = v;
end

% ── Pre-compute v0 per unique seed (avoids per-row potential evaluation) ────
% cr3bp_potential on ~1000 seeds is O(1000x) faster than on
% millions of rows.  Per-row v0 is then a cheap index lookup: v0(iSeed).
pot_u = cr3bp_potential(S.SeedsUpper(:,1), S.SeedsUpper(:,2), S.mu);
v0_upper = sqrt(max(2*pot_u.U(:) - S.CJ, 0));   % [Nseeds_upper, 1]

if isfield(S,'SeedsLower') && ~isempty(S.SeedsLower)
    pot_l = cr3bp_potential(S.SeedsLower(:,1), S.SeedsLower(:,2), S.mu);
    v0_lower = sqrt(max(2*pot_l.U(:) - S.CJ, 0));
else
    v0_lower = v0_upper;
end

% ── Process FRS_upper rows ─────────────────────────────────────────────────
nu = double(S.Step4.rows_FRS_upper.n);
if nu > 0
    [ids_u, dv_u, t_u, ix_u, iy_u, it_u] = local_fp_rows( ...
        S.Step4.rows_FRS_upper, nu, v0_upper, ...
        delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
else
    ids_u=zeros(0,1); dv_u=zeros(0,1); t_u=zeros(0,1);
    ix_u=zeros(0,1);  iy_u=zeros(0,1); it_u=zeros(0,1);
end

% ── Process FRS_lower rows ─────────────────────────────────────────────────
nl = double(S.Step4.rows_FRS_lower.n);
if nl > 0
    [ids_l, dv_l, t_l, ix_l, iy_l, it_l] = local_fp_rows( ...
        S.Step4.rows_FRS_lower, nl, v0_lower, ...
        delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
else
    ids_l=zeros(0,1); dv_l=zeros(0,1); t_l=zeros(0,1);
    ix_l=zeros(0,1);  iy_l=zeros(0,1); it_l=zeros(0,1);
end

% ── Aggregate FRS voxels ───────────────────────────────────────────────────
[F.uid_frs, F.dv_min_frs, F.t_mean_frs] = local_fp_agg( ...
    [ids_u; ids_l], [dv_u; dv_l], [t_u; t_l]);

% ── BRS: mirror FRS_upper ──────────────────────────────────────────────────
if ~isempty(ix_u)
    biy_u = Ny - iy_u + 1;
    bit_u = double(it_lut(it_u));
    ok_u  = bit_u > 0;
    ids_brs_u = sub2ind([Ny,Nx,Nt], biy_u(ok_u), ix_u(ok_u), ...
        max(1, min(Nt, bit_u(ok_u))));
    dv_bu = dv_u(ok_u);  t_bu = t_u(ok_u);
else
    ids_brs_u = zeros(0,1);  dv_bu = zeros(0,1);  t_bu = zeros(0,1);
end

% ── BRS: mirror FRS_lower ──────────────────────────────────────────────────
if ~isempty(ix_l)
    biy_l = Ny - iy_l + 1;
    bit_l = double(it_lut(it_l));
    ok_l  = bit_l > 0;
    ids_brs_l = sub2ind([Ny,Nx,Nt], biy_l(ok_l), ix_l(ok_l), ...
        max(1, min(Nt, bit_l(ok_l))));
    dv_bl = dv_l(ok_l);  t_bl = t_l(ok_l);
else
    ids_brs_l = zeros(0,1);  dv_bl = zeros(0,1);  t_bl = zeros(0,1);
end

% ── Aggregate BRS voxels ───────────────────────────────────────────────────
[F.uid_brs, F.dv_min_brs, F.t_mean_brs] = local_fp_agg( ...
    [ids_brs_u; ids_brs_l], [dv_bu; dv_bl], [t_bu; t_bl]);

F.CJ   = S.CJ;
F.mu   = S.mu;
F.name = S.name;
end

% ─────────────────────────────────────────────────────────────────────────────
function [ids, dv_mps, t_days, ix_out, iy_out, it_out] = local_fp_rows( ...
        rows, n, v0_per_seed, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days)
%LOCAL_FP_ROWS  Extract voxel IDs, DV (m/s), and |TOF| (days) for n packed rows.
% v0_per_seed  [Nseeds,1] — pre-computed sqrt(max(2U-CJ,0)) per seed position.
%   Caller computes this once with cr3bp_potential on the seeds
%   matrix (~hundreds of evals) so this function avoids per-row pot evaluation.
ix_out = double(rows.ix(1:n));
iy_out = double(rows.iy(1:n));
it_out = double(rows.it(1:n));
ids    = sub2ind([Ny, Nx, Nt], iy_out, ix_out, it_out);

iSeed = double(rows.iSeed(1:n));
iHead = double(rows.iHead(1:n));
t_nd  = double(rows.t(1:n));

% Vectorised delta-angle lookup
lin   = sub2ind([Ns, max_h], iSeed, iHead);
delta = delta_mat(lin);

% Per-row v0 via seed lookup (O(n) index, no potential evaluation)
v0     = v0_per_seed(iSeed);
dv_mps = 2 * v0(:) .* sin(abs(delta(:)) / 2) * VU_mps;
t_days = abs(t_nd(:)) * TU_days;
end

% ─────────────────────────────────────────────────────────────────────────────
function [uid, dv_min, t_mean] = local_fp_agg(ids, dv, t)
%LOCAL_FP_AGG  Aggregate per-voxel min DV and mean TOF from raw row arrays.
if isempty(ids)
    uid = zeros(0,1);  dv_min = zeros(0,1);  t_mean = zeros(0,1);
    return;
end
[uid, ~, ic] = unique(ids(:));
dv_min = accumarray(ic, dv(:), [], @min);          % @min has a native built-in path
t_mean = accumarray(ic, t(:)) ./ accumarray(ic, ones(numel(ic), 1)); % @mean does not → use sum/count
end

% ─────────────────────────────────────────────────────────────────────────────
function local_save_checkpoint(fpath, done_mask, ...
        DVmatrix_sweep, TOFmatrix_sweep, winners_sweep, source_sweep)
%LOCAL_SAVE_CHECKPOINT  Atomically save sweep progress.
tmp = [fpath '.tmp'];
save(tmp, ...
    'done_mask', 'DVmatrix_sweep', 'TOFmatrix_sweep', 'winners_sweep', 'source_sweep', ...
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
