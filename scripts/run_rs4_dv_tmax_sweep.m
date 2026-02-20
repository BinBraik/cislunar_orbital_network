%% RUN_RS4_DV_TMAX_SWEEP
% Sweep DVmatrix (and TOFmatrix) over all (DV_cap_nd, Tmax) combinations.
%
% Strategy:
%   1. Load base atlases once  (Tmax = pi, DV_cap = 0.2).
%   2. For each (DV_cap, Tmax) cell:
%        a. Check rs3_results for a finished run whose config matches
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
CHECKPOINT_FILE = fullfile(repoRoot, 'rs3_sweep_checkpoint.mat');
OUTPUT_DIR      = fullfile(repoRoot, 'rs3_sweep_results');
SWEEP_CACHE_DIR = fullfile(repoRoot, 'rs3_cache_sweep');   % derived caches (if ever saved)

% ── Sweep grid ────────────────────────────────────────────────────────────────
DV_cap_list = (0.025 : 0.025 : 0.200)';   % 8 values
Tmax_list   = [pi/4; pi/3; pi/2; 2*pi/3; 3*pi/4; pi];   % 6 values

% Short labels for Excel sheet names  (≤31 chars; keep them compact)
Tmax_labels = {'Tpi4', 'Tpi3', 'Tpi2', 'T2pi3', 'T3pi4', 'Tpi'};

% ══════════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION  (must match your full cached atlases exactly)
% ══════════════════════════════════════════════════════════════════════════════
cfg_base = rs3_cfg_defaults();

cfg_base.families.list      = families;
cfg_base.families.test_only = false;

cfg_base.grid.dx     = 0.001;
cfg_base.grid.dy     = 0.001;
cfg_base.grid.dtheta = deg2rad(0.5);

cfg_base.seed.ds_seed = 0.01;

cfg_base.propag.Tmax    = pi;    % base — largest Tmax in the sweep
cfg_base.fan.DV_cap_nd  = 0.2;  % base — largest DV_cap in the sweep
cfg_base.fan.dtheta_fan = deg2rad(0.5);

cfg_base.propag.absTol = 1e-8;
cfg_base.propag.relTol = 1e-8;
cfg_base.propag.v2tol  = 1e-8;

cfg_base.log.step_len_factor = 0.75;
cfg_base.log.maxstep_factor  = 0.75;

cfg_base.cache.enable  = true;
cfg_base.cache.dir     = fullfile(repoRoot, 'rs3_cache');
cfg_base.cache.rebuild = false;

% Suppress all figure/file output from sub-functions
cfg_base.io.save_figs   = false;
cfg_base.io.save_fig    = false;
cfg_base.io.fig_visible = 'off';

cfg_base.plot.rs4.overlap_xy   = false;
cfg_base.plot.rs4.overlap_xyz  = false;
cfg_base.plot.rs4.combo_xy     = false;
cfg_base.plot.rs4.combo_xyz    = false;
cfg_base.plot.rs4.bounds_lb    = false;
cfg_base.plot.rs4.bounds_ub    = false;
cfg_base.plot.rs4.bounds_proxy = false;

if exist('rs3_cfg_validate', 'file') == 2
    rs3_cfg_validate(cfg_base);
end

% ── Derived constants ─────────────────────────────────────────────────────────
VU_mps = cfg_base.units.VU_mps;
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
%  SCAN rs3_results FOR REUSABLE COMPLETED RUNS
% ══════════════════════════════════════════════════════════════════════════════
fprintf('[sweep] Scanning rs3_results for matching completed runs...\n');
existing = local_scan_results(cfg_base.io.out_root, families, ...
    DV_cap_list, Tmax_list, cfg_base);

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
    grid3_base = rs3_grid_make(cfg_base);
    Sall_base  = cell(N, 1);
    for i = 1:N
        fprintf('[sweep]   family %d/%d: %s\n', i, N, families{i});
        [Sall_base{i}, ~] = rs3_prepare_or_load_family(families{i}, cfg_base, grid3_base);
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
    % ══════════════════════════════════════════════════════════════════════════
    for di = 1:nDV
        for dj = 1:nTmax

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
            cfg_sub                = cfg_base;
            cfg_sub.propag.Tmax    = tmax;
            cfg_sub.fan.DV_cap_nd  = dv_cap;
            % Note: derived atlases are NOT saved to disk (no rs3_cache_save_family call).
            % SWEEP_CACHE_DIR is set here for completeness only.
            cfg_sub.cache.dir      = SWEEP_CACHE_DIR;

            % ── Derive in-memory subset atlases ───────────────────────────────
            fprintf('[sweep] Deriving subset atlases...\n');
            Sall_sub = cell(N, 1);
            for i = 1:N
                Sall_sub{i} = rs3_atlas_derive_subset(Sall_base{i}, cfg_sub);
            end

            % ── Pre-extract per-pair atlas pointers (required for parfor) ─────
            SA_arr = cell(nPairs, 1);
            SB_arr = cell(nPairs, 1);
            for p = 1:nPairs
                SA_arr{p} = Sall_sub{pairI(p)};
                SB_arr{p} = Sall_sub{pairJ(p)};
            end
            clear Sall_sub   % free memory — pairs hold only what they need

            % ── Run pair loop ─────────────────────────────────────────────────
            pair_minDV   = nan(nPairs, 1);
            pair_DVlb    = nan(nPairs, 1);
            pair_DVpatch = nan(nPairs, 1);
            pair_TOF     = nan(nPairs, 1);
            pair_voxelId = nan(nPairs, 1);

            if N_WORKERS > 0
                % parfor over pairs — SA_arr/SB_arr are sliced, cfg_sub/VU_mps broadcast
                parfor p = 1:nPairs
                    [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
                        pair_TOF(p), pair_voxelId(p)] = ...
                        local_run_pair(SA_arr{p}, SB_arr{p}, cfg_sub, VU_mps);
                end
            else
                for p = 1:nPairs
                    fprintf('[sweep]   pair %d/%d: %-30s → %s\n', ...
                        p, nPairs, families{pairI(p)}, families{pairJ(p)});
                    [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
                        pair_TOF(p), pair_voxelId(p)] = ...
                        local_run_pair(SA_arr{p}, SB_arr{p}, cfg_sub, VU_mps);
                end
            end

            clear SA_arr SB_arr   % free pair-atlas memory

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
%  WRITE EXCEL WORKBOOKS
% ══════════════════════════════════════════════════════════════════════════════
fprintf('[sweep] Writing Excel files...\n');

local_write_matrix_excel(OUTPUT_DIR, 'sweep_DVmatrix.xlsx', ...
    families, DV_cap_list, Tmax_list, Tmax_labels, DVmatrix_sweep);

local_write_matrix_excel(OUTPUT_DIR, 'sweep_TOFmatrix.xlsx', ...
    families, DV_cap_list, Tmax_list, Tmax_labels, TOFmatrix_sweep);

local_write_winners_excel(OUTPUT_DIR, 'sweep_winners.xlsx', ...
    DV_cap_list, Tmax_list, Tmax_labels, winners_sweep);

fprintf('\n[sweep] ══════════ SWEEP COMPLETE ══════════\n');
fprintf('  Output dir : %s\n', OUTPUT_DIR);
fprintf('  Final .mat : %s\n', outMat);
fprintf('  DVmatrix   : sweep_DVmatrix.xlsx\n');
fprintf('  TOFmatrix  : sweep_TOFmatrix.xlsx\n');
fprintf('  Winners    : sweep_winners.xlsx\n');

% ══════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────────────
function [minDV, dvlb, dvpatch, tof, voxelId] = ...
        local_run_pair(SA, SB, cfg, VU_mps)
%LOCAL_RUN_PAIR  Compute DV proxy and TOF for one pair — no files, no figures.
minDV = NaN;  dvlb = NaN;  dvpatch = NaN;  tof = NaN;  voxelId = NaN;
try
    O = rs4_overlap_pair(SA, SB, cfg);
    if isempty(O.ids), return; end

    V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg);
    if isempty(V.x), return; end

    % Inline DV-proxy calculation (mirrors rs4_overlap_visualize_bounds, no plots)
    CJstar = min(SA.CJ, SB.CJ);
    pot    = rs3_core_cr3bp_U_and_derivs(V.x(:), V.y(:), SA.mu);
    v_box  = sqrt(max(2 * pot.U - CJstar, 0));
    dv_patch_vec = 2 * v_box .* sin(abs(SA.grid3.dtheta) / 2) * VU_mps;
    dv_lb_vec    = V.dv_turn_mps_min_A + V.dv_turn_mps_min_B;
    dv_proxy     = dv_lb_vec + dv_patch_vec;

    valid = isfinite(dv_proxy);
    if ~any(valid), return; end

    idxValid   = find(valid);
    [~, iLoc]  = min(dv_proxy(idxValid));
    iWin       = idxValid(iLoc);

    minDV   = dv_proxy(iWin);
    dvlb    = dv_lb_vec(iWin);
    dvpatch = dv_patch_vec(iWin);

    if iWin >= 1 && iWin <= numel(V.ids)
        voxelId = V.ids(iWin);
        tof     = V.t_days_mean_A(iWin) + V.t_days_mean_B(iWin);
    end
catch ME
    warning('[sweep:pair] %s', ME.message);
end
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
%LOCAL_SCAN_RESULTS  Walk rs3_results, find batch_summary_workspace.mat files
% whose config matches (grid, tolerances, log) and whose (Tmax, DV_cap) falls
% on a sweep grid point.  Returns a nDV×nTmax cell: empty or data struct.
nDV   = numel(DV_cap_list);
nTmax = numel(Tmax_list);
existing = cell(nDV, nTmax);

if ~isfolder(resultsRoot)
    fprintf('[sweep]   rs3_results not found — skipping scan.\n');
    return;
end

hits = dir(fullfile(resultsRoot, '**', 'batch_summary_workspace.mat'));
if isempty(hits)
    fprintf('[sweep]   No batch_summary_workspace.mat found in rs3_results.\n');
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
