%% RUN_RS4_DC_SWEEP_TMAX
% Two-phase parametric sweep: footprint overlap (Phase 1) + differential
% correction (Phase 2) over all (DV_cap_nd, Tmax) combinations.
%
% Strategy:
%   1. Load base atlases once  (Tmax = pi, DV_cap = 0.2) with dense PO trace.
%   2. For each (DV_cap, Tmax) cell:
%        a. Derive in-memory subset atlases via rs3_atlas_derive_subset.
%        b. Phase 1: build extended footprints (compact, ~5-25 MB/family),
%           run pairwise intersection in parfor → one "ticket" per pair
%           (winner voxel ID + argmin row identifiers).
%        c. Phase 2: for each pair, re-integrate only the two winning arcs
%           and run rs4_diffcorr. Uses sliced parfor so each worker receives
%           only 2 subset atlases — no broadcast of all families.
%        d. Assemble post-DC DV/TOF matrices and extended winners table.
%   3. Checkpoint after EVERY completed cell (safe to kill and requeue).
%   4. Write final .mat when all cells are done.
%
% Parallelism knobs:
%   N_WORKERS_P1 — parfor workers for footprint build + pair intersection
%   N_WORKERS_P2 — parfor workers for Phase 2 DC  (set 0 for serial DC)
%
% Outputs written to OUTPUT_DIR (rs3_dc_sweep_results/):
%   sweep_dc_checkpoint.mat   — live progress (atomic write after every cell)
%   sweep_dc_results.mat      — final authoritative data store
%
% Output variables in sweep_dc_results.mat:
%   DVmatrix_sweep      [N×N per cell]  post-DC DV     (NaN = no overlap / not converged)
%   TOFmatrix_sweep     [N×N per cell]  post-DC TOF
%   DVmatrix_sweep_proxy [N×N per cell] pre-DC proxy DV (Phase 1 only)
%   TOFmatrix_sweep_proxy [N×N per cell] pre-DC proxy TOF
%   winners_sweep       {nDV×nTmax}     extended winners tables (see below)
%   DV_cap_list, Tmax_list, Tmax_labels, families — sweep axes
%   done_mask, source_sweep              — bookkeeping
%
% winners_sweep{di,dj} table columns:
%   FamilyA, FamilyB
%   minDVproxy_mps, DVlb_mps, DVpatch_ub_mps, EstimatedTOF_days, VoxelId  (Phase 1)
%   DV_total_dc_mps, TOF_dc_days, converged, exitflag, dv_change_mps       (Phase 2)
%
% NOTE: DVmatrix_sweep uses NaN for non-converged pairs (not the proxy),
% so net_build_graph treats them as no-edge.  DVmatrix_sweep_proxy always
% holds the proxy value for comparison.
%
% Drop-in compatible with run_network_centrality_sweep.m (same var names).

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

% Phase 1 workers: footprint build + pair intersection (very low memory).
%   0 → serial;  N → parfor
N_WORKERS_P1 = 4;

% Phase 2 workers: DC per pair (each worker receives 2 subset atlases).
%   0 → serial (safest for memory-limited nodes);  N → parfor
N_WORKERS_P2 = 4;

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

% I/O paths
CHECKPOINT_FILE = fullfile(repoRoot, 'rs3_dc_sweep_results', 'sweep_dc_checkpoint.mat');
OUTPUT_DIR      = fullfile(repoRoot, 'rs3_dc_sweep_results');

% ── Sweep grid ────────────────────────────────────────────────────────────────
DV_cap_list = linspace(0.025, 0.200, 20)';   % 20 values  (0.025 → 0.200)
Tmax_list   = linspace(pi/4,  pi,    20)';   % 20 values  (pi/4  → pi)

Tmax_labels = arrayfun(@(t) strrep(sprintf('Tp%0.3f', t/pi), '.', 'p'), ...
                       Tmax_list, 'UniformOutput', false);

% ══════════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION  (must match cached atlases exactly)
% ══════════════════════════════════════════════════════════════════════════════
cfg = rs3_cfg_defaults();

cfg.families.list      = families;
cfg.families.test_only = false;

cfg.grid.dx               = 0.0005;
cfg.grid.dy               = 0.0005;
cfg.grid.dtheta           = deg2rad(0.5);
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
cfg.cache.dir     = fullfile(repoRoot, 'rs3_cache');
cfg.cache.rebuild = false;

% Suppress all figure/file output
cfg.io.save_figs   = false;
cfg.io.save_fig    = false;
cfg.io.fig_visible = 'off';
cfg.plot.rs4.overlap_xy   = false;
cfg.plot.rs4.overlap_xyz  = false;
cfg.plot.rs4.bounds_lb    = false;
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = false;
% ---- Differential Correction (rs4_diffcorr knobs) ----
cfg.diffcorr.tol_patch     = 1e-4;   % normalized convergence threshold
cfg.diffcorr.tol_converged = 1e-4;   % report CONVERGED if ||r_sc|| <= this (can be >= tol_patch)
cfg.diffcorr.display       = 'off';  % 'off' | 'iter' | 'final'
cfg.diffcorr.MaxIterations = 300;    % fmincon iteration budget (raise to 600 for hard cases)
cfg.diffcorr.MaxFunEvals   = 8000;   % fmincon function eval budget
cfg.diffcorr.N_po_dt       = 0.003;  % PO knot spacing [ND] — auto-scales to orbit period
cfg.diffcorr.N_po_min      = 1001;   % minimum knot count (floor)

if exist('rs3_cfg_validate', 'file') == 2
    rs3_cfg_validate(cfg);
end

% ── Derived constants ─────────────────────────────────────────────────────────
VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
N       = numel(families);
nDV     = numel(DV_cap_list);
nTmax   = numel(Tmax_list);
nPairs  = N * (N - 1) / 2;

% Pre-compute flat pair-index vectors
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

% Ensure output directory exists
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

fprintf('[dc_sweep] Sweep grid: %d DV_cap × %d Tmax = %d cells.\n', nDV, nTmax, nDV*nTmax);
fprintf('[dc_sweep] Checkpoint : %s\n', CHECKPOINT_FILE);
fprintf('[dc_sweep] Phase 1 workers: %d   Phase 2 workers: %d\n\n', N_WORKERS_P1, N_WORKERS_P2);

% ══════════════════════════════════════════════════════════════════════════════
%  LOAD / INITIALISE CHECKPOINT
% ══════════════════════════════════════════════════════════════════════════════
done_mask             = false(nDV, nTmax);
DVmatrix_sweep        = cell(nDV, nTmax);
TOFmatrix_sweep       = cell(nDV, nTmax);
DVmatrix_sweep_proxy  = cell(nDV, nTmax);
TOFmatrix_sweep_proxy = cell(nDV, nTmax);
winners_sweep         = cell(nDV, nTmax);
source_sweep          = cell(nDV, nTmax);

if isfile(CHECKPOINT_FILE)
    fprintf('[dc_sweep] Checkpoint found — loading...\n');
    try
        ck = load(CHECKPOINT_FILE, ...
            'done_mask', 'DVmatrix_sweep', 'TOFmatrix_sweep', ...
            'DVmatrix_sweep_proxy', 'TOFmatrix_sweep_proxy', ...
            'winners_sweep', 'source_sweep');
        if isequal(size(ck.done_mask), [nDV, nTmax])
            done_mask             = ck.done_mask;
            DVmatrix_sweep        = ck.DVmatrix_sweep;
            TOFmatrix_sweep       = ck.TOFmatrix_sweep;
            DVmatrix_sweep_proxy  = ck.DVmatrix_sweep_proxy;
            TOFmatrix_sweep_proxy = ck.TOFmatrix_sweep_proxy;
            winners_sweep         = ck.winners_sweep;
            source_sweep          = ck.source_sweep;
            fprintf('[dc_sweep] Checkpoint loaded. %d/%d cells done.\n', ...
                sum(done_mask(:)), nDV*nTmax);
        else
            warning('[dc_sweep] Checkpoint grid size mismatch — starting fresh.');
        end
    catch ME
        warning('[dc_sweep] Could not load checkpoint (%s) — starting fresh.', ME.message);
    end
end

% Early exit if everything is done
if all(done_mask(:))
    fprintf('[dc_sweep] All cells already done — skipping computation.\n');
else

% ══════════════════════════════════════════════════════════════════════════════
%  LOAD BASE ATLASES  (Tmax = pi, DV_cap = 0.2) with dense PO trace for DC
% ══════════════════════════════════════════════════════════════════════════════
fprintf('[dc_sweep] Loading base atlases (Tmax=π, DV_cap=0.2) with Xpo...\n');
grid3_base = rs3_grid_make(cfg);
Sall_base  = cell(N, 1);
relTol = cfg.propag.relTol;
absTol = cfg.propag.absTol;
for i = 1:N
    fprintf('[dc_sweep]   family %d/%d: %s\n', i, N, families{i});
    [Sall_base{i}, ~] = rs3_prepare_or_load_family(families{i}, cfg, grid3_base);
    % Ensure dense PO trace is present (needed by rs4_diffcorr, inherited by subsets)
    Sall_base{i} = local_ensure_xpo(Sall_base{i}, relTol, absTol, 1001);
end
fprintf('[dc_sweep] Base atlases loaded.\n\n');

% ══════════════════════════════════════════════════════════════════════════════
%  PARALLEL POOL SETUP
% ══════════════════════════════════════════════════════════════════════════════
nW = max(N_WORKERS_P1, N_WORKERS_P2);
if nW > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', nW);
        fprintf('[dc_sweep] Started parpool with %d workers.\n', nW);
    else
        fprintf('[dc_sweep] Using existing parpool (%d workers).\n', pool.NumWorkers);
    end
else
    fprintf('[dc_sweep] N_WORKERS = 0 → fully serial.\n');
end

% ══════════════════════════════════════════════════════════════════════════════
%  MAIN SWEEP LOOP — serial over (DV_cap, Tmax) cells
%  Reversed order (large→small) so workers hit peak memory on first cell only.
% ══════════════════════════════════════════════════════════════════════════════
for di = nDV:-1:1
    for dj = nTmax:-1:1

        if done_mask(di, dj)
            fprintf('[dc_sweep] (%d,%d) DV=%.3f Tmax=%s — already done.\n', ...
                di, dj, DV_cap_list(di), local_tmax_str(Tmax_list(dj)));
            continue;
        end

        dv_cap = DV_cap_list(di);
        tmax   = Tmax_list(dj);
        cellNo = (di - 1) * nTmax + dj;

        fprintf('\n[dc_sweep] ══ cell %d/%d  |  DV_cap=%.3f  Tmax=%s ══\n', ...
            cellNo, nDV*nTmax, dv_cap, local_tmax_str(tmax));
        tCell = tic;

        % ── Build cfg for this cell ───────────────────────────────────────────
        cfg_sub               = cfg;
        cfg_sub.propag.Tmax   = tmax;
        cfg_sub.fan.DV_cap_nd = dv_cap;

        % ── Derive subset atlases ─────────────────────────────────────────────
        fprintf('[dc_sweep] Deriving subset atlases...\n');
        Sall_sub = cell(N, 1);
        if N_WORKERS_P1 > 0
            parfor i = 1:N
                Sall_sub{i} = rs3_atlas_derive_subset(Sall_base{i}, cfg_sub);
            end
        else
            for i = 1:N
                Sall_sub{i} = rs3_atlas_derive_subset(Sall_base{i}, cfg_sub);
            end
        end

        % ── PHASE 1: Build extended footprints ────────────────────────────────
        fprintf('[dc_sweep] Phase 1: building extended footprints...\n');
        Fall_sub = cell(N, 1);
        if N_WORKERS_P1 > 0
            parfor i = 1:N
                Fall_sub{i} = local_compute_footprint_dc( ...
                    Sall_sub{i}, grid3_base, VU_mps, TU_days);
            end
        else
            for i = 1:N
                Fall_sub{i} = local_compute_footprint_dc( ...
                    Sall_sub{i}, grid3_base, VU_mps, TU_days);
            end
        end

        % ── PHASE 1: Pair intersection → tickets ──────────────────────────────
        fprintf('[dc_sweep] Phase 1: pair intersection (%d pairs)...\n', nPairs);
        FA_arr = cell(nPairs, 1);
        FB_arr = cell(nPairs, 1);
        for p2 = 1:nPairs
            FA_arr{p2} = Fall_sub{pairI(p2)};
            FB_arr{p2} = Fall_sub{pairJ(p2)};
        end
        clear Fall_sub

        tickets = cell(nPairs, 1);
        if N_WORKERS_P1 > 0
            parfor p2 = 1:nPairs
                tickets{p2} = local_run_pair_dc( ...
                    FA_arr{p2}, FB_arr{p2}, grid3_base, cfg_sub, VU_mps);
            end
        else
            for p2 = 1:nPairs
                tickets{p2} = local_run_pair_dc( ...
                    FA_arr{p2}, FB_arr{p2}, grid3_base, cfg_sub, VU_mps);
            end
        end
        clear FA_arr FB_arr

        % ── PHASE 2: DC per pair ──────────────────────────────────────────────
        % Pre-slice subset atlases for parfor (each worker gets only 2 atlases).
        fprintf('[dc_sweep] Phase 2: DC (%d pairs, %d workers)...\n', ...
            nPairs, N_WORKERS_P2);

        SA_arr = cell(nPairs, 1);
        SB_arr = cell(nPairs, 1);
        for p2 = 1:nPairs
            SA_arr{p2} = Sall_sub{pairI(p2)};
            SB_arr{p2} = Sall_sub{pairJ(p2)};
        end
        clear Sall_sub

        pair_dv_dc    = nan(nPairs, 1);
        pair_tof_dc   = nan(nPairs, 1);
        pair_converged = false(nPairs, 1);
        pair_exitflag = nan(nPairs, 1);

        cfg_dc = cfg_sub;   % broadcast-friendly copy

        if N_WORKERS_P2 > 0
            parfor p2 = 1:nPairs
                if isnan(tickets{p2}.voxelId), continue; end %#ok<PFBNS>
                [pair_dv_dc(p2), pair_tof_dc(p2), ...
                    pair_converged(p2), pair_exitflag(p2)] = ...
                    local_dc_pair(SA_arr{p2}, SB_arr{p2}, tickets{p2}, ...
                                  cfg_dc, grid3_base, VU_mps, TU_days);
            end
        else
            for p2 = 1:nPairs
                if isnan(tickets{p2}.voxelId), continue; end
                fprintf('[dc_sweep]   pair %d/%d: %-25s → %s\n', ...
                    p2, nPairs, families{pairI(p2)}, families{pairJ(p2)});
                [pair_dv_dc(p2), pair_tof_dc(p2), ...
                    pair_converged(p2), pair_exitflag(p2)] = ...
                    local_dc_pair(SA_arr{p2}, SB_arr{p2}, tickets{p2}, ...
                                  cfg_dc, grid3_base, VU_mps, TU_days);
            end
        end
        clear SA_arr SB_arr

        % ── Extract Phase 1 proxy scalars from tickets ────────────────────────
        pair_minDV   = nan(nPairs, 1);
        pair_DVlb    = nan(nPairs, 1);
        pair_DVpatch = nan(nPairs, 1);
        pair_TOF     = nan(nPairs, 1);
        pair_voxelId = nan(nPairs, 1);
        for p2 = 1:nPairs
            tk = tickets{p2};
            pair_minDV(p2)   = tk.dv_proxy_mps;
            pair_DVlb(p2)    = tk.dvlb_mps;
            pair_DVpatch(p2) = tk.dvpatch_mps;
            pair_TOF(p2)     = tk.tof_days;
            pair_voxelId(p2) = tk.voxelId;
        end
        clear tickets

        % ── Assemble N×N matrices ─────────────────────────────────────────────
        minDVproxyMat   = nan(N, N);
        TOFproxyMat     = nan(N, N);
        minDVdcMat      = nan(N, N);
        TOFdcMat        = nan(N, N);

        for p2 = 1:nPairs
            i = pairI(p2);  j = pairJ(p2);

            minDVproxyMat(i, j) = pair_minDV(p2);
            minDVproxyMat(j, i) = pair_minDV(p2);
            TOFproxyMat(i, j)   = pair_TOF(p2);
            TOFproxyMat(j, i)   = pair_TOF(p2);

            if pair_converged(p2)
                minDVdcMat(i, j) = pair_dv_dc(p2);
                minDVdcMat(j, i) = pair_dv_dc(p2);
                TOFdcMat(i, j)   = pair_tof_dc(p2);
                TOFdcMat(j, i)   = pair_tof_dc(p2);
            end
            % Non-converged → NaN (deliberate: failed DC is not passed to network)
        end

        % ── Build extended winners table ──────────────────────────────────────
        pair_famA      = families(pairI)';
        pair_famB      = families(pairJ)';
        dv_change      = pair_dv_dc - pair_minDV;   % negative = proxy overestimated

        T_win = table(pair_famA(:), pair_famB(:), ...
            pair_minDV, pair_DVlb, pair_DVpatch, pair_TOF, pair_voxelId, ...
            pair_dv_dc, pair_tof_dc, pair_converged, pair_exitflag, dv_change, ...
            'VariableNames', { ...
                'FamilyA', 'FamilyB', ...
                'minDVproxy_mps', 'DVlb_mps', 'DVpatch_ub_mps', ...
                'EstimatedTOF_days', 'VoxelId', ...
                'DV_total_dc_mps', 'TOF_dc_days', ...
                'converged', 'exitflag', 'dv_change_mps'});

        % ── Store in sweep arrays ─────────────────────────────────────────────
        DVmatrix_sweep{di, dj}        = minDVdcMat;
        TOFmatrix_sweep{di, dj}       = TOFdcMat;
        DVmatrix_sweep_proxy{di, dj}  = minDVproxyMat;
        TOFmatrix_sweep_proxy{di, dj} = TOFproxyMat;
        winners_sweep{di, dj}         = T_win;
        source_sweep{di, dj}          = 'computed';
        done_mask(di, dj)             = true;

        % ── Convergence report for this cell ──────────────────────────────────
        has_overlap = isfinite(pair_minDV);
        n_overlap   = sum(has_overlap);
        n_conv      = sum(pair_converged);
        fprintf('[dc_sweep] Cell (%d,%d): %d pairs with overlap, %d converged (%.0f%%), %.1fs\n', ...
            di, dj, n_overlap, n_conv, 100*n_conv/max(1,n_overlap), toc(tCell));

        % ── Atomic checkpoint ─────────────────────────────────────────────────
        local_save_checkpoint(CHECKPOINT_FILE, done_mask, ...
            DVmatrix_sweep, TOFmatrix_sweep, ...
            DVmatrix_sweep_proxy, TOFmatrix_sweep_proxy, ...
            winners_sweep, source_sweep);
        fprintf('[dc_sweep] Checkpoint saved.\n');

    end  % dj
end  % di

clear Sall_base grid3_base

end  % if ~all(done_mask)

% ══════════════════════════════════════════════════════════════════════════════
%  SAVE FINAL .MAT
% ══════════════════════════════════════════════════════════════════════════════
outMat = fullfile(OUTPUT_DIR, 'sweep_dc_results.mat');
save(outMat, ...
    'DV_cap_list', 'Tmax_list', 'Tmax_labels', 'families', ...
    'DVmatrix_sweep', 'TOFmatrix_sweep', ...
    'DVmatrix_sweep_proxy', 'TOFmatrix_sweep_proxy', ...
    'winners_sweep', 'source_sweep', 'done_mask', ...
    '-v7.3');
fprintf('\n[dc_sweep] Final .mat saved:\n  %s\n', outMat);

% Print overall convergence summary
n_cells_done = sum(done_mask(:));
fprintf('\n[dc_sweep] ══════════ SWEEP COMPLETE ══════════\n');
fprintf('  Cells completed : %d / %d\n', n_cells_done, nDV*nTmax);
fprintf('  Output dir      : %s\n', OUTPUT_DIR);
fprintf('  DVmatrix_sweep  → post-DC DV  (NaN = no overlap or not converged)\n');
fprintf('  DVmatrix_sweep_proxy → pre-DC proxy DV (Phase 1 only)\n');
fprintf('  Use run_network_centrality_sweep.m with sweep_dc_results.mat\n');
fprintf('  to compare DC-corrected vs proxy-based network metrics.\n');

% ══════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────────────
function F = local_compute_footprint_dc(S, grid3, VU_mps, TU_days)
%LOCAL_COMPUTE_FOOTPRINT_DC  Extended footprint with argmin row identifiers.
%
% Same as local_compute_footprint in run_rs4_dv_tmax_sweep, but also stores
% per-voxel argmin (iSeed, iHead, halfFlag, t_nd) for Phase 2 DC warm start.
%
% Fields added vs. the proxy footprint:
%   F.iSeed_argmin_frs / F.iHead_argmin_frs / F.halfFlag_argmin_frs / F.t_nd_argmin_frs
%   Same 4 fields for BRS (referring to the FRS row that was mirrored).

Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

% Theta mirror LUT (identical to rs4_overlap_pair)
thm    = rs3_wrapToPi(pi - grid3.th_centers(:));
lut    = discretize(thm, grid3.th_edges);
lut(isnan(lut)) = 0;
it_lut = uint16(lut);

% Delta-angle lookup matrix (vectorised, avoids per-row cell indexing)
dlists = S.Step4.delta_lists;
Ns     = numel(dlists);
max_h  = max(1, max(cellfun(@numel, dlists)));
delta_mat = zeros(Ns, max_h);
for s = 1:Ns
    v = double(dlists{s});
    delta_mat(s, 1:numel(v)) = v;
end

% v0 per unique seed (avoids per-row potential evaluation)
pot_u    = rs3_core_cr3bp_U_and_derivs(S.SeedsUpper(:,1), S.SeedsUpper(:,2), S.mu);
v0_upper = sqrt(max(2*pot_u.U(:) - S.CJ, 0));

if isfield(S, 'SeedsLower') && ~isempty(S.SeedsLower)
    pot_l    = rs3_core_cr3bp_U_and_derivs(S.SeedsLower(:,1), S.SeedsLower(:,2), S.mu);
    v0_lower = sqrt(max(2*pot_l.U(:) - S.CJ, 0));
else
    v0_lower = v0_upper;
end

% ── Process FRS_upper rows ────────────────────────────────────────────────────
nu = double(S.Step4.rows_FRS_upper.n);
if nu > 0
    [ids_u, dv_u, t_u, ix_u, iy_u, it_u, iSeed_u, iHead_u, t_nd_u] = local_fp_rows_dc( ...
        S.Step4.rows_FRS_upper, nu, v0_upper, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
    half_u = repmat(int8(1), nu, 1);
else
    ids_u=[]; dv_u=[]; t_u=[]; ix_u=[]; iy_u=[]; it_u=[];
    iSeed_u=[]; iHead_u=[]; t_nd_u=[]; half_u=int8([]);
end

% ── Process FRS_lower rows ────────────────────────────────────────────────────
nl = double(S.Step4.rows_FRS_lower.n);
if nl > 0
    [ids_l, dv_l, t_l, ix_l, iy_l, it_l, iSeed_l, iHead_l, t_nd_l] = local_fp_rows_dc( ...
        S.Step4.rows_FRS_lower, nl, v0_lower, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
    half_l = repmat(int8(-1), nl, 1);
else
    ids_l=[]; dv_l=[]; t_l=[]; ix_l=[]; iy_l=[]; it_l=[];
    iSeed_l=[]; iHead_l=[]; t_nd_l=[]; half_l=int8([]);
end

% ── Aggregate FRS voxels ──────────────────────────────────────────────────────
[F.uid_frs, F.dv_min_frs, F.t_mean_frs, ...
 F.iSeed_argmin_frs, F.iHead_argmin_frs, F.halfFlag_argmin_frs, F.t_nd_argmin_frs] = ...
    local_fp_agg_dc([ids_u; ids_l], [dv_u; dv_l], [t_u; t_l], ...
                    [iSeed_u; iSeed_l], [iHead_u; iHead_l], ...
                    [half_u;  half_l],  [t_nd_u;  t_nd_l]);

% ── BRS: mirror FRS_upper ─────────────────────────────────────────────────────
if ~isempty(ix_u)
    biy_u   = Ny - iy_u + 1;
    bit_u   = double(it_lut(it_u));
    ok_u    = bit_u > 0;
    ids_brs_u = sub2ind([Ny,Nx,Nt], biy_u(ok_u), ix_u(ok_u), max(1, min(Nt, bit_u(ok_u))));
    dv_bu  = dv_u(ok_u);    t_bu    = t_u(ok_u);    t_nd_bu  = t_nd_u(ok_u);
    iS_bu  = iSeed_u(ok_u); iH_bu   = iHead_u(ok_u); half_bu  = half_u(ok_u);
else
    ids_brs_u=[]; dv_bu=[]; t_bu=[]; t_nd_bu=[]; iS_bu=[]; iH_bu=[]; half_bu=int8([]);
end

% ── BRS: mirror FRS_lower ─────────────────────────────────────────────────────
if ~isempty(ix_l)
    biy_l   = Ny - iy_l + 1;
    bit_l   = double(it_lut(it_l));
    ok_l    = bit_l > 0;
    ids_brs_l = sub2ind([Ny,Nx,Nt], biy_l(ok_l), ix_l(ok_l), max(1, min(Nt, bit_l(ok_l))));
    dv_bl  = dv_l(ok_l);    t_bl    = t_l(ok_l);    t_nd_bl  = t_nd_l(ok_l);
    iS_bl  = iSeed_l(ok_l); iH_bl   = iHead_l(ok_l); half_bl  = half_l(ok_l);
else
    ids_brs_l=[]; dv_bl=[]; t_bl=[]; t_nd_bl=[]; iS_bl=[]; iH_bl=[]; half_bl=int8([]);
end

% ── Aggregate BRS voxels ──────────────────────────────────────────────────────
% Argmin row identifiers for BRS refer to the FRS row that was mirrored.
% Phase 2 uses halfFlag to pick SeedsUpper (+1) or SeedsLower (-1) from SB.
[F.uid_brs, F.dv_min_brs, F.t_mean_brs, ...
 F.iSeed_argmin_brs, F.iHead_argmin_brs, F.halfFlag_argmin_brs, F.t_nd_argmin_brs] = ...
    local_fp_agg_dc([ids_brs_u; ids_brs_l], [dv_bu; dv_bl], [t_bu; t_bl], ...
                    [iS_bu; iS_bl], [iH_bu; iH_bl], [half_bu; half_bl], [t_nd_bu; t_nd_bl]);

F.CJ   = S.CJ;
F.mu   = S.mu;
F.name = S.name;
end

% ─────────────────────────────────────────────────────────────────────────────
function [ids, dv_mps, t_days, ix_out, iy_out, it_out, iSeed_out, iHead_out, t_nd_out] = ...
        local_fp_rows_dc(rows, n, v0_per_seed, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days)
%LOCAL_FP_ROWS_DC  Extract voxel IDs, DV, TOF, and row identifiers.
% Extends local_fp_rows (from run_rs4_dv_tmax_sweep) with iSeed, iHead, t_nd outputs.

ix_out    = double(rows.ix(1:n));
iy_out    = double(rows.iy(1:n));
it_out    = double(rows.it(1:n));
ids       = sub2ind([Ny, Nx, Nt], iy_out, ix_out, it_out);
iSeed_out = double(rows.iSeed(1:n));
iHead_out = double(rows.iHead(1:n));
t_nd_out  = double(rows.t(1:n));   % signed nondim integration time (for ODE warm start)

lin    = sub2ind([Ns, max_h], iSeed_out, iHead_out);
delta  = delta_mat(lin);
v0     = v0_per_seed(iSeed_out);
dv_mps = 2 * v0(:) .* sin(abs(delta(:)) / 2) * VU_mps;
t_days = abs(t_nd_out(:)) * TU_days;
end

% ─────────────────────────────────────────────────────────────────────────────
function [uid, dv_min, t_mean, iSeed_win, iHead_win, half_win, t_nd_win] = ...
        local_fp_agg_dc(ids, dv, t, iSeed_arr, iHead_arr, half_arr, t_nd_arr)
%LOCAL_FP_AGG_DC  Per-voxel aggregation with argmin row tracking.
% Returns per-voxel min DV, mean TOF, and the row (iSeed,iHead,halfFlag,t_nd)
% that achieved the minimum DV for that voxel.

if isempty(ids)
    uid=zeros(0,1); dv_min=zeros(0,1); t_mean=zeros(0,1);
    iSeed_win=uint16(zeros(0,1)); iHead_win=uint16(zeros(0,1));
    half_win=int8(zeros(0,1));    t_nd_win=single(zeros(0,1));
    return;
end
[uid, ~, ic] = unique(ids(:));

% Argmin per voxel group: sort by (ic, dv) ascending → first per group = argmin
[~, srt] = sortrows([ic(:), dv(:)]);
[~, first_in_group] = unique(ic(srt), 'first');
win = srt(first_in_group);   % indices into original input arrays

dv_min    = dv(win);
t_mean    = accumarray(ic, t(:)) ./ accumarray(ic, ones(numel(ic), 1));
iSeed_win = uint16(iSeed_arr(win));
iHead_win = uint16(iHead_arr(win));
half_win  = int8(half_arr(win));
t_nd_win  = single(t_nd_arr(win));
end

% ─────────────────────────────────────────────────────────────────────────────
function ticket = local_run_pair_dc(FA, FB, grid3, cfg, VU_mps)
%LOCAL_RUN_PAIR_DC  Phase 1 pair intersection returning full ticket.
%
% Same overlap logic as local_run_pair (run_rs4_dv_tmax_sweep) plus argmin
% row identifiers from the extended footprints for Phase 2 DC warm start.

ticket = struct('voxelId',NaN, 'dv_proxy_mps',NaN, 'dvlb_mps',NaN, ...
                'dvpatch_mps',NaN, 'tof_days',NaN, ...
                'iSeed_A',0,      'iHead_A',0,      'halfFlag_A',int8(0), 't_nd_A',single(0), ...
                'iSeed_B',0,      'iHead_B',0,      'halfFlag_B',int8(0), 't_nd_B',single(0));
try
    % 1. Intersect FRS(A) ∩ BRS(B) (both sorted unique)
    idsO = intersect(FA.uid_frs, FB.uid_brs);
    if isempty(idsO), return; end

    % 2. Unpack grid indices
    Ny = numel(grid3.y_centers);
    Nx = numel(grid3.x_centers);
    Nt = numel(grid3.th_centers);
    [iy, ix, ~] = ind2sub([Ny, Nx, Nt], idsO);

    % 3. Keep + primary-buffer filter (mirrors rs4_overlap_pair)
    bufFrac = 0.05;
    if isfield(cfg,'overlap') && isfield(cfg.overlap,'primary_buffer_frac') && ...
            ~isempty(cfg.overlap.primary_buffer_frac)
        bufFrac = cfg.overlap.primary_buffer_frac;
    end
    RE = cfg.sys.RE_nd;
    RM = cfg.sys.RM_nd;
    mu = FA.mu;

    if isfield(grid3,'Keep') && ~isempty(grid3.Keep)
        keepXY = logical(grid3.Keep);
        if ~isequal(size(keepXY), [Ny, Nx]), keepXY = keepXY.'; end
        okKeep = keepXY(sub2ind([Ny, Nx], iy, ix));
    else
        okKeep = true(numel(idsO), 1);
    end
    x = grid3.x_centers(ix);
    y = grid3.y_centers(iy);
    okE = hypot(x + mu, y)     > (1 + bufFrac) * RE;
    okM = hypot(x - (1-mu), y) > (1 + bufFrac) * RM;
    ok  = okKeep(:) & okE(:) & okM(:);

    idsO = idsO(ok);
    if isempty(idsO), return; end
    ix = ix(ok);  iy = iy(ok);

    % 4. Lookup pre-computed per-voxel DV and TOF
    [~, locA] = ismember(idsO, FA.uid_frs);
    [~, locB] = ismember(idsO, FB.uid_brs);
    dv_min_A  = FA.dv_min_frs(locA);
    dv_min_B  = FB.dv_min_brs(locB);
    t_mean_A  = FA.t_mean_frs(locA);
    t_mean_B  = FB.t_mean_brs(locB);

    % 5. DV proxy (identical formula to run_rs4_dv_tmax_sweep)
    x_ok = grid3.x_centers(ix);
    y_ok = grid3.y_centers(iy);
    CJstar = min(FA.CJ, FB.CJ);
    pot    = rs3_core_cr3bp_U_and_derivs(x_ok(:), y_ok(:), mu);
    v_box  = sqrt(max(2 * pot.U - CJstar, 0));
    dv_patch_vec = 2 * v_box .* sin(abs(grid3.dtheta) / 2) * VU_mps;
    dv_lb_vec    = dv_min_A(:) + dv_min_B(:);
    dv_proxy     = dv_lb_vec + dv_patch_vec;

    valid = isfinite(dv_proxy);
    if ~any(valid), return; end

    idxValid  = find(valid);
    [~, iLoc] = min(dv_proxy(idxValid));
    iWin      = idxValid(iLoc);

    % 6. Fill scalar metrics
    ticket.voxelId      = idsO(iWin);
    ticket.dv_proxy_mps = dv_proxy(iWin);
    ticket.dvlb_mps     = dv_lb_vec(iWin);
    ticket.dvpatch_mps  = dv_patch_vec(iWin);
    ticket.tof_days     = t_mean_A(iWin) + t_mean_B(iWin);

    % 7. Argmin row identifiers for Phase 2
    wA = locA(iWin);
    ticket.iSeed_A    = double(FA.iSeed_argmin_frs(wA));
    ticket.iHead_A    = double(FA.iHead_argmin_frs(wA));
    ticket.halfFlag_A = double(FA.halfFlag_argmin_frs(wA));
    ticket.t_nd_A     = double(FA.t_nd_argmin_frs(wA));

    wB = locB(iWin);
    ticket.iSeed_B    = double(FB.iSeed_argmin_brs(wB));
    ticket.iHead_B    = double(FB.iHead_argmin_brs(wB));
    ticket.halfFlag_B = double(FB.halfFlag_argmin_brs(wB));
    ticket.t_nd_B     = double(FB.t_nd_argmin_brs(wB));

catch ME
    warning('[dc_sweep:p1] %s', ME.message);
end
end

% ─────────────────────────────────────────────────────────────────────────────
function [dv_dc, tof_dc, converged, exitflag] = ...
        local_dc_pair(SA, SB, ticket, cfg, grid3, VU_mps, TU_days)
%LOCAL_DC_PAIR  Phase 2 DC for one pair. Returns scalar metrics only.
% Reconstructs the warm-start T struct from the ticket, then calls rs4_diffcorr.

dv_dc = NaN;  tof_dc = NaN;  converged = false;  exitflag = NaN;
try
    T  = local_traj_from_ticket(SA, SB, ticket, cfg, grid3, VU_mps, TU_days);
    if isempty(T), return; end
    Tc = rs4_diffcorr(T, SA, SB, cfg);
    dv_dc    = Tc.DV_total_mps;
    tof_dc   = Tc.tof_A_days + Tc.tof_B_days;
    converged = Tc.converged;
    exitflag = Tc.exitflag;
catch ME
    warning('[dc_sweep:p2] %s', ME.message);
end
end

% ─────────────────────────────────────────────────────────────────────────────
function T = local_traj_from_ticket(SA, SB, ticket, cfg, grid3, VU_mps, TU_days)
%LOCAL_TRAJ_FROM_TICKET  Reconstruct T struct (warm start for rs4_diffcorr).
%
% Replaces rs4_overlap_pair + rs4_overlap_extract_voxel_info +
%          rs4_overlap_visualize_bounds + rs4_voxel_traj_extract
% with a direct, targeted arc reconstruction using only the winning rows.
%
% The ticket identifies:
%   Side A: argmin FRS row of SA  (iSeed_A, iHead_A, halfFlag_A, t_nd_A)
%   Side B: argmin FRS row of SB  (iSeed_B, iHead_B, halfFlag_B, t_nd_B)
%           → integrated as FRS, then mirrored to BRS
%
% Note: iHead values reference the atlas that was used to build the footprint
% (subset atlas in the sweep, full atlas in run_rs4_dc_sweep). Phase 2 must
% use the same atlas to ensure iHead → delta_lists mapping is consistent.

T = [];

Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

[iy0, ix0, it0] = ind2sub([Ny, Nx, Nt], ticket.voxelId);
T.vid = ticket.voxelId;
T.xc  = grid3.x_centers(ix0);
T.yc  = grid3.y_centers(iy0);
T.thc = grid3.th_centers(it0);

opts = odeset('RelTol', cfg.propag.relTol, 'AbsTol', cfg.propag.absTol);

% ── Side A: FRS of SA ────────────────────────────────────────────────────────
iSeed_A = ticket.iSeed_A;
iHead_A = ticket.iHead_A;
halfA   = ticket.halfFlag_A;
t_nd_A  = ticket.t_nd_A;

if halfA >= 0   % upper half (or unspecified)
    seed_A = SA.SeedsUpper(iSeed_A, :);
else
    if isfield(SA,'SeedsLower') && size(SA.SeedsLower,1) >= iSeed_A
        seed_A = SA.SeedsLower(iSeed_A, :);
    else
        seed_A = SA.SeedsUpper(iSeed_A, :);
    end
end
delta_A = double(SA.Step4.delta_lists{iSeed_A}(iHead_A));
IC_A    = [seed_A(1); seed_A(2); rs3_wrapToPi(seed_A(3) + delta_A)];

t_end_A = max(abs(t_nd_A), 1e-4);
tA_span = linspace(0, t_end_A, 501)';
sol_A   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,SA.CJ,SA.mu,false), ...
                  [0, t_end_A], IC_A, opts);
XA = deval(sol_A, tA_span)';

pot_A = rs3_core_cr3bp_U_and_derivs(seed_A(1), seed_A(2), SA.mu);
v0_A  = sqrt(max(2*pot_A.U - SA.CJ, 0));
T.DV_turn_A_mps = 2 * v0_A * sin(abs(delta_A)/2) * VU_mps;
T.IC_A   = IC_A;     T.seed_A  = seed_A;
T.t_A    = t_end_A;  T.tA_vec  = tA_span;  T.XA = XA;
T.iSeed_A = iSeed_A; T.iHead_A = iHead_A;  T.halfFlag_A = halfA;

% ── Side B: FRS of SB, then mirror → BRS ─────────────────────────────────────
iSeed_B = ticket.iSeed_B;
iHead_B = ticket.iHead_B;
halfB   = ticket.halfFlag_B;
t_nd_B  = ticket.t_nd_B;

if halfB >= 0
    seed_B_frs = SB.SeedsUpper(iSeed_B, :);
else
    if isfield(SB,'SeedsLower') && size(SB.SeedsLower,1) >= iSeed_B
        seed_B_frs = SB.SeedsLower(iSeed_B, :);
    else
        seed_B_frs = SB.SeedsUpper(iSeed_B, :);
    end
end
delta_B  = double(SB.Step4.delta_lists{iSeed_B}(iHead_B));
IC_B_frs = [seed_B_frs(1); seed_B_frs(2); rs3_wrapToPi(seed_B_frs(3) + delta_B)];

t_end_B = max(abs(t_nd_B), 1e-4);
tB_span = linspace(0, t_end_B, 501)';
sol_B   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,SB.CJ,SB.mu,false), ...
                  [0, t_end_B], IC_B_frs, opts);
XB_frs = deval(sol_B, tB_span)';

% Mirror FRS → BRS: (x, y, θ) → (x, -y, π-θ)
x_B  = XB_frs(:,1);
y_B  = -XB_frs(:,2);
th_B = rs3_wrapToPi(pi - XB_frs(:,3));

pot_B = rs3_core_cr3bp_U_and_derivs(seed_B_frs(1), seed_B_frs(2), SB.mu);
v0_B  = sqrt(max(2*pot_B.U - SB.CJ, 0));
T.DV_turn_B_mps  = 2 * v0_B * sin(abs(delta_B)/2) * VU_mps;
T.IC_B_frs       = IC_B_frs;  T.seed_B_frs  = seed_B_frs;
T.t_B            = t_end_B;   T.tB_vec      = tB_span;
T.x_B = x_B;  T.y_B = y_B;  T.th_B = th_B;
T.iSeed_B = iSeed_B;  T.iHead_B = iHead_B;  T.halfFlag_B_frs = halfB;
T.from_lower_B = (halfB < 0);

% ── Closest approach to voxel center ─────────────────────────────────────────
xc = T.xc;  yc = T.yc;
dA_vec = hypot(XA(:,1) - xc, XA(:,2) - yc);
[T.dA_nd, T.i_star] = min(dA_vec);
T.th_A_star = XA(T.i_star, 3);

dB_vec = hypot(x_B - xc, y_B - yc);
[T.dB_nd, T.j_star] = min(dB_vec);
T.th_B_star = th_B(T.j_star);

T.delta_th_rad = abs(rs3_circ_diff(T.th_A_star, T.th_B_star));

% ── DV patch and totals ───────────────────────────────────────────────────────
CJstar = min(SA.CJ, SB.CJ);
pot_c  = rs3_core_cr3bp_U_and_derivs(xc, yc, SA.mu);
T.v_box_center_nd   = sqrt(max(2*pot_c.U - CJstar, 0));
T.DV_patch_nd       = 2 * T.v_box_center_nd * sin(T.delta_th_rad / 2);
T.DV_patch_mps      = T.DV_patch_nd * VU_mps;
T.DV_total_true_mps = T.DV_turn_A_mps + T.DV_patch_mps + T.DV_turn_B_mps;
T.DV_proxy_mps      = ticket.dv_proxy_mps;
T.tof_A_days        = abs(t_nd_A) * TU_days;
T.tof_B_days        = abs(t_nd_B) * TU_days;
end

% ─────────────────────────────────────────────────────────────────────────────
function local_save_checkpoint(fpath, done_mask, ...
        DVmatrix_sweep, TOFmatrix_sweep, ...
        DVmatrix_sweep_proxy, TOFmatrix_sweep_proxy, ...
        winners_sweep, source_sweep)
%LOCAL_SAVE_CHECKPOINT  Atomic checkpoint save.
tmp = [fpath '.tmp'];
save(tmp, ...
    'done_mask', 'DVmatrix_sweep', 'TOFmatrix_sweep', ...
    'DVmatrix_sweep_proxy', 'TOFmatrix_sweep_proxy', ...
    'winners_sweep', 'source_sweep', '-v7.3');
if isfile(fpath), delete(fpath); end
movefile(tmp, fpath);
end

% ─────────────────────────────────────────────────────────────────────────────
function S = local_ensure_xpo(S, relTol, absTol, N_po)
%LOCAL_ENSURE_XPO  Re-integrate PO if Xpo was stripped from cache.
if isfield(S, 'Xpo') && ~isempty(S.Xpo) && ...
   isfield(S, 't_dense') && ~isempty(S.t_dense)
    return;
end
opts    = odeset('RelTol', relTol, 'AbsTol', absTol);
solPO   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,S.CJ,S.mu,true), ...
                  [0, S.Tf_PO], S.X0, opts);
t_dense = linspace(0, S.Tf_PO, N_po)';
S.t_dense = t_dense;
S.Xpo     = deval(solPO, t_dense)';
end

% ─────────────────────────────────────────────────────────────────────────────
function s = local_tmax_str(tmax)
known_val = [pi/4,   pi/3,   pi/2,   2*pi/3, 3*pi/4, pi  ];
known_str = {'pi/4', 'pi/3', 'pi/2', '2pi/3','3pi/4','pi' };
[~, idx] = min(abs(tmax - known_val));
if abs(tmax - known_val(idx)) < 1e-9
    s = known_str{idx};
else
    s = sprintf('%.5g', tmax);
end
end
