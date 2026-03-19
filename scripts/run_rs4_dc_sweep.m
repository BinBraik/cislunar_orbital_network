%% RUN_RS4_DC_SWEEP
% Two-phase batch differential correction over every pair of the 13
% periodic-orbit families.
%
% PHASE 1  — Footprint-based overlap (tiny memory, parfor-safe):
%   Build compact per-family footprints (~5-25 MB each).
%   Run pairwise intersection in parfor → one "ticket" per pair containing
%   the winner voxel ID and the argmin row identifiers (iSeed, iHead,
%   halfFlag, t_nd) for both sides.  No full atlas data in workers.
%
% PHASE 2  — Targeted DC (2 families at a time):
%   For each pair, use the ticket to re-integrate only the two winning arcs
%   (bypassing rs4_overlap_pair / rs4_overlap_extract_voxel_info / etc.),
%   then run rs4_diffcorr.  Each worker holds exactly 2 family atlases.
%
% Two modes (set NUM_WORKERS below):
%
%   SERIAL  (NUM_WORKERS <= 1)
%     Phase 1 runs serially with families in memory.
%     Phase 2 runs serially with checkpoint after every pair.
%
%   PARALLEL (NUM_WORKERS >= 2)
%     Phase 1 parfor uses only footprints (broadcast ~25 MB each, not GB).
%     Phase 2 uses process-based parpool; each worker loads only its 2
%     needed families from disk — no large struct broadcast, no OOM.
%     Checkpoint saved after Phase 2 completes (Phase 1 is fast).
%
% Outputs  (saved to rs3_results/<tag>/rs4_dc_sweep/rs4_dc_sweep_results.mat)
%   before_mat   [nPairs × 13]  pre-DC metrics
%   after_mat    [nPairs × 14]  post-DC metrics + solver metadata
%   traj_cell    {nPairs × 1}   struct per pair — sufficient for plot replay
%   pairs_ij     [nPairs × 2]   family index pairs (i < j)
%   families     {13 × 1}       family name strings
%   BEFORE_COLS / AFTER_COLS    column-name cell arrays
%   T_summary                   table: per-pair proxy vs DC comparison
%   cfg                         config used for this run

clear; clc;

% ── path setup ────────────────────────────────────────────────────────────────
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% ════════════════════════════════════════════════════════════════════════════
% USER SETTINGS
% ════════════════════════════════════════════════════════════════════════════
NUM_WORKERS = 0;   % 0 or 1 = serial with Phase 2 checkpoint
                   % >= 2    = Phase 1 parfor + Phase 2 process-based parfor

families = {
    'Lyapunov L1',
    'Lyapunov L2',
    'Cycler 21',
    'Cycler 11a',
    'Cycler 11b',
    'Cycler 32',
    'Resonant 2to1 Stable',
    'Resonant 2to1 Unstable',
    'Resonant 3to1 Stable',
    'Resonant 3to1 Unstable',
    'Resonant 5to2 Stable',
    'Resonant 5to2 Unstable',
    'Distant Prograde Orbit'
    };

cfg = rs3_cfg_defaults();
cfg.cache.enable          = true;
cfg.cache.rebuild         = false;

cfg.io.save_figs          = false;
cfg.io.save_fig           = false;
cfg.io.fig_visible        = 'off';

% ---- Differential Correction (rs4_diffcorr knobs) ----
cfg.diffcorr.tol_patch     = 1e-4;   % normalized convergence threshold
cfg.diffcorr.tol_converged = 1e-4;   % report CONVERGED if ||r_sc|| <= this (can be >= tol_patch)
cfg.diffcorr.display       = 'off';  % 'off' | 'iter' | 'final'
cfg.diffcorr.MaxIterations = 300;    % fmincon iteration budget (raise to 600 for hard cases)
cfg.diffcorr.MaxFunEvals   = 8000;   % fmincon function eval budget
cfg.diffcorr.N_po_dt       = 0.003;  % PO knot spacing [ND] — auto-scales to orbit period
cfg.diffcorr.N_po_min      = 1001;   % minimum knot count (floor)

% Cache settings — must match what was used to build the caches
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

rs3_cfg_validate(cfg);

N      = numel(families);
VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
relTol  = cfg.propag.relTol;
absTol  = cfg.propag.absTol;

% ── enumerate pairs ───────────────────────────────────────────────────────────
pairs_ij = zeros(0, 2);
for ii = 1:N
    for jj = ii+1:N
        pairs_ij(end+1, :) = [ii, jj]; %#ok<AGROW>
    end
end
nPairs = size(pairs_ij, 1);
pairI  = pairs_ij(:, 1);
pairJ  = pairs_ij(:, 2);

use_par = NUM_WORKERS >= 2 && ...
          license('test', 'Distrib_Computing_Toolbox') && ...
          nPairs > 1;

% ── output folder ─────────────────────────────────────────────────────────────
dcRoot = fullfile(cfg.io.out_root, cfg.io.tag, 'rs4_dc_sweep');
if ~exist(dcRoot, 'dir'), mkdir(dcRoot); end
fprintf('[dc_sweep] Output   : %s\n', dcRoot);
fprintf('[dc_sweep] Mode     : %s\n\n', ...
    ternary(use_par, sprintf('PARALLEL (%d workers)', min(NUM_WORKERS,nPairs)), 'SERIAL'));

% ════════════════════════════════════════════════════════════════════════════
% 1. WARM CACHE + BUILD PHASE 1 FOOTPRINTS  (serial, main process)
%    rs3_prepare_or_load_family uses Java for MD5 — must run in main process.
%    Footprints (~25 MB/family) replace the full-atlas broadcast that caused OOM.
% ════════════════════════════════════════════════════════════════════════════
grid3 = rs3_grid_make(cfg);

fprintf('[dc_sweep] Warming %d family caches and building footprints...\n', N);
cache_fpaths = cell(N, 1);
Fall         = cell(N, 1);
Sall         = cell(N, 1);   % kept only in serial mode for Phase 2

for k = 1:N
    fprintf('  [%2d/%2d]  %s ...', k, N, families{k});
    tL = tic;
    [Stmp, ~] = rs3_prepare_or_load_family(families{k}, cfg, grid3);
    info = rs3_cache_get_path(families{k}, Stmp.mu, Stmp.CJ, cfg);
    cache_fpaths{k} = info.fpath;

    % Ensure dense PO trace (needed by rs4_diffcorr)
    Stmp = local_ensure_xpo(Stmp, relTol, absTol, 1001);

    % Build extended footprint (Phase 1 — tiny, ~25 MB each)
    Fall{k} = local_compute_footprint_dc(Stmp, grid3, VU_mps, TU_days);

    Sall{k} = Stmp;   % keep for Phase 2 (thread pool shares memory, no copy)
    fprintf('  %.1fs  footprint: %d FRS voxels, %d BRS voxels\n', toc(tL), ...
        numel(Fall{k}.uid_frs), numel(Fall{k}.uid_brs));
    clear Stmp info;
end

fprintf('[dc_sweep] Footprints built. Families in memory for Phase 2.\n\n');

% ════════════════════════════════════════════════════════════════════════════
% 2. PHASE 1  — Pair intersection via footprints  (fast, low memory)
%    All 78 pairs processed in parfor using only footprints.
%    Each footprint is ~5-25 MB; workers never see full atlas structs.
% ════════════════════════════════════════════════════════════════════════════
fprintf('[dc_sweep] Phase 1: footprint pair intersection (%d pairs)...\n', nPairs);
tP1 = tic;

FA_arr = cell(nPairs, 1);
FB_arr = cell(nPairs, 1);
for p = 1:nPairs
    FA_arr{p} = Fall{pairI(p)};
    FB_arr{p} = Fall{pairJ(p)};
end
clear Fall

tickets = cell(nPairs, 1);
if use_par
    % Lightweight pool (threads or shared) for Phase 1 — footprints only
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', min(NUM_WORKERS, nPairs));
    end
    parfor p = 1:nPairs
        tickets{p} = local_run_pair_dc(FA_arr{p}, FB_arr{p}, grid3, cfg, VU_mps);
    end
else
    for p = 1:nPairs
        tickets{p} = local_run_pair_dc(FA_arr{p}, FB_arr{p}, grid3, cfg, VU_mps);
    end
end
clear FA_arr FB_arr

n_overlap = sum(cellfun(@(tk) isfinite(tk.voxelId), tickets));
fprintf('[dc_sweep] Phase 1 done in %.1fs — %d/%d pairs have overlap.\n\n', ...
    toc(tP1), n_overlap, nPairs);

% ════════════════════════════════════════════════════════════════════════════
% 3. PRE-ALLOCATE OUTPUT
% ════════════════════════════════════════════════════════════════════════════
BEFORE_COLS = { ...
    'pair_idx',      'fam_i',        'fam_j',         'vid', ...
    'DV_turn_A_mps', 'DV_patch_mps', 'DV_turn_B_mps', 'DV_total_mps', ...
    'dA_nd',         'dB_nd',        'delta_th_deg', ...
    'tof_A_days',    'tof_B_days'};                              % 13 cols

AFTER_COLS = { ...
    'pair_idx',      'fam_i',        'fam_j',         'vid', ...
    'DV_turn_A_mps', 'DV_patch_mps', 'DV_turn_B_mps', 'DV_total_mps', ...
    'delta_th_deg',  'tof_A_days',   'tof_B_days', ...
    'exitflag',      'converged',    'r_norm'};                   % 14 cols

nBC = numel(BEFORE_COLS);
nAC = numel(AFTER_COLS);

before_mat = NaN(nPairs, nBC);
after_mat  = NaN(nPairs, nAC);
traj_cell  = cell(nPairs, 1);

% ════════════════════════════════════════════════════════════════════════════
% 4. PHASE 2  — Targeted DC  (2 families at a time per worker)
%    Workers only re-integrate the 2 argmin arcs from the ticket — no call
%    to rs4_overlap_pair or rs4_overlap_extract_voxel_info.
% ════════════════════════════════════════════════════════════════════════════
ckptFile = fullfile(dcRoot, 'rs4_dc_sweep_checkpoint.mat');

if use_par
    % ── PARALLEL Phase 2 ─────────────────────────────────────────────────────
    % Reuse the existing 'local' (thread) pool from Phase 1 — identical pattern
    % to run_rs4_dc_sweep_tmax.m.  Slice SA_arr/SB_arr so MATLAB sends only
    % each worker's 2 atlas slices via shared memory; no disk I/O, no OOM.
    fprintf('[dc_sweep] Phase 2: DC in parallel (%d workers)...\n', min(NUM_WORKERS,nPairs));

    SA_arr = cell(nPairs, 1);
    SB_arr = cell(nPairs, 1);
    for p = 1:nPairs
        SA_arr{p} = Sall{pairI(p)};
        SB_arr{p} = Sall{pairJ(p)};
    end
    clear Sall;

    cfg_p   = cfg;
    grid3_p = grid3;
    fams    = families;
    pij     = pairs_ij;
    VU      = VU_mps;
    TU      = TU_days;
    tks     = tickets;

    parfor p = 1:nPairs
        ii = pij(p, 1);
        jj = pij(p, 2);

        if ~isfinite(tks{p}.voxelId)
            % No overlap — leave NaN rows
            before_mat(p, 1) = p; %#ok<PFPIE>  % sentinel
            continue;
        end

        [br, ar, tr] = local_process_pair_dc( ...
            SA_arr{p}, SB_arr{p}, fams{ii}, fams{jj}, ii, jj, p, cfg_p, grid3_p, tks{p}, VU, TU);

        before_mat(p, :) = br; %#ok<PFPIE>
        after_mat(p, :)  = ar; %#ok<PFPIE>
        traj_cell{p}     = tr; %#ok<PFPIE>
    end

else
    % ── SERIAL Phase 2  (with checkpoint) ────────────────────────────────────
    fprintf('[dc_sweep] Phase 2: DC serial with checkpoint...\n');
    fprintf('[dc_sweep] Checkpoint: %s\n\n', ckptFile);

    if exist(ckptFile, 'file')
        ck = load(ckptFile, 'before_mat', 'after_mat', 'traj_cell');
        before_mat = ck.before_mat;
        after_mat  = ck.after_mat;
        traj_cell  = ck.traj_cell;
        n_done = sum(isfinite(before_mat(:, 1)));
        fprintf('[dc_sweep] Resumed from checkpoint (%d / %d done).\n\n', n_done, nPairs);
        clear ck;
    end

    for p = 1:nPairs
        ii = pairs_ij(p, 1);
        jj = pairs_ij(p, 2);

        if isfinite(before_mat(p, 1))
            fprintf('  [%2d/%2d]  SKIP: %s → %s\n', ...
                p, nPairs, families{ii}, families{jj});
            continue;
        end

        if ~isfinite(tickets{p}.voxelId)
            fprintf('  [%2d/%2d]  no overlap: %s → %s\n', ...
                p, nPairs, families{ii}, families{jj});
            before_mat(p, 1) = p;   % sentinel: processed, no overlap
            save(ckptFile, 'before_mat', 'after_mat', 'traj_cell', '-v7.3');
            continue;
        end

        fprintf('  [%2d/%2d]  %s → %s ...', p, nPairs, families{ii}, families{jj});
        tStart = tic;

        [before_mat(p,:), after_mat(p,:), traj_cell{p}] = ...
            local_process_pair_dc(Sall{ii}, Sall{jj}, ...
                                  families{ii}, families{jj}, ii, jj, p, ...
                                  cfg, grid3, tickets{p}, VU_mps, TU_days);

        if ~isnan(after_mat(p, strcmp(AFTER_COLS,'converged')))
            conv_val = after_mat(p, strcmp(AFTER_COLS,'converged'));
            dv_bef   = before_mat(p, strcmp(BEFORE_COLS,'DV_total_mps'));
            dv_aft   = after_mat(p,  strcmp(AFTER_COLS, 'DV_total_mps'));
            fprintf('  conv=%d  DV: %.1f→%.1f m/s  (Δ%+.1f)  (%.1fs)\n', ...
                conv_val, dv_bef, dv_aft, dv_aft-dv_bef, toc(tStart));
        else
            fprintf('  error  (%.1fs)\n', toc(tStart));
        end

        save(ckptFile, 'before_mat', 'after_mat', 'traj_cell', '-v7.3');
    end
end

% ════════════════════════════════════════════════════════════════════════════
% 5. TABULATE RESULTS + DELTA-DV SUMMARY
% ════════════════════════════════════════════════════════════════════════════
col_ef  = find(strcmp(AFTER_COLS,  'exitflag'),  1);
col_cv  = find(strcmp(AFTER_COLS,  'converged'), 1);
col_dv_bef = find(strcmp(BEFORE_COLS, 'DV_total_mps'), 1);
col_dv_aft = find(strcmp(AFTER_COLS,  'DV_total_mps'), 1);
col_prox   = find(strcmp(BEFORE_COLS, 'dA_nd'), 1);   % proxy stored via ticket

ef_col    = after_mat(:, col_ef);
cv_col    = after_mat(:, col_cv);
dv_before = before_mat(:, col_dv_bef);
dv_after  = after_mat(:, col_dv_aft);
valid     = isfinite(ef_col);
n_dc      = sum(valid);
conv_mask = valid & (cv_col == 1);

dv_change = dv_after - dv_before;  % negative = proxy overestimated, DC improved

fprintf('\n[dc_sweep] ===== RESULTS =====\n');
fprintf('  Total pairs              : %d\n', nPairs);
fprintf('  Pairs with overlap + DC  : %d\n', n_dc);
fprintf('  Pairs with no overlap    : %d\n', nPairs - n_dc);

if n_dc > 0
    n_conv  = sum(conv_mask);
    fprintf('  Converged (res<tol)      : %d  (%.1f%%)\n', n_conv, 100*n_conv/n_dc);

    if n_conv > 0
        dv_ch_conv = dv_change(conv_mask);
        thr = 1.0;   % m/s threshold for "unchanged"
        n_better  = sum(dv_ch_conv < -thr);
        n_same    = sum(abs(dv_ch_conv) <= thr);
        n_worse   = sum(dv_ch_conv > thr);

        fprintf('\n  Of the %d converged pairs:\n', n_conv);
        fprintf('    Improved  (ΔDV < -%.0f m/s)  : %d', thr, n_better);
        if n_better > 0
            fprintf('  (avg %+.1f m/s, best %+.1f m/s)', ...
                mean(dv_ch_conv(dv_ch_conv < -thr)), min(dv_ch_conv));
        end
        fprintf('\n    Unchanged (|ΔDV| ≤ %.0f m/s)  : %d\n', thr, n_same);
        fprintf('    Degraded  (ΔDV >  %.0f m/s)  : %d', thr, n_worse);
        if n_worse > 0
            fprintf('  (avg %+.1f m/s, worst %+.1f m/s)', ...
                mean(dv_ch_conv(dv_ch_conv > thr)), max(dv_ch_conv));
        end
        fprintf('\n');

        if n_worse > 0
            fprintf('\n  Pairs where DC degraded DV (sorted worst→best):\n');
            idx_worse = find(conv_mask & dv_change > thr);
            [~, srt] = sort(dv_change(idx_worse), 'descend');
            for kk = 1:min(5, numel(srt))
                pp = idx_worse(srt(kk));
                fprintf('    %-28s → %-28s  %+.1f m/s  (exitflag=%d)\n', ...
                    families{pairs_ij(pp,1)}, families{pairs_ij(pp,2)}, ...
                    dv_change(pp), ef_col(pp));
            end
        end
    end

    fprintf('\n  Exit flag breakdown:\n');
    for ef = [-2, -1, 0, 1, 2, 3]
        n_ef = sum(ef_col(valid) == ef);
        if n_ef > 0
            fprintf('    exitflag %2d : %d  (%.1f%%)\n', ef, n_ef, 100*n_ef/n_dc);
        end
    end
end

% ── Build summary table ────────────────────────────────────────────────────────
fam_A_col = families(pairs_ij(:,1))';
fam_B_col = families(pairs_ij(:,2))';
proxy_dv  = arrayfun(@(k) tickets{k}.dv_proxy_mps, 1:nPairs)';

T_summary = table(fam_A_col(:), fam_B_col(:), ...
    proxy_dv, dv_before, dv_after, dv_change, cv_col, ef_col, ...
    'VariableNames', {'FamilyA','FamilyB', ...
        'DV_proxy_mps','DV_before_mps','DV_after_mps', ...
        'DV_change_mps','converged','exitflag'});

% Sort by post-DC DV (best first; NaN at bottom)
T_summary = sortrows(T_summary, 'DV_after_mps', 'MissingPlacement', 'last');

% ════════════════════════════════════════════════════════════════════════════
% 6. SAVE RESULTS
% ════════════════════════════════════════════════════════════════════════════
resultsFile = fullfile(dcRoot, 'rs4_dc_sweep_results.mat');
save(resultsFile, ...
    'before_mat', 'after_mat', 'traj_cell', ...
    'pairs_ij',   'families',  'BEFORE_COLS', 'AFTER_COLS', ...
    'T_summary',  'cfg', ...
    '-v7.3');
fprintf('\n[dc_sweep] Saved → %s\n', resultsFile);

% Save T_summary as CSV for quick inspection
csvFile = fullfile(dcRoot, 'rs4_dc_change_summary.csv');
writetable(T_summary, csvFile);
fprintf('[dc_sweep] Summary CSV → %s\n', csvFile);

if ~use_par && exist(ckptFile, 'file')
    delete(ckptFile);
    fprintf('[dc_sweep] Checkpoint removed.\n');
end

fprintf('[dc_sweep] Done.\n');

% ════════════════════════════════════════════════════════════════════════════
% LOCAL FUNCTIONS
% ════════════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────────────
function [before_row, after_row, traj] = local_process_pair_dc( ...
        SA, SB, famA, famB, ii, jj, p, cfg, grid3, ticket, VU_mps, TU_days)
%LOCAL_PROCESS_PAIR_DC  Phase 2 DC for one pair using ticket warm start.
%
% Replaces local_process_pair: instead of calling rs4_overlap_pair +
% rs4_overlap_extract_voxel_info + rs4_overlap_visualize_bounds +
% rs4_voxel_traj_extract on the full atlas, it calls local_traj_from_ticket
% to re-integrate only the two argmin arcs identified in Phase 1.
% This eliminates all large intermediate allocation in Phase 2 workers.

before_row = NaN(1, 13);
after_row  = NaN(1, 14);
traj       = [];

try
    T  = local_traj_from_ticket(SA, SB, ticket, cfg, grid3, VU_mps, TU_days);
    if isempty(T), return; end
    Tc = rs4_diffcorr(T, SA, SB, cfg);

    % ── before row (13 cols) ─────────────────────────────────────────────────
    before_row = [ ...
        p, ii, jj, T.vid, ...
        T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps, T.DV_total_true_mps, ...
        T.dA_nd, T.dB_nd, rad2deg(T.delta_th_rad), ...
        T.tof_A_days, T.tof_B_days];

    % ── after row (14 cols) ──────────────────────────────────────────────────
    r_norm = norm(Tc.r_final);
    after_row = [ ...
        p, ii, jj, T.vid, ...
        Tc.DV_turn_A_mps, Tc.DV_patch_mps, Tc.DV_turn_B_mps, Tc.DV_total_mps, ...
        rad2deg(Tc.delta_th_rad), Tc.tof_A_days, Tc.tof_B_days, ...
        Tc.exitflag, double(Tc.converged), r_norm];

    % ── traj struct for plot replay ───────────────────────────────────────────
    traj = struct();
    traj.famA            = famA;
    traj.famB            = famB;
    traj.vid             = T.vid;
    % before DC
    traj.T_IC_A          = T.IC_A;
    traj.T_IC_B_frs      = T.IC_B_frs;
    traj.T_seed_A        = T.seed_A;
    traj.T_seed_B_frs    = T.seed_B_frs;
    traj.T_XA            = T.XA;
    traj.T_x_B           = T.x_B;
    traj.T_y_B           = T.y_B;
    traj.T_th_B          = T.th_B;
    traj.T_DV_turn_A_mps = T.DV_turn_A_mps;
    traj.T_DV_patch_mps  = T.DV_patch_mps;
    traj.T_DV_turn_B_mps = T.DV_turn_B_mps;
    traj.T_DV_total_mps  = T.DV_total_true_mps;
    traj.T_tof_A_days    = T.tof_A_days;
    traj.T_tof_B_days    = T.tof_B_days;
    traj.T_i_star        = T.i_star;
    traj.T_j_star        = T.j_star;
    traj.T_xc            = T.xc;
    traj.T_yc            = T.yc;
    traj.T_dA_nd         = T.dA_nd;
    traj.T_dB_nd         = T.dB_nd;
    traj.T_delta_th_rad  = T.delta_th_rad;
    traj.T_DV_proxy_mps  = T.DV_proxy_mps;
    % after DC
    traj.IC_A            = Tc.IC_A;
    traj.IC_B_frs        = Tc.IC_B_frs;
    traj.seed_A          = Tc.seed_A;
    traj.seed_B_frs      = Tc.seed_B_frs;
    traj.XA              = Tc.XA;
    traj.x_B             = Tc.x_B;
    traj.y_B             = Tc.y_B;
    traj.th_B            = Tc.th_B;
    traj.xp              = Tc.xp;
    traj.yp              = Tc.yp;
    traj.DV_turn_A_mps   = Tc.DV_turn_A_mps;
    traj.DV_patch_mps    = Tc.DV_patch_mps;
    traj.DV_turn_B_mps   = Tc.DV_turn_B_mps;
    traj.DV_total_mps    = Tc.DV_total_mps;
    traj.tof_A_days      = Tc.tof_A_days;
    traj.tof_B_days      = Tc.tof_B_days;
    traj.delta_th_rad    = Tc.delta_th_rad;
    traj.r_final         = Tc.r_final;
    traj.exitflag        = Tc.exitflag;
    traj.converged       = Tc.converged;
    traj.r_norm_scaled   = r_norm;

catch ME
    fprintf('    [pair %d] ERROR %s→%s: %s\n', p, famA, famB, ME.message);
end
end

% ─────────────────────────────────────────────────────────────────────────────
function F = local_compute_footprint_dc(S, grid3, VU_mps, TU_days)
%LOCAL_COMPUTE_FOOTPRINT_DC  Extended footprint with argmin row identifiers.
% Identical to the version in run_rs4_dc_sweep_tmax.m.

Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

thm    = rs3_wrapToPi(pi - grid3.th_centers(:));
lut    = discretize(thm, grid3.th_edges);
lut(isnan(lut)) = 0;
it_lut = uint16(lut);

dlists = S.Step4.delta_lists;
Ns     = numel(dlists);
max_h  = max(1, max(cellfun(@numel, dlists)));
delta_mat = zeros(Ns, max_h);
for s = 1:Ns
    v = double(dlists{s});
    delta_mat(s, 1:numel(v)) = v;
end

pot_u    = rs3_core_cr3bp_U_and_derivs(S.SeedsUpper(:,1), S.SeedsUpper(:,2), S.mu);
v0_upper = sqrt(max(2*pot_u.U(:) - S.CJ, 0));
if isfield(S, 'SeedsLower') && ~isempty(S.SeedsLower)
    pot_l    = rs3_core_cr3bp_U_and_derivs(S.SeedsLower(:,1), S.SeedsLower(:,2), S.mu);
    v0_lower = sqrt(max(2*pot_l.U(:) - S.CJ, 0));
else
    v0_lower = v0_upper;
end

nu = double(S.Step4.rows_FRS_upper.n);
if nu > 0
    [ids_u, dv_u, t_u, ix_u, iy_u, it_u, iSeed_u, iHead_u, t_nd_u] = local_fp_rows_dc( ...
        S.Step4.rows_FRS_upper, nu, v0_upper, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
    half_u = repmat(int8(1), nu, 1);
else
    ids_u=[]; dv_u=[]; t_u=[]; ix_u=[]; iy_u=[]; it_u=[];
    iSeed_u=[]; iHead_u=[]; t_nd_u=[]; half_u=int8([]);
end

nl = double(S.Step4.rows_FRS_lower.n);
if nl > 0
    [ids_l, dv_l, t_l, ix_l, iy_l, it_l, iSeed_l, iHead_l, t_nd_l] = local_fp_rows_dc( ...
        S.Step4.rows_FRS_lower, nl, v0_lower, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
    half_l = repmat(int8(-1), nl, 1);
else
    ids_l=[]; dv_l=[]; t_l=[]; ix_l=[]; iy_l=[]; it_l=[];
    iSeed_l=[]; iHead_l=[]; t_nd_l=[]; half_l=int8([]);
end

[F.uid_frs, F.dv_min_frs, F.t_mean_frs, ...
 F.iSeed_argmin_frs, F.iHead_argmin_frs, F.halfFlag_argmin_frs, F.t_nd_argmin_frs] = ...
    local_fp_agg_dc([ids_u; ids_l], [dv_u; dv_l], [t_u; t_l], ...
                    [iSeed_u; iSeed_l], [iHead_u; iHead_l], ...
                    [half_u;  half_l],  [t_nd_u;  t_nd_l]);

if ~isempty(ix_u)
    biy_u = Ny - iy_u + 1;
    bit_u = double(it_lut(it_u));
    ok_u  = bit_u > 0;
    ids_brs_u = sub2ind([Ny,Nx,Nt], biy_u(ok_u), ix_u(ok_u), max(1,min(Nt,bit_u(ok_u))));
    dv_bu=dv_u(ok_u); t_bu=t_u(ok_u); t_nd_bu=t_nd_u(ok_u);
    iS_bu=iSeed_u(ok_u); iH_bu=iHead_u(ok_u); half_bu=half_u(ok_u);
else
    ids_brs_u=[]; dv_bu=[]; t_bu=[]; t_nd_bu=[]; iS_bu=[]; iH_bu=[]; half_bu=int8([]);
end

if ~isempty(ix_l)
    biy_l = Ny - iy_l + 1;
    bit_l = double(it_lut(it_l));
    ok_l  = bit_l > 0;
    ids_brs_l = sub2ind([Ny,Nx,Nt], biy_l(ok_l), ix_l(ok_l), max(1,min(Nt,bit_l(ok_l))));
    dv_bl=dv_l(ok_l); t_bl=t_l(ok_l); t_nd_bl=t_nd_l(ok_l);
    iS_bl=iSeed_l(ok_l); iH_bl=iHead_l(ok_l); half_bl=half_l(ok_l);
else
    ids_brs_l=[]; dv_bl=[]; t_bl=[]; t_nd_bl=[]; iS_bl=[]; iH_bl=[]; half_bl=int8([]);
end

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

ix_out    = double(rows.ix(1:n));
iy_out    = double(rows.iy(1:n));
it_out    = double(rows.it(1:n));
ids       = sub2ind([Ny, Nx, Nt], iy_out, ix_out, it_out);
iSeed_out = double(rows.iSeed(1:n));
iHead_out = double(rows.iHead(1:n));
t_nd_out  = double(rows.t(1:n));
lin       = sub2ind([Ns, max_h], iSeed_out, iHead_out);
delta     = delta_mat(lin);
v0        = v0_per_seed(iSeed_out);
dv_mps    = 2 * v0(:) .* sin(abs(delta(:)) / 2) * VU_mps;
t_days    = abs(t_nd_out(:)) * TU_days;
end

% ─────────────────────────────────────────────────────────────────────────────
function [uid, dv_min, t_mean, iSeed_win, iHead_win, half_win, t_nd_win] = ...
        local_fp_agg_dc(ids, dv, t, iSeed_arr, iHead_arr, half_arr, t_nd_arr)

if isempty(ids)
    uid=zeros(0,1); dv_min=zeros(0,1); t_mean=zeros(0,1);
    iSeed_win=uint16(zeros(0,1)); iHead_win=uint16(zeros(0,1));
    half_win=int8(zeros(0,1)); t_nd_win=single(zeros(0,1));
    return;
end
[uid, ~, ic] = unique(ids(:));
[~, srt] = sortrows([ic(:), dv(:)]);
[~, first_in_group] = unique(ic(srt), 'first');
win = srt(first_in_group);
dv_min    = dv(win);
t_mean    = accumarray(ic, t(:)) ./ accumarray(ic, ones(numel(ic), 1));
iSeed_win = uint16(iSeed_arr(win));
iHead_win = uint16(iHead_arr(win));
half_win  = int8(half_arr(win));
t_nd_win  = single(t_nd_arr(win));
end

% ─────────────────────────────────────────────────────────────────────────────
function ticket = local_run_pair_dc(FA, FB, grid3, cfg, VU_mps)

ticket = struct('voxelId',NaN, 'dv_proxy_mps',NaN, 'dvlb_mps',NaN, ...
                'dvpatch_mps',NaN, 'tof_days',NaN, ...
                'iSeed_A',0, 'iHead_A',0, 'halfFlag_A',int8(0), 't_nd_A',single(0), ...
                'iSeed_B',0, 'iHead_B',0, 'halfFlag_B',int8(0), 't_nd_B',single(0));
try
    idsO = intersect(FA.uid_frs, FB.uid_brs);
    if isempty(idsO), return; end

    Ny = numel(grid3.y_centers);
    Nx = numel(grid3.x_centers);
    Nt = numel(grid3.th_centers);
    [iy, ix, ~] = ind2sub([Ny, Nx, Nt], idsO);

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
    ok = okKeep(:) & ...
         (hypot(x + mu, y)     > (1+bufFrac)*RE) & ...
         (hypot(x - (1-mu), y) > (1+bufFrac)*RM);

    idsO = idsO(ok);
    if isempty(idsO), return; end
    ix = ix(ok);  iy = iy(ok);

    [~, locA] = ismember(idsO, FA.uid_frs);
    [~, locB] = ismember(idsO, FB.uid_brs);
    x_ok = grid3.x_centers(ix);
    y_ok = grid3.y_centers(iy);
    CJstar = min(FA.CJ, FB.CJ);
    pot    = rs3_core_cr3bp_U_and_derivs(x_ok(:), y_ok(:), mu);
    v_box  = sqrt(max(2*pot.U - CJstar, 0));
    dv_patch_vec = 2*v_box .* sin(abs(grid3.dtheta)/2) * VU_mps;
    dv_lb_vec    = FA.dv_min_frs(locA(:)) + FB.dv_min_brs(locB(:));
    dv_proxy     = dv_lb_vec + dv_patch_vec;

    valid = isfinite(dv_proxy);
    if ~any(valid), return; end
    idxValid = find(valid);
    [~, iLoc] = min(dv_proxy(idxValid));
    iWin = idxValid(iLoc);

    ticket.voxelId      = idsO(iWin);
    ticket.dv_proxy_mps = dv_proxy(iWin);
    ticket.dvlb_mps     = dv_lb_vec(iWin);
    ticket.dvpatch_mps  = dv_patch_vec(iWin);
    ticket.tof_days     = FA.t_mean_frs(locA(iWin)) + FB.t_mean_brs(locB(iWin));

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
function T = local_traj_from_ticket(SA, SB, ticket, cfg, grid3, VU_mps, TU_days)
%LOCAL_TRAJ_FROM_TICKET  Reconstruct T struct for rs4_diffcorr from ticket.
% Identical to the version in run_rs4_dc_sweep_tmax.m.

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

% ── Side A ────────────────────────────────────────────────────────────────────
iSeed_A = ticket.iSeed_A;  iHead_A = ticket.iHead_A;
halfA   = ticket.halfFlag_A; t_nd_A = ticket.t_nd_A;

if halfA >= 0
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
T.DV_turn_A_mps = 2*v0_A*sin(abs(delta_A)/2)*VU_mps;
T.IC_A=IC_A; T.seed_A=seed_A; T.t_A=t_end_A; T.tA_vec=tA_span; T.XA=XA;
T.iSeed_A=iSeed_A; T.iHead_A=iHead_A; T.halfFlag_A=halfA;

% ── Side B ────────────────────────────────────────────────────────────────────
iSeed_B = ticket.iSeed_B;  iHead_B = ticket.iHead_B;
halfB   = ticket.halfFlag_B; t_nd_B = ticket.t_nd_B;

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

x_B  = XB_frs(:,1);
y_B  = -XB_frs(:,2);
th_B = rs3_wrapToPi(pi - XB_frs(:,3));

pot_B = rs3_core_cr3bp_U_and_derivs(seed_B_frs(1), seed_B_frs(2), SB.mu);
v0_B  = sqrt(max(2*pot_B.U - SB.CJ, 0));
T.DV_turn_B_mps=2*v0_B*sin(abs(delta_B)/2)*VU_mps;
T.IC_B_frs=IC_B_frs; T.seed_B_frs=seed_B_frs; T.t_B=t_end_B; T.tB_vec=tB_span;
T.x_B=x_B; T.y_B=y_B; T.th_B=th_B;
T.iSeed_B=iSeed_B; T.iHead_B=iHead_B; T.halfFlag_B_frs=halfB; T.from_lower_B=(halfB<0);

% ── Closest approach ─────────────────────────────────────────────────────────
xc=T.xc; yc=T.yc;
dA_vec=hypot(XA(:,1)-xc, XA(:,2)-yc); [T.dA_nd,T.i_star]=min(dA_vec);
T.th_A_star=XA(T.i_star,3);
dB_vec=hypot(x_B-xc, y_B-yc); [T.dB_nd,T.j_star]=min(dB_vec);
T.th_B_star=th_B(T.j_star);
T.delta_th_rad=abs(rs3_circ_diff(T.th_A_star, T.th_B_star));

% ── DV metrics ───────────────────────────────────────────────────────────────
CJstar=min(SA.CJ,SB.CJ);
pot_c=rs3_core_cr3bp_U_and_derivs(xc,yc,SA.mu);
T.v_box_center_nd=sqrt(max(2*pot_c.U-CJstar,0));
T.DV_patch_nd=2*T.v_box_center_nd*sin(T.delta_th_rad/2);
T.DV_patch_mps=T.DV_patch_nd*VU_mps;
T.DV_total_true_mps=T.DV_turn_A_mps+T.DV_patch_mps+T.DV_turn_B_mps;
T.DV_proxy_mps=ticket.dv_proxy_mps;
T.tof_A_days=abs(t_nd_A)*TU_days;
T.tof_B_days=abs(t_nd_B)*TU_days;
end

% ─────────────────────────────────────────────────────────────────────────────
function S = local_ensure_xpo(S, relTol, absTol, N_po)
if isfield(S,'Xpo') && ~isempty(S.Xpo) && isfield(S,'t_dense') && ~isempty(S.t_dense)
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
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
