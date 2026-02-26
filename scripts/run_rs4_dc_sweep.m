%% RUN_RS4_DC_SWEEP
% Self-contained batch differential correction over every pair of the 13
% periodic-orbit families.
%
% This script follows the exact same pipeline as run_rs4_voxel_trajectories.m
% but loops over all N*(N-1)/2 pairs:
%
%   1. Pre-load ALL families into memory (main process, serial)
%   2. For each pair (i, j) — serial or parfor:
%        rs4_overlap_pair               -> O
%        rs4_overlap_extract_voxel_info -> V
%        rs4_overlap_visualize_bounds   -> B  (no plots)
%        rs4_voxel_traj_extract         -> T  (before-DC arcs + true DV)
%        rs4_diffcorr                   -> Tc (after-DC arcs)
%
% All family structs live in memory — NO file I/O inside the loop.
%
% Outputs  (saved to rs3_results/<tag>/rs4_dc_sweep/rs4_dc_sweep_results.mat)
%   before_mat   [nPairs x 13]  pre-DC metrics
%   after_mat    [nPairs x 14]  post-DC metrics + solver metadata
%   traj_cell    {nPairs x 1}   struct per pair — sufficient for plot replay
%   pairs_ij     [nPairs x 2]   family index pairs (i < j)
%   families     {13 x 1}       family name strings
%   BEFORE_COLS / AFTER_COLS    column-name cell arrays
%   cfg                         config used for this run

clear; clc;

% ---- path setup ----------------------------------------------------------
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% =========================================================================
% USER SETTINGS
% =========================================================================
NUM_WORKERS = 0;   % 0 or 1 = serial (with checkpoint), >=2 = parallel
                   % parallel uses process-based pool (not threads)

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

% No figures anywhere in this runner
cfg.io.save_figs          = false;
cfg.io.save_fig           = false;
cfg.io.fig_visible        = 'off';
cfg.diffcorr.display      = 'off';   % suppress fmincon iteration output
cfg.plot.rs4.bounds_lb    = false;   % no plots inside visualize_bounds
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = false;

% -------------------------------------------------------------------------
% CACHE SETTINGS — MUST MATCH EXACTLY WHAT WAS USED TO BUILD THE FAMILY
% CACHES.  If these differ from the cache build, rs3_prepare_or_load_family
% will reject the cache and recompute the full FRS/BRS from scratch.
% -------------------------------------------------------------------------
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

N = numel(families);

% ---- output folder -------------------------------------------------------
dcRoot = fullfile(cfg.io.out_root, cfg.io.tag, 'rs4_dc_sweep');
if ~exist(dcRoot, 'dir'), mkdir(dcRoot); end
fprintf('[dc_sweep] Output : %s\n\n', dcRoot);

% =========================================================================
% 1. BUILD GRID + PRE-LOAD ALL FAMILIES (serial, main process)
%    This is the same pattern as run_rs4_voxel_trajectories.m lines 60-63.
%    Everything happens in the main process where Java and V7.3 work fine.
% =========================================================================
grid3 = rs3_grid_make(cfg);

fprintf('[dc_sweep] Loading %d families into memory ...\n', N);
Sall = cell(N, 1);
for k = 1:N
    fprintf('  [%2d/%2d]  %s ...', k, N, families{k});
    tLoad = tic;
    [Sall{k}, ~] = rs3_prepare_or_load_family(families{k}, cfg, grid3);
    Sall{k} = local_ensure_xpo(Sall{k}, cfg.propag.relTol, cfg.propag.absTol, 1001);
    fprintf('  done (%.1fs)  Nseeds_upper=%d\n', toc(tLoad), size(Sall{k}.SeedsUpper, 1));
end
fprintf('[dc_sweep] All families loaded.\n\n');

% =========================================================================
% 2. ENUMERATE ALL N*(N-1)/2 PAIRS
% =========================================================================
pairs_ij = zeros(0, 2);
for ii = 1:N
    for jj = ii+1:N
        pairs_ij(end+1, :) = [ii, jj]; %#ok<AGROW>
    end
end
nPairs = size(pairs_ij, 1);
fprintf('[dc_sweep] %d pairs to process  (%d families).\n\n', nPairs, N);

% =========================================================================
% 3. PRE-ALLOCATE OUTPUT
% =========================================================================
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

nBC = numel(BEFORE_COLS);   % 13
nAC = numel(AFTER_COLS);    % 14

before_mat = NaN(nPairs, nBC);
after_mat  = NaN(nPairs, nAC);
traj_cell  = cell(nPairs, 1);

% =========================================================================
% 4. MAIN LOOP
% =========================================================================
ckptFile = fullfile(dcRoot, 'rs4_dc_sweep_checkpoint.mat');

use_par = NUM_WORKERS >= 2 && ...
          license('test', 'Distrib_Computing_Toolbox') && ...
          nPairs > 1;

if use_par
    % ------------------------------------------------------------------
    % PARALLEL — process-based pool (not threads!)
    % Family structs are broadcast from main process — no file I/O in
    % workers.  Process-based pools support Java + V7.3 if ever needed.
    % ------------------------------------------------------------------
    nW = min(NUM_WORKERS, nPairs);
    fprintf('[dc_sweep] Mode: PARALLEL (%d process-based workers).\n', nW);
    fprintf('           Family structs are in memory — no file I/O in workers.\n\n');

    % Ensure we have a process-based pool of the right size
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= nW || ...
            ~strcmp(pool.Cluster.Profile, 'Processes')
        delete(gcp('nocreate'));
        parpool('Processes', nW);
    end

    % Broadcast variables (read-only inside parfor)
    Sall_par = Sall;
    cfg_par  = cfg;
    dcRoot_par = dcRoot;
    fam_cell = families;
    pij = pairs_ij;

    parfor p = 1:nPairs
        ii = pij(p, 1);
        jj = pij(p, 2);
        [br, ar, tr] = local_process_pair( ...
            Sall_par{ii}, Sall_par{jj}, ...
            fam_cell{ii}, fam_cell{jj}, ii, jj, p, ...
            cfg_par, dcRoot_par);
        before_mat(p, :) = br; %#ok<PFPIE>
        after_mat(p, :)  = ar; %#ok<PFPIE>
        traj_cell{p}     = tr; %#ok<PFPIE>
    end

else
    % ------------------------------------------------------------------
    % SERIAL — checkpoint after every pair for safe resume.
    % ------------------------------------------------------------------
    fprintf('[dc_sweep] Mode: SERIAL (checkpoint: %s).\n\n', ckptFile);

    % Resume from checkpoint if available
    if exist(ckptFile, 'file')
        ck = load(ckptFile, 'before_mat', 'after_mat', 'traj_cell');
        before_mat = ck.before_mat;
        after_mat  = ck.after_mat;
        traj_cell  = ck.traj_cell;
        n_done = sum(isfinite(before_mat(:, 1)));
        fprintf('[dc_sweep] Resumed from checkpoint (%d / %d done).\n\n', ...
            n_done, nPairs);
        clear ck;
    end

    for p = 1:nPairs
        ii = pairs_ij(p, 1);
        jj = pairs_ij(p, 2);

        % col 1 = pair_idx — finite means this pair was already tried
        if isfinite(before_mat(p, 1))
            fprintf('  [%2d/%2d]  SKIP (done): %s -> %s\n', ...
                p, nPairs, families{ii}, families{jj});
            continue;
        end

        fprintf('  [%2d/%2d]  %s -> %s ...', p, nPairs, families{ii}, families{jj});
        tStart = tic;

        [before_mat(p,:), after_mat(p,:), traj_cell{p}] = ...
            local_process_pair(Sall{ii}, Sall{jj}, ...
                               families{ii}, families{jj}, ii, jj, p, ...
                               cfg, dcRoot);

        % Write pair_idx sentinel even for no-overlap pairs so they are
        % not retried on resume.
        if ~isfinite(before_mat(p, 1))
            before_mat(p, 1) = p;
        end

        % Progress summary
        tr = traj_cell{p};
        if ~isempty(tr)
            fprintf('  conv=%d  DV: %.1f -> %.1f m/s  (%.1fs)\n', ...
                tr.converged, tr.T_DV_total_mps, tr.DV_total_mps, toc(tStart));
        else
            fprintf('  no overlap / no winner  (%.1fs)\n', toc(tStart));
        end

        % Checkpoint
        save(ckptFile, 'before_mat', 'after_mat', 'traj_cell', '-v7.3');
    end
end

% =========================================================================
% 5. TABULATE RESULTS
% =========================================================================
ac_ef = find(strcmp(AFTER_COLS, 'exitflag'),  1);
ac_cv = find(strcmp(AFTER_COLS, 'converged'), 1);

ef_col = after_mat(:, ac_ef);
cv_col = after_mat(:, ac_cv);
valid  = isfinite(ef_col);
n_dc   = sum(valid);

fprintf('\n[dc_sweep] ===== RESULTS =====\n');
fprintf('  Total pairs            : %d\n', nPairs);
fprintf('  Pairs with overlap+DC  : %d\n', n_dc);
fprintf('  Pairs with no overlap  : %d\n', nPairs - n_dc);
if n_dc > 0
    fprintf('  Converged (res<tol)    : %d  (%.1f%%)\n', ...
        sum(cv_col(valid) == 1), 100 * mean(cv_col(valid) == 1));
    fprintf('\n  Exit flag breakdown:\n');
    for ef = [-2, -1, 0, 1, 2, 3]
        n_ef = sum(ef_col(valid) == ef);
        if n_ef > 0
            fprintf('    exitflag %2d : %d  (%.1f%%)\n', ef, n_ef, 100*n_ef/n_dc);
        end
    end
end

% =========================================================================
% 6. SAVE RESULTS
% =========================================================================
resultsFile = fullfile(dcRoot, 'rs4_dc_sweep_results.mat');
save(resultsFile, ...
    'before_mat', 'after_mat', 'traj_cell', ...
    'pairs_ij',   'families',  'BEFORE_COLS', 'AFTER_COLS', 'cfg', ...
    '-v7.3');
fprintf('\n[dc_sweep] Saved -> %s\n', resultsFile);

% Remove checkpoint on clean serial finish
if ~use_par && exist(ckptFile, 'file')
    delete(ckptFile);
    fprintf('[dc_sweep] Checkpoint removed.\n');
end

fprintf('[dc_sweep] Done.\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [before_row, after_row, traj] = local_process_pair( ...
        SA, SB, famA, famB, ii, jj, p, cfg, dcRoot)
%LOCAL_PROCESS_PAIR  Full pipeline for one pair.
%   Receives family structs SA, SB directly from the caller (in-memory).
%   NO file I/O — no load(), no rs3_md5, no Java calls.
%   Exactly mirrors run_rs4_voxel_trajectories.m lines 66-112.

before_row = NaN(1, 13);
after_row  = NaN(1, 14);
traj       = [];

try
    % ---- Overlap (same as run_rs4_voxel_trajectories.m line 66) ----------
    O = rs4_overlap_pair(SA, SB, cfg);
    if isempty(O.ids)
        return;
    end

    % ---- Voxel metadata (line 67) ----------------------------------------
    V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg);

    % ---- Winner voxel B (line 70) — CORRECT signature: V, SA, SB, O, ... -
    pairTag = sprintf('pair_%02d_%02d', ii, jj);
    B = rs4_overlap_visualize_bounds(V, SA, SB, O, cfg, dcRoot, pairTag);

    if ~isfield(B, 'imin') || ~isfinite(B.imin) || ...
            isempty(B.dv_proxy) || ~any(isfinite(B.dv_proxy))
        return;
    end

    % ---- True trajectory (line 76) ---------------------------------------
    T = rs4_voxel_traj_extract(SA, SB, V, B, cfg);

    % ---- Differential correction (line 97) -------------------------------
    Tc = rs4_diffcorr(T, SA, SB, cfg);

    % ---- Pack before row (13 cols) ---------------------------------------
    before_row = [ ...
        p, ii, jj, T.vid, ...
        T.DV_turn_A_mps,  T.DV_patch_mps,  T.DV_turn_B_mps,  T.DV_total_true_mps, ...
        T.dA_nd, T.dB_nd, rad2deg(T.delta_th_rad), ...
        T.tof_A_days, T.tof_B_days];

    % ---- Pack after row (14 cols) ----------------------------------------
    r_norm = norm(Tc.r_final);
    after_row = [ ...
        p, ii, jj, T.vid, ...
        Tc.DV_turn_A_mps, Tc.DV_patch_mps, Tc.DV_turn_B_mps, Tc.DV_total_mps, ...
        rad2deg(Tc.delta_th_rad), Tc.tof_A_days, Tc.tof_B_days, ...
        Tc.exitflag, double(Tc.converged), r_norm];

    % ---- Trajectory struct for plot replay -------------------------------
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
    fprintf('    [pair %d] ERROR %s -> %s:\n      %s\n', ...
        p, famA, famB, ME.message);
end
end

% -------------------------------------------------------------------------

function S = local_ensure_xpo(S, relTol, absTol, N_po)
% Re-integrate the periodic orbit if Xpo was stripped from cache.
    if isfield(S, 'Xpo') && ~isempty(S.Xpo) && ...
       isfield(S, 't_dense') && ~isempty(S.t_dense)
        return;
    end
    fprintf('    [ensure_xpo] rebuilding Xpo for "%s" ...\n', S.name);
    opts    = odeset('RelTol', relTol, 'AbsTol', absTol);
    solPO   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,S.CJ,S.mu,true), ...
                     [0, S.Tf_PO], S.X0, opts);
    t_dense = linspace(0, S.Tf_PO, N_po)';
    S.t_dense = t_dense;
    S.Xpo     = deval(solPO, t_dense)';
end
