%% RUN_RS4_DC_SWEEP
% Self-contained batch differential correction over every pair of the 13
% periodic-orbit families.
%
% Does NOT rely on cached pair_result.mat files from
% run_rs4_all_pairs_summary.m.  Computes everything from scratch:
%
%   For each pair (i, j):
%     1. rs4_overlap_pair               -> O  (overlap voxel IDs)
%     2. rs4_overlap_extract_voxel_info -> V  (per-voxel DV metadata)
%     3. rs4_overlap_visualize_bounds   -> B  (winner voxel, no plots)
%     4. rs4_voxel_traj_extract         -> T  (before-DC arcs + true DV)
%     5. rs4_diffcorr                   -> Tc (after-DC arcs)
%
% Parallelism (USE_PARFOR = true):
%   Each parfor worker loads its two family structs from cache
%   independently — the heavy Step4.rows_FRS_* arrays are NOT broadcast.
%   Requires Parallel Computing Toolbox.
%
% Serial mode (USE_PARFOR = false):
%   Saves a checkpoint after every pair so a run can be safely interrupted
%   and resumed by re-running the script.
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
USE_PARFOR = true;    % false = serial with per-pair checkpoint saves

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
cfg.cache.rebuild         = false;

% No figures anywhere in this runner
cfg.io.save_figs          = false;
cfg.io.save_fig           = false;
cfg.io.fig_visible        = 'off';
cfg.diffcorr.display      = 'off';   % suppress fmincon iteration output
cfg.plot.rs4.bounds_lb    = false;   % no plots inside visualize_bounds
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = false;

% Grid / propagation settings (must match the family caches)
cfg.grid.dx               = 0.01;
cfg.grid.dy               = 0.01;
cfg.grid.dtheta           = deg2rad(2);
cfg.seed.ds_seed          = 0.02;
cfg.propag.Tmax           = pi;
cfg.fan.DV_cap_nd         = 0.2;
cfg.fan.dtheta_fan        = deg2rad(1.0);
cfg.propag.absTol         = 1e-9;
cfg.propag.relTol         = 1e-9;

N = numel(families);

% ---- output folder -------------------------------------------------------
dcRoot = fullfile(cfg.io.out_root, cfg.io.tag, 'rs4_dc_sweep');
if ~exist(dcRoot, 'dir'), mkdir(dcRoot); end
fprintf('[dc_sweep] Output : %s\n\n', dcRoot);

% =========================================================================
% 1. VERIFY / WARM CACHE
%    Load each family once (serial) so cache files definitely exist before
%    parallel workers try to read them concurrently.
% =========================================================================
grid3 = rs3_grid_make(cfg);
fprintf('[dc_sweep] Verifying %d family caches ...\n', N);
for k = 1:N
    [Stmp, ~] = rs3_prepare_or_load_family(families{k}, cfg, grid3);
    fprintf('  [%2d/%2d]  %-30s  Nseeds_upper=%d\n', k, N, ...
        families{k}, size(Stmp.SeedsUpper, 1));
    clear Stmp;
end
fprintf('[dc_sweep] Cache verified.\n\n');

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

use_par = USE_PARFOR && ...
          license('test', 'Distrib_Computing_Toolbox') && ...
          nPairs > 1;

if use_par
    % ------------------------------------------------------------------
    % PARALLEL — each worker loads its two families from cache.
    % Only tiny things are broadcast: cfg, grid3, families, dcRoot.
    % ------------------------------------------------------------------
    fprintf('[dc_sweep] Mode: PARALLEL (parfor).\n');
    fprintf('           Each worker loads family caches independently.\n\n');

    fam_cell   = families;   % cell of strings — tiny
    cfg_par    = cfg;
    grid3_par  = grid3;
    dcRoot_par = dcRoot;

    parfor p = 1:nPairs
        ii = pairs_ij(p, 1);
        jj = pairs_ij(p, 2);
        [br, ar, tr] = local_process_pair( ...
            fam_cell{ii}, fam_cell{jj}, ii, jj, p, cfg_par, grid3_par, dcRoot_par);
        before_mat(p, :) = br; %#ok<PFPIE>
        after_mat(p, :)  = ar; %#ok<PFPIE>
        traj_cell{p}     = tr; %#ok<PFPIE>
    end

else
    % ------------------------------------------------------------------
    % SERIAL — checkpoint after every pair for safe resume.
    % Re-run the script to continue from where it stopped.
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
            fprintf('  [%2d/%2d]  SKIP: %s -> %s\n', ...
                p, nPairs, families{ii}, families{jj});
            continue;
        end

        fprintf('  [%2d/%2d]  %s -> %s\n', p, nPairs, families{ii}, families{jj});
        tStart = tic;

        [before_mat(p,:), after_mat(p,:), traj_cell{p}] = ...
            local_process_pair(families{ii}, families{jj}, ii, jj, p, ...
                               cfg, grid3, dcRoot);

        % Write pair_idx sentinel even for no-overlap pairs so they are
        % not retried on resume (they are fast but still wasteful).
        if ~isfinite(before_mat(p, 1))
            before_mat(p, 1) = p;
        end

        % Progress summary
        tr = traj_cell{p};
        if ~isempty(tr)
            fprintf('    -> conv=%d  DV: %.1f -> %.1f m/s  (%.1fs)\n', ...
                tr.converged, tr.T_DV_total_mps, tr.DV_total_mps, toc(tStart));
        else
            fprintf('    -> no overlap / no winner  (%.1fs)\n', toc(tStart));
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

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [before_row, after_row, traj] = local_process_pair( ...
        famA, famB, ii, jj, p, cfg, grid3, dcRoot)
% Full pipeline for one pair.
% Loads family structs from cache here (not passed in) so this function
% can safely run inside a parfor worker without broadcasting the heavy
% Step4.rows_FRS_* arrays.

before_row = NaN(1, 13);
after_row  = NaN(1, 14);
traj       = [];

try
    % ---- Load family structs from cache ----------------------------------
    [SA] = rs3_prepare_or_load_family(famA, cfg, grid3);
    SA   = local_ensure_xpo(SA, cfg.propag.relTol, cfg.propag.absTol, 1001);
    [SB] = rs3_prepare_or_load_family(famB, cfg, grid3);
    SB   = local_ensure_xpo(SB, cfg.propag.relTol, cfg.propag.absTol, 1001);

    % ---- Overlap ---------------------------------------------------------
    O = rs4_overlap_pair(SA, SB, cfg);
    if isempty(O.ids)
        return;
    end

    % ---- Voxel metadata --------------------------------------------------
    V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg);

    % ---- Winner voxel (B) — correct signature: (V, SA, SB, O, cfg, ...) --
    pairTag = sprintf('pair_%02d_%02d', ii, jj);
    B = rs4_overlap_visualize_bounds(V, SA, SB, O, cfg, dcRoot, pairTag);

    if ~isfield(B, 'imin') || ~isfinite(B.imin) || ...
            isempty(B.dv_proxy) || ~any(isfinite(B.dv_proxy))
        return;
    end

    % ---- True trajectory -------------------------------------------------
    T  = rs4_voxel_traj_extract(SA, SB, V, B, cfg);

    % ---- Differential correction -----------------------------------------
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
    opts    = odeset('RelTol', relTol, 'AbsTol', absTol);
    solPO   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,S.CJ,S.mu,true), ...
                     [0, S.Tf_PO], S.X0, opts);
    t_dense = linspace(0, S.Tf_PO, N_po)';
    S.t_dense = t_dense;
    S.Xpo     = deval(solPO, t_dense)';
end
