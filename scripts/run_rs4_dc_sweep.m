%% RUN_RS4_DC_SWEEP
% Self-contained batch differential correction over every pair of the 13
% periodic-orbit families.
%
% Pipeline per pair (identical to run_rs4_voxel_trajectories.m):
%   rs4_overlap_pair               -> O
%   rs4_overlap_extract_voxel_info -> V
%   rs4_overlap_visualize_bounds   -> B  (no plots)
%   rs4_voxel_traj_extract         -> T  (before-DC arcs + true DV)
%   rs4_diffcorr                   -> Tc (after-DC arcs)
%
% Two modes (set NUM_WORKERS below):
%
%   SERIAL  (NUM_WORKERS <= 1)
%     All 13 families loaded into memory once.  Passed directly to the for
%     loop — no file I/O per pair.  Checkpoint saved after every pair so a
%     run can be safely interrupted and resumed by re-running the script.
%
%   PARALLEL (NUM_WORKERS >= 2)
%     Uses a PROCESS-BASED pool (parpool 'Processes') so workers have full
%     MATLAB: Java is available, V7.3 .mat files load fine.
%     Each worker loads only its 2 needed families from disk — no large
%     struct broadcast, no OOM.
%     Cache paths are pre-computed in the serial warm-up phase (where Java
%     is available) and passed to workers as plain strings.
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
NUM_WORKERS = 0;   % 0 or 1 = serial with checkpoint
                   % >= 2    = process-based parfor (each worker loads 2 families)

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
cfg.plot.rs4.bounds_lb    = false;
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = false;

% -------------------------------------------------------------------------
% CACHE SETTINGS — must match exactly what was used to build the caches
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

N      = numel(families);
relTol = cfg.propag.relTol;
absTol = cfg.propag.absTol;

% ---- enumerate pairs early (needed to set use_par) -----------------------
pairs_ij = zeros(0, 2);
for ii = 1:N
    for jj = ii+1:N
        pairs_ij(end+1, :) = [ii, jj]; %#ok<AGROW>
    end
end
nPairs = size(pairs_ij, 1);

use_par = NUM_WORKERS >= 2 && ...
          license('test', 'Distrib_Computing_Toolbox') && ...
          nPairs > 1;

% ---- output folder -------------------------------------------------------
dcRoot = fullfile(cfg.io.out_root, cfg.io.tag, 'rs4_dc_sweep');
if ~exist(dcRoot, 'dir'), mkdir(dcRoot); end
fprintf('[dc_sweep] Output   : %s\n', dcRoot);
fprintf('[dc_sweep] Mode     : %s\n\n', ...
    ternary(use_par, sprintf('PARALLEL (%d workers)', min(NUM_WORKERS,nPairs)), 'SERIAL'));

% =========================================================================
% 1. WARM CACHE + GET FILE PATHS  (serial, main process)
%    rs3_prepare_or_load_family uses Java for MD5 — must run here.
%    Captures the exact .mat path for each family.
%
%    SERIAL MODE:  keep the loaded struct in Sall{k}  (in-memory pass-through)
%    PARALLEL MODE: discard the struct — workers will load per-pair from disk
% =========================================================================
grid3 = rs3_grid_make(cfg);

fprintf('[dc_sweep] Warming %d family caches ...\n', N);
cache_fpaths = cell(N, 1);
Sall         = cell(N, 1);   % used only in serial mode

for k = 1:N
    fprintf('  [%2d/%2d]  %s ...', k, N, families{k});
    tL = tic;
    [Stmp, ~] = rs3_prepare_or_load_family(families{k}, cfg, grid3);
    info = rs3_cache_get_path(families{k}, Stmp.mu, Stmp.CJ, cfg);
    cache_fpaths{k} = info.fpath;
    if ~use_par
        Sall{k} = local_ensure_xpo(Stmp, relTol, absTol, 1001);
    end
    fprintf('  %.1fs  Nseeds_upper=%d\n', toc(tL), size(Stmp.SeedsUpper, 1));
    clear Stmp info;
end

if use_par
    clear Sall;   % not needed — keep memory free before launching workers
    fprintf('[dc_sweep] Cache warmed. Family structs released (workers load per-pair).\n\n');
else
    fprintf('[dc_sweep] Families in memory.\n\n');
end

fprintf('[dc_sweep] %d pairs to process.\n\n', nPairs);

% =========================================================================
% 2. PRE-ALLOCATE OUTPUT
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
% 3. MAIN LOOP
% =========================================================================
ckptFile = fullfile(dcRoot, 'rs4_dc_sweep_checkpoint.mat');

if use_par
    % ------------------------------------------------------------------
    % PARALLEL — process-based pool.
    %   Broadcast: only small data (path strings, cfg, family names).
    %   Each worker loads ONLY its 2 families from disk.
    %   Process-based workers support Java + V7.3 MAT loading.
    % ------------------------------------------------------------------
    nW = min(NUM_WORKERS, nPairs);

    % Create a fresh process-based pool
    delete(gcp('nocreate'));
    parpool('Processes', nW);

    % tiny broadcasts
    fpaths   = cache_fpaths;
    cfg_p    = cfg;
    dcR      = dcRoot;
    fams     = families;
    pij      = pairs_ij;
    rT       = relTol;
    aT       = absTol;

    parfor p = 1:nPairs
        ii = pij(p, 1);
        jj = pij(p, 2);

        % Each worker loads only 2 families (no OOM from broadcasting)
        % Process-based workers can load V7.3 and use Java fine.
        dA = load(fpaths{ii}, 'S');
        SA = local_ensure_xpo(dA.S, rT, aT, 1001);
        clear dA;

        dB = load(fpaths{jj}, 'S');
        SB = local_ensure_xpo(dB.S, rT, aT, 1001);
        clear dB;

        [br, ar, tr] = local_process_pair( ...
            SA, SB, fams{ii}, fams{jj}, ii, jj, p, cfg_p, dcR);

        before_mat(p, :) = br; %#ok<PFPIE>
        after_mat(p, :)  = ar; %#ok<PFPIE>
        traj_cell{p}     = tr; %#ok<PFPIE>
    end

else
    % ------------------------------------------------------------------
    % SERIAL — families already in memory, checkpoint after every pair.
    % ------------------------------------------------------------------
    fprintf('[dc_sweep] Checkpoint: %s\n\n', ckptFile);

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

        if isfinite(before_mat(p, 1))
            fprintf('  [%2d/%2d]  SKIP: %s -> %s\n', ...
                p, nPairs, families{ii}, families{jj});
            continue;
        end

        fprintf('  [%2d/%2d]  %s -> %s ...', p, nPairs, families{ii}, families{jj});
        tStart = tic;

        % Families already in memory — no file I/O here
        [before_mat(p,:), after_mat(p,:), traj_cell{p}] = ...
            local_process_pair(Sall{ii}, Sall{jj}, ...
                               families{ii}, families{jj}, ii, jj, p, ...
                               cfg, dcRoot);

        if ~isfinite(before_mat(p, 1))
            before_mat(p, 1) = p;   % sentinel for no-overlap pairs
        end

        tr = traj_cell{p};
        if ~isempty(tr)
            fprintf('  conv=%d  DV: %.1f->%.1f m/s  (%.1fs)\n', ...
                tr.converged, tr.T_DV_total_mps, tr.DV_total_mps, toc(tStart));
        else
            fprintf('  no overlap  (%.1fs)\n', toc(tStart));
        end

        save(ckptFile, 'before_mat', 'after_mat', 'traj_cell', '-v7.3');
    end
end

% =========================================================================
% 4. TABULATE RESULTS
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
% 5. SAVE RESULTS
% =========================================================================
resultsFile = fullfile(dcRoot, 'rs4_dc_sweep_results.mat');
save(resultsFile, ...
    'before_mat', 'after_mat', 'traj_cell', ...
    'pairs_ij',   'families',  'BEFORE_COLS', 'AFTER_COLS', 'cfg', ...
    '-v7.3');
fprintf('\n[dc_sweep] Saved -> %s\n', resultsFile);

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
%   SA, SB already in memory (pre-loaded by caller or by worker).
%   Mirrors run_rs4_voxel_trajectories.m lines 66-112 exactly.

before_row = NaN(1, 13);
after_row  = NaN(1, 14);
traj       = [];

try
    O = rs4_overlap_pair(SA, SB, cfg);
    if isempty(O.ids), return; end

    V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg);

    pairTag = sprintf('pair_%02d_%02d', ii, jj);
    B = rs4_overlap_visualize_bounds(V, SA, SB, O, cfg, dcRoot, pairTag);

    if ~isfield(B,'imin') || ~isfinite(B.imin) || ...
            isempty(B.dv_proxy) || ~any(isfinite(B.dv_proxy))
        return;
    end

    T  = rs4_voxel_traj_extract(SA, SB, V, B, cfg);
    Tc = rs4_diffcorr(T, SA, SB, cfg);

    % --- before row (13 cols) ---
    before_row = [ ...
        p, ii, jj, T.vid, ...
        T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps, T.DV_total_true_mps, ...
        T.dA_nd, T.dB_nd, rad2deg(T.delta_th_rad), ...
        T.tof_A_days, T.tof_B_days];

    % --- after row (14 cols) ---
    r_norm = norm(Tc.r_final);
    after_row = [ ...
        p, ii, jj, T.vid, ...
        Tc.DV_turn_A_mps, Tc.DV_patch_mps, Tc.DV_turn_B_mps, Tc.DV_total_mps, ...
        rad2deg(Tc.delta_th_rad), Tc.tof_A_days, Tc.tof_B_days, ...
        Tc.exitflag, double(Tc.converged), r_norm];

    % --- traj struct for plot replay ---
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
    fprintf('    [pair %d] ERROR %s->%s: %s\n', p, famA, famB, ME.message);
end
end

% -------------------------------------------------------------------------

function S = local_ensure_xpo(S, relTol, absTol, N_po)
% Re-integrate the periodic orbit if Xpo was stripped from cache.
    if isfield(S,'Xpo') && ~isempty(S.Xpo) && ...
       isfield(S,'t_dense') && ~isempty(S.t_dense)
        return;
    end
    opts    = odeset('RelTol', relTol, 'AbsTol', absTol);
    solPO   = ode113(@(t,X) rs3_core_reduced_cr3bp_model(t,X,S.CJ,S.mu,true), ...
                     [0, S.Tf_PO], S.X0, opts);
    t_dense = linspace(0, S.Tf_PO, N_po)';
    S.t_dense = t_dense;
    S.Xpo     = deval(solPO, t_dense)';
end

% -------------------------------------------------------------------------

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
