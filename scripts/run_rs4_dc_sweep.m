%% RUN_RS4_DC_SWEEP
% Batch differential correction over every valid pair of the 13 periodic-
% orbit families, starting from the overlap / pair-result caches produced
% by run_rs4_all_pairs_summary.m.
%
% Workflow
%   1. Load all 13 family structs from cache; rebuild dense Xpo for each.
%   2. Enumerate valid pairs (filtered by minDVproxyMat from the batch
%      summary so pairs with no overlap are skipped without I/O).
%   3. Serial pre-loop: load O, compute V for each pair.
%      (Needs SA.Step4.rows_FRS_* — the large arrays — only here.)
%   4. Strip rows_FRS_* from family structs to make slim broadcast copies.
%   5. parfor (or serial-with-checkpoint): for each pair
%        load B  →  rs4_voxel_traj_extract  →  rs4_diffcorr
%   6. Tabulate exit flags and convergence.
%   7. Save results.mat + (in serial) checkpoint file.
%
% Outputs  (saved to <dcRoot>/rs4_dc_sweep_results.mat)
%   before_mat   [nPairs × 13]  numeric, pre-DC metrics
%   after_mat    [nPairs × 14]  numeric, post-DC metrics + solver metadata
%   traj_cell    {nPairs × 1}   struct per pair, sufficient for plot replay
%   pairs_ij     [nPairs × 2]   family index pairs (i < j)
%   families     {13 × 1}       family name strings
%   BEFORE_COLS / AFTER_COLS    column-name cell arrays

clear; clc;

% ---- path setup (mirrors run_rs4_all_pairs_summary.m) -------------------
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% =========================================================================
% USER SETTINGS
% Set PAIR_OUTPUT_ROOT to the folder that contains the per-pair sub-folders
% (e.g. .../rs3_results/20250110_143022/rs4_pairs_13fam).
% Leave empty to auto-detect the newest matching folder.
% =========================================================================
PAIR_OUTPUT_ROOT = '';   % <-- fill in if auto-detection misses it

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

% ---- config (must match what run_rs4_all_pairs_summary.m used) ----------
cfg = rs3_cfg_defaults();
cfg.cache.rebuild         = false;
cfg.io.save_figs          = false;
cfg.io.save_fig           = false;
cfg.io.fig_visible        = 'off';
cfg.diffcorr.display      = 'off';   % suppress fmincon output in batch
cfg.grid.dx               = 0.01;
cfg.grid.dy               = 0.01;
cfg.grid.dtheta           = deg2rad(2);
cfg.seed.ds_seed          = 0.02;
cfg.propag.Tmax           = pi;
cfg.fan.DV_cap_nd         = 0.2;
cfg.fan.dtheta_fan        = deg2rad(1.0);
cfg.propag.absTol         = 1e-9;
cfg.propag.relTol         = 1e-9;

N        = numel(families);
batchTag = sprintf('rs4_pairs_%dfam', N);   % 'rs4_pairs_13fam'

% ---- resolve pair output root -------------------------------------------
if isempty(PAIR_OUTPUT_ROOT)
    resultsRoot = fullfile(repoRoot, 'rs3_results');
    d = dir(resultsRoot);
    d = d([d.isdir] & ~strncmp({d.name}, '.', 1));
    [~, ord] = sort([d.datenum], 'descend');
    d = d(ord);
    for kd = 1:numel(d)
        candidate = fullfile(resultsRoot, d(kd).name, batchTag);
        if exist(candidate, 'dir')
            PAIR_OUTPUT_ROOT = candidate;
            break;
        end
    end
    if isempty(PAIR_OUTPUT_ROOT)
        error('[dc_sweep] Cannot auto-detect pair output root. Set PAIR_OUTPUT_ROOT manually.');
    end
end
fprintf('[dc_sweep] Pair output root : %s\n', PAIR_OUTPUT_ROOT);

% DC results go into a sibling folder so they don't pollute the pair cache
dcRoot = fullfile(fileparts(PAIR_OUTPUT_ROOT), 'rs4_dc_sweep');
if ~exist(dcRoot, 'dir'), mkdir(dcRoot); end
fprintf('[dc_sweep] DC results dir   : %s\n\n', dcRoot);

% =========================================================================
% 1. LOAD ALL 13 FAMILIES AND REBUILD Xpo
% =========================================================================
grid3 = rs3_grid_make(cfg);
fprintf('[dc_sweep] Loading %d family structs ...\n', N);

family_list = cell(N, 1);
for k = 1:N
    [family_list{k}, ~] = rs3_prepare_or_load_family(families{k}, cfg, grid3);
    family_list{k} = local_ensure_xpo(family_list{k}, ...
        cfg.propag.relTol, cfg.propag.absTol, 1001);
    fprintf('  [%2d/%2d]  %-30s  Nseeds_upper=%d\n', k, N, ...
        families{k}, size(family_list{k}.SeedsUpper, 1));
end

% =========================================================================
% 2. ENUMERATE VALID PAIRS
% =========================================================================
% Load batch summary to skip no-overlap pairs without touching disk.
% IMPORTANT: also load 'families' from the summary so the matrix row/column
% indices are correctly remapped to our local families list, regardless of
% the order they were saved in.
has_batch_summary = false;
dvFilter = zeros(N, N);   % default: try all pairs (0 = finite → passes filter)
batchSummaryFile = fullfile(PAIR_OUTPUT_ROOT, 'Summary', 'batch_summary_workspace.mat');
if exist(batchSummaryFile, 'file')
    bw = load(batchSummaryFile, 'minDVproxyMat', 'families');
    bw_mat = bw.minDVproxyMat;

    if isfield(bw, 'families') && ~isequal(bw.families(:), families(:))
        % The saved family order may differ → remap to local order.
        [~, loc] = ismember(families, bw.families);
        dvFilter = NaN(N, N);
        for ki = 1:N
            if loc(ki) == 0, continue; end
            for kj = 1:N
                if loc(kj) == 0, continue; end
                dvFilter(ki, kj) = bw_mat(loc(ki), loc(kj));
            end
        end
        fprintf('[dc_sweep] Loaded minDVproxyMat (remapped family order).\n');
    else
        dvFilter = bw_mat;
        fprintf('[dc_sweep] Loaded minDVproxyMat from batch summary.\n');
    end

    n_finite = sum(isfinite(dvFilter(:)));
    n_upper  = sum(sum(isfinite(triu(dvFilter, 1))));
    fprintf('[dc_sweep] dvFilter: %d finite entries total, %d in upper triangle.\n', ...
        n_finite, n_upper);
    has_batch_summary = true;
else
    fprintf('[dc_sweep] batch_summary_workspace.mat not found — will attempt all pairs.\n');
end

pairs_ij   = zeros(0, 2);
n_no_ovlap = 0;
n_no_file  = 0;
n_no_win   = 0;
for ii = 1:N
    for jj = ii+1:N
        % Skip pairs confirmed to have no overlap (check both triangles for safety).
        if has_batch_summary && ~isfinite(dvFilter(ii,jj)) && ~isfinite(dvFilter(jj,ii))
            n_no_ovlap = n_no_ovlap + 1;
            continue;
        end
        pairTag  = sprintf('%s__TO__%s', families{ii}, families{jj});
        pairSafe = rs3_sanitize_fname(pairTag);
        pairDir  = fullfile(PAIR_OUTPUT_ROOT, pairSafe);
        prFile   = fullfile(pairDir, ['rs4_' pairSafe '_pair_result.mat']);
        if ~exist(prFile, 'file')
            n_no_file = n_no_file + 1;
            continue;
        end
        pr = load(prFile);
        if ~isfield(pr, 'winnerMeta') || ~isfield(pr, 'B')
            n_no_win = n_no_win + 1;
            continue;
        end
        pairs_ij(end+1, :) = [ii, jj]; %#ok<AGROW>
    end
end
nPairs = size(pairs_ij, 1);
fprintf('\n[dc_sweep] Pair enumeration: %d valid | %d no-overlap | %d missing file | %d no winner\n', ...
    nPairs, n_no_ovlap, n_no_file, n_no_win);
fprintf('[dc_sweep] %d valid pairs to process  (%d total possible).\n\n', ...
    nPairs, N*(N-1)/2);

if nPairs == 0
    error('[dc_sweep] No valid pairs found. Check PAIR_OUTPUT_ROOT.');
end

% =========================================================================
% 3. SERIAL PRE-LOOP: compute V and load B for every pair
%    This step needs SA.Step4.rows_FRS_*, which are large.
%    We do it once here, serially, before parfor.
% =========================================================================
fprintf('[dc_sweep] Pre-computing voxel metadata V for all pairs ...\n');
V_cell = cell(nPairs, 1);
B_cell = cell(nPairs, 1);

for p = 1:nPairs
    ii = pairs_ij(p, 1);  jj = pairs_ij(p, 2);
    pairTag  = sprintf('%s__TO__%s', families{ii}, families{jj});
    pairSafe = rs3_sanitize_fname(pairTag);
    pairDir  = fullfile(PAIR_OUTPUT_ROOT, pairSafe);
    ovFile   = fullfile(pairDir, ['rs4_' pairSafe '_overlap.mat']);
    prFile   = fullfile(pairDir, ['rs4_' pairSafe '_pair_result.mat']);
    try
        ov = load(ovFile, 'O');
        pr = load(prFile, 'B');
        V_cell{p} = rs4_overlap_extract_voxel_info(family_list{ii}, family_list{jj}, ov.O, cfg);
        B_cell{p} = pr.B;
        fprintf('  [%2d/%2d]  %-26s → %s   (%d voxels)\n', p, nPairs, ...
            families{ii}, families{jj}, numel(V_cell{p}.ids));
    catch ME
        fprintf('  [%2d/%2d]  WARNING: V/B load failed for %s → %s: %s\n', ...
            p, nPairs, families{ii}, families{jj}, ME.message);
    end
end

% =========================================================================
% 4. SLIM FAMILY STRUCTS FOR PARFOR BROADCAST
%    Strip rows_FRS_upper/lower from Step4 — not needed after V is computed,
%    and they can be hundreds of MB per family.
% =========================================================================
slim_list = cell(N, 1);
for k = 1:N
    S = family_list{k};
    if isfield(S, 'Step4')
        if isfield(S.Step4, 'rows_FRS_upper'), S.Step4 = rmfield(S.Step4, 'rows_FRS_upper'); end
        if isfield(S.Step4, 'rows_FRS_lower'), S.Step4 = rmfield(S.Step4, 'rows_FRS_lower'); end
    end
    slim_list{k} = S;
end
clear family_list;   % free the heavy structs

% =========================================================================
% 5. PRE-ALLOCATE OUTPUT ARRAYS
% =========================================================================
BEFORE_COLS = { ...
    'pair_idx',  'fam_i',        'fam_j',       'vid', ...
    'DV_turn_A_mps', 'DV_patch_mps', 'DV_turn_B_mps', 'DV_total_mps', ...
    'dA_nd',     'dB_nd',        'delta_th_deg', ...
    'tof_A_days','tof_B_days'};

AFTER_COLS = { ...
    'pair_idx',  'fam_i',        'fam_j',       'vid', ...
    'DV_turn_A_mps', 'DV_patch_mps', 'DV_turn_B_mps', 'DV_total_mps', ...
    'delta_th_deg',  'tof_A_days',   'tof_B_days', ...
    'exitflag',  'converged',    'r_norm_scaled'};

nBC = numel(BEFORE_COLS);
nAC = numel(AFTER_COLS);

before_mat = NaN(nPairs, nBC);
after_mat  = NaN(nPairs, nAC);
traj_cell  = cell(nPairs, 1);

% =========================================================================
% 6. PARALLEL / SERIAL SETUP
% =========================================================================
use_par = license('test', 'Distrib_Computing_Toolbox') && nPairs > 1;
if use_par
    try
        if isempty(gcp('nocreate'))
            parpool('local');
        end
    catch ME
        warning('[dc_sweep] parpool failed (%s); falling back to serial.', ME.message);
        use_par = false;
    end
end
fprintf('\n[dc_sweep] Running in %s mode.\n\n', local_mode_str(use_par));

% =========================================================================
% 7. MAIN LOOP
% =========================================================================
chkFile = fullfile(dcRoot, 'rs4_dc_sweep_checkpoint.mat');

if use_par
    % ------------------------------------------------------------------
    % PARALLEL  — slim_list and cfg broadcast; V/B/pairs sliced per iter
    % ------------------------------------------------------------------
    parfor p = 1:nPairs
        [br, ar, tr] = local_process_pair( ...
            pairs_ij(p,:), p, families, slim_list, V_cell{p}, B_cell{p}, cfg, nBC, nAC);
        before_mat(p,:) = br;
        after_mat(p,:)  = ar;
        traj_cell{p}    = tr;
    end

else
    % ------------------------------------------------------------------
    % SERIAL  — print progress, checkpoint after every pair
    % ------------------------------------------------------------------
    for p = 1:nPairs
        ii = pairs_ij(p,1);  jj = pairs_ij(p,2);
        fprintf('[dc_sweep] %3d/%d  %-26s → %s\n', p, nPairs, ...
            families{ii}, families{jj});

        [br, ar, tr] = local_process_pair( ...
            pairs_ij(p,:), p, families, slim_list, V_cell{p}, B_cell{p}, cfg, nBC, nAC);
        before_mat(p,:) = br;
        after_mat(p,:)  = ar;
        traj_cell{p}    = tr;

        save(chkFile, 'before_mat', 'after_mat', 'traj_cell', 'pairs_ij', 'p', '-v7.3');
    end
end

% =========================================================================
% 8. EXIT-FLAG TABULATION
% =========================================================================
ef_col   = find(strcmp(AFTER_COLS, 'exitflag'),    1);
conv_col = find(strcmp(AFTER_COLS, 'converged'),   1);

ef_all   = after_mat(:, ef_col);
conv_all = after_mat(:, conv_col);
solved   = ~isnan(ef_all);

ef_vals   = ef_all(solved);
conv_vals = conv_all(solved);

fprintf('\n[dc_sweep] =========== EXIT FLAG SUMMARY ===========\n');
fprintf('[dc_sweep] %d of %d pairs solved.\n\n', sum(solved), nPairs);
fprintf('  %-12s  %6s   %16s   %7s\n', 'exitflag', 'count', 'converged / total', '%total');
fprintf('  %s\n', repmat('-', 1, 52));

for ef = sort(unique(ef_vals))'
    mask   = ef_vals == ef;
    n_conv = sum(conv_vals(mask) == 1);
    fprintf('  %+3d           %4d    %4d / %4d            %5.1f%%\n', ...
        ef, sum(mask), n_conv, sum(mask), 100*sum(mask)/sum(solved));
end

fprintf('  %s\n', repmat('-', 1, 52));
fprintf('  TOTAL          %4d    %4d / %4d            %5.1f%%\n', ...
    sum(solved), sum(conv_vals==1), sum(solved), 100*sum(conv_vals==1)/sum(solved));
fprintf('  Failed (NaN)   %4d\n', sum(~solved));

% =========================================================================
% 9. SAVE RESULTS
% =========================================================================
resFile = fullfile(dcRoot, 'rs4_dc_sweep_results.mat');
save(resFile, ...
    'before_mat', 'after_mat', 'traj_cell', ...
    'pairs_ij', 'families', ...
    'BEFORE_COLS', 'AFTER_COLS', ...
    'cfg', '-v7.3');
fprintf('\n[dc_sweep] Results saved → %s\n', resFile);

if ~use_par && exist(chkFile, 'file')
    delete(chkFile);
end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [before_row, after_row, traj] = local_process_pair( ...
        ij, p, families, slim_list, V, B, cfg, nBC, nAC)
% Run DC for one pair and pack results.  Called from both parfor and for.

    before_row = NaN(1, nBC);
    after_row  = NaN(1, nAC);
    traj       = [];

    if isempty(V) || isempty(B)
        return;
    end

    ii = ij(1);  jj = ij(2);
    SA = slim_list{ii};
    SB = slim_list{jj};

    try
        % --- warm-start trajectory ----------------------------------------
        T = rs4_voxel_traj_extract(SA, SB, V, B, cfg);

        % --- before-DC row ------------------------------------------------
        before_row = [ ...
            p, ii, jj, T.vid, ...
            T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps, T.DV_total_true_mps, ...
            T.dA_nd, T.dB_nd, rad2deg(T.delta_th_rad), ...
            T.tof_A_days, T.tof_B_days];

        % --- differential correction --------------------------------------
        Tc = rs4_diffcorr(T, SA, SB, cfg);

        % --- after-DC row -------------------------------------------------
        r_sc = norm(Tc.r_final ./ [SA.grid3.dx; SA.grid3.dy; SA.grid3.dtheta]);
        after_row = [ ...
            p, ii, jj, T.vid, ...
            Tc.DV_turn_A_mps, Tc.DV_patch_mps, Tc.DV_turn_B_mps, Tc.DV_total_mps, ...
            rad2deg(Tc.delta_th_rad), Tc.tof_A_days, Tc.tof_B_days, ...
            Tc.exitflag, double(Tc.converged), r_sc];

        % --- trajectory struct for later plotting -------------------------
        traj = struct();
        % identification
        traj.famA_idx        = ii;
        traj.famB_idx        = jj;
        traj.famA            = families{ii};
        traj.famB            = families{jj};
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
        % before-DC visualisation fields (for plot replay)
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
        traj.exitflag        = Tc.exitflag;
        traj.converged       = Tc.converged;
        traj.r_norm_scaled   = r_sc;
        % after-DC visualisation fields (for plot replay)
        traj.delta_th_rad    = Tc.delta_th_rad;
        traj.r_final         = Tc.r_final;

    catch ME
        warning('[dc_sweep] pair (%d,%d) %s → %s failed: %s', ...
            ii, jj, families{ii}, families{jj}, ME.message);
    end
end

% -------------------------------------------------------------------------

function S = local_ensure_xpo(S, relTol, absTol, N_po)
    if isfield(S,'Xpo') && ~isempty(S.Xpo) && ...
       isfield(S,'t_dense') && ~isempty(S.t_dense)
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

% -------------------------------------------------------------------------

function s = local_mode_str(use_par)
    if use_par, s = 'PARALLEL'; else, s = 'SERIAL (with checkpoint)'; end
end
