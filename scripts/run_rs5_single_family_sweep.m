%% RUN_RS5_SINGLE_FAMILY_SWEEP
% Transfer sweep from ONE origin family to a list of target families.
%
% For each (origin → target) pair this script:
%   1. Computes the forward/backward reachable-set overlap
%   2. Extracts the argmin-DV trajectory pair for the best overlap voxel
%   3. Runs differential correction to enforce patch continuity
%   4. Assembles a SINGLE concatenated trajectory  (departure → patch → arrival)
%   5. Saves figures with visible departure / arrival direction arrows
%   6. Saves a self-contained result.mat with everything needed to
%      re-integrate or re-plot the trajectories later without re-running
%      the full pipeline
%
% Output layout:
%   <out_root>/<tag>/<origin>/
%       sweep_summary.mat          —  all results + cfg in one file
%       <origin>_TO_<target>/
%           result.mat             —  T, Tc, traj_raw, traj_dc, regen recipe
%           rs5_*_before_dc.png    —  trajectory before DC  (with arrows)
%           rs5_*_after_dc.png     —  trajectory after DC   (with arrows)
%
% Parallel execution:  set use_parallel = true if RAM > nWorkers * 2 * ~500 MB.
% Each parfor worker receives a full copy of SA and the matching SB struct.
% With the default 2 targets this is very manageable; for large sweeps over
% many families monitor memory before enabling.
%
% Available families (from rs3_core_family_ic):
%   'Lyapunov L1'           'Lyapunov L2'
%   'Cycler 21'             'Cycler 11a'           'Cycler 11b'   'Cycler 32'
%   'Resonant 2to1 Stable'  'Resonant 2to1 Unstable'
%   'Resonant 3to1 Stable'  'Resonant 3to1 Unstable'
%   'Resonant 5to2 Stable'  'Resonant 5to2 Unstable'
%   'Distant Prograde Orbit'

clear; clc;

% ===================== USER KNOBS =====================
famOrigin  = 'Cycler 11a';
famTargets = {                          % all families except the origin
    'Lyapunov L1',
    'Lyapunov L2',
    'Cycler 21',
    'Cycler 11b',
    'Cycler 32',
    'Resonant 2to1 Stable',
    'Resonant 2to1 Unstable',
    'Resonant 3to1 Stable',
    'Resonant 3to1 Unstable',
    'Resonant 5to2 Stable',
    'Resonant 5to2 Unstable',
    'Distant Prograde Orbit',
};

use_parallel = false;  % true requires nWorkers * 2 * ~500 MB free RAM

cfg = rs3_cfg_defaults();

cfg.io.save_figs   = true;
cfg.io.fig_visible = 'on';

cfg.cache.enable  = true;
cfg.cache.rebuild = false;

% Grid — must match whichever cache was built
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

% Suppress intermediate bound heatmaps (we only need B.imin from that step)
cfg.plot.rs4.bounds_lb    = false;
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = false;

rs3_cfg_validate(cfg);

% ===================== OUTPUT ROOT =====================
origin_tag = rs3_sanitize_fname(famOrigin);
outroot    = fullfile(cfg.io.out_root, cfg.io.tag, origin_tag);
if ~exist(outroot, 'dir'), mkdir(outroot); end

fprintf('\n[sweep] Origin: %s\n', famOrigin);
fprintf('[sweep] Targets: %s\n', strjoin(famTargets, ', '));
fprintf('[sweep] Output:  %s\n\n', outroot);

% ===================== GRID + FAMILIES =====================
grid3 = rs3_grid_make(cfg);

fprintf('[sweep] Loading origin family ...\n');
[SA, ~] = rs3_prepare_or_load_family(famOrigin, cfg, grid3);

nT = numel(famTargets);
fprintf('[sweep] Loading %d target famil%s ...\n', nT, 'y/ies');
SB_cell = cell(nT, 1);
for iT = 1:nT
    fprintf('[sweep]   [%d/%d] %s\n', iT, nT, famTargets{iT});
    [SB_cell{iT}, ~] = rs3_prepare_or_load_family(famTargets{iT}, cfg, grid3);
end

% ===================== SWEEP =====================
results = cell(nT, 1);

if use_parallel && nT > 1
    % Each worker gets a broadcast copy of SA + one SB.
    % mkdir inside parfor needs try/catch to handle concurrent creation.
    parfor iT = 1:nT
        results{iT} = local_run_pair(SA, SB_cell{iT}, cfg, outroot);   %#ok<PFBNS>
    end
else
    for iT = 1:nT
        fprintf('\n[sweep] ===== Pair %d/%d : %s  ->  %s =====\n', ...
            iT, nT, famOrigin, famTargets{iT});
        results{iT} = local_run_pair(SA, SB_cell{iT}, cfg, outroot);
    end
end

% ===================== CONSOLE SUMMARY =====================
fprintf('\n[sweep] ======== SWEEP SUMMARY  ( origin: %s ) ========\n', famOrigin);
fprintf('%-28s  %8s  %8s  %8s  %9s  %9s  %7s  %s\n', ...
    'Target', 'DV_A', 'DV_ptch', 'DV_B', 'DV_raw', 'DV_dc', 'TOF(d)', 'DC?');
fprintf('%s\n', repmat('-', 1, 95));
for iT = 1:nT
    R = results{iT};
    if isempty(R) || ~isfield(R, 'Tc')
        err_msg = '';
        if isfield(R, 'error'), err_msg = R.error; end
        fprintf('%-28s  FAILED  %s\n', famTargets{iT}, err_msg);
    else
        fprintf('%-28s  %8.1f  %8.1f  %8.1f  %9.1f  %9.1f  %7.1f  %s\n', ...
            famTargets{iT}, ...
            R.T.DV_turn_A_mps,  R.T.DV_patch_mps, R.T.DV_turn_B_mps, ...
            R.T.DV_total_true_mps, R.Tc.DV_total_mps, ...
            R.traj_dc.tof_total_days, ...
            local_yesno(R.Tc.converged));
    end
end
fprintf('%s\n\n', repmat('=', 1, 95));

% ===================== SAVE SUMMARY =====================
save(fullfile(outroot, 'sweep_summary.mat'), ...
    'results', 'famOrigin', 'famTargets', 'cfg', '-v7.3');
fprintf('[sweep] Saved: %s/sweep_summary.mat\n', outroot);
fprintf('[sweep] Done.\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function R = local_run_pair(SA, SB, cfg, outroot)
%LOCAL_RUN_PAIR  Full rs4 pipeline for one (SA, SB) pair + save outputs.

famA     = SA.name;
famB     = SB.name;
pair_tag = [rs3_sanitize_fname(famA) '_TO_' rs3_sanitize_fname(famB)];
pair_dir = fullfile(outroot, pair_tag);

% mkdir is not parfor-safe without error handling
try
    if ~exist(pair_dir, 'dir'), mkdir(pair_dir); end
catch
end

R        = struct();
R.famA   = famA;
R.famB   = famB;
R.tag    = pair_tag;

try
    % ------------------------------------------------------------------
    % Step 1-3: overlap + voxel metadata + best voxel index
    % ------------------------------------------------------------------
    O = rs4_overlap_pair(SA, SB, cfg);
    V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg);
    B = rs4_overlap_visualize_bounds(V, SA, SB, O, cfg, pair_dir, pair_tag);

    fprintf('[%s]  Best voxel DVproxy = %.3f m/s  at (%.4f, %.4f)\n', ...
        pair_tag, B.min_dvproxy, B.x_at_min, B.y_at_min);

    % ------------------------------------------------------------------
    % Step 4: extract argmin-DV trajectory (re-integrate + true DV)
    % ------------------------------------------------------------------
    T = rs4_voxel_traj_extract(SA, SB, V, B, cfg);

    fprintf('[%s]  Raw DV: %.3f + %.3f + %.3f = %.3f m/s\n', pair_tag, ...
        T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps, T.DV_total_true_mps);

    % ------------------------------------------------------------------
    % Step 5: differential correction
    % ------------------------------------------------------------------
    fprintf('[%s]  Running differential correction ...\n', pair_tag);
    Tc = rs4_diffcorr(T, SA, SB, cfg);

    fprintf('[%s]  DC DV:  %.3f + %.3f + %.3f = %.3f m/s  (exitflag=%d, converged=%d)\n', ...
        pair_tag, Tc.DV_turn_A_mps, Tc.DV_patch_mps, Tc.DV_turn_B_mps, ...
        Tc.DV_total_mps, Tc.exitflag, Tc.converged);

    % ------------------------------------------------------------------
    % Step 6: build full (concatenated) trajectory structs
    % ------------------------------------------------------------------
    traj_raw = rs5_build_full_traj(T);   % before DC
    traj_dc  = rs5_build_full_traj(Tc);  % after  DC

    % ------------------------------------------------------------------
    % Step 7: visualise with departure / arrival arrows
    % ------------------------------------------------------------------
    rs5_visualize_full_traj(T,  traj_raw, SA, SB, cfg, pair_dir, [pair_tag '_before_dc']);
    rs5_visualize_full_traj(Tc, traj_dc,  SA, SB, cfg, pair_dir, [pair_tag '_after_dc']);

    % ------------------------------------------------------------------
    % Step 8: pack result and save
    % ------------------------------------------------------------------
    R.T        = T;
    R.Tc       = Tc;
    R.traj_raw = traj_raw;   % full traj before DC
    R.traj_dc  = traj_dc;    % full traj after  DC

    % Regen recipe: minimal data to re-integrate either trajectory from
    % scratch given only the family structs and cfg.
    R.regen.famA          = famA;
    R.regen.famB          = famB;
    % Uncorrected arc parameters
    R.regen.IC_A          = T.IC_A;       % [x;y;th] departure IC on PO_A
    R.regen.t_A           = T.t_A;        % forward integration time (nd)
    R.regen.IC_B_frs      = T.IC_B_frs;  % [x;y;th] FRS IC on PO_B
    R.regen.t_B           = T.t_B;        % forward integration time (nd)
    % Corrected arc parameters
    R.regen.alpha_A_dc    = Tc.alpha_A;   % PO phase at departure (rad)
    R.regen.delta_A_dc    = Tc.delta_A;   % heading kick at departure (rad)
    R.regen.t_A_dc        = Tc.t_A;
    R.regen.alpha_B_dc    = Tc.alpha_B;
    R.regen.delta_B_dc    = Tc.delta_B;
    R.regen.t_B_dc        = Tc.t_B;
    % Grid / propagation settings used (to check cache compatibility later)
    R.regen.cfg_grid      = cfg.grid;
    R.regen.cfg_propag    = cfg.propag;
    R.regen.cfg_fan       = cfg.fan;
    R.regen.cfg_seed      = cfg.seed;

    save(fullfile(pair_dir, 'result.mat'), 'R', 'cfg', '-v7.3');
    fprintf('[%s]  Saved: result.mat\n', pair_tag);

catch ME
    warning('[%s] FAILED: %s\n  %s', pair_tag, ME.message, ...
        ME.getReport('basic'));
    R.error = ME.message;
end
end

% -------------------------------------------------------------------------

function s = local_yesno(b)
if b, s = 'YES'; else, s = 'no'; end
end
