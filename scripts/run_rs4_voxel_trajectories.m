%% RUN_RS4_VOXEL_TRAJECTORIES  —  extract, correct, and plot true-DV trajectory pair
%
% Re-integrates the argmin-DV arc pair for the best overlap voxel of a single
% family pair and applies differential correction.  This gives the true
% DV_total (replacing the proxy upper-bound) and produces side-by-side before/
% after DC trajectory figures.  Use this script to inspect and verify a specific
% transfer in detail before running the full 78-pair batch.
%
% Workflow:
%   1. Load both family atlases from cache (or rebuild if needed).
%   2. Compute overlap O and per-voxel metadata V (same as run_rs4_overlap_and_visuals).
%   3. Find the best-DV voxel (B struct with B.imin) via rs4_overlap_visualize_bounds.
%   4. Re-integrate the two winning arcs → T struct with true DV values.
%      DV_total_true = DV_turn_A_min + DV_patch_true + DV_turn_B_min
%      where DV_patch_true = 2·v_box·sin(|Δθ|/2) at the voxel centre.
%   5. Run differential correction → Tc struct (tightened DV + corrected arcs).
%   6. Generate cislunar trajectory figure and before/after DC comparison figure.
%
% Prerequisites:
%   - Cached atlases for famA and famB must exist in rs3_cache/.
%   - cfg below (grid / fan / propag) must match those cached atlases exactly.
%
% User knobs:
%   famA, famB              — the two orbit families to analyse
%   cfg.grid/fan/propag     — must match cached atlases
%   cfg.diffcorr.*          — DC solver settings:
%     tol_patch             — normalised convergence tolerance for fmincon
%     tol_converged         — threshold for the "CONVERGED" status label
%     display               — 'iter' shows iteration-by-iteration output;
%                             'off' for silent; 'final' for summary only
%     MaxIterations         — fmincon iteration budget (raise to 600 for hard cases)
%     MaxFunEvals           — fmincon function-evaluation budget
%     N_po_dt / N_po_min    — PO spline knot density (finer = slower but smoother)
%
% Outputs written to rs3_results/<timestamp>/:
%   rs4_<tag>_voxel_traj.png/fig           — re-integrated arc pair figure
%   rs4_<tag>_diffcorr_compare.png/fig     — before/after DC side-by-side figure

clear; clc;

% ===================== USER KNOBS =====================
famA = 'Lyapunov L1';
famB = 'Resonant 2to1 Unstable';   % <-- change this pair
tag  = [famA '_TO_' famB];

cfg = rs3_cfg_defaults();

cfg.io.save_figs   = true;
cfg.io.fig_visible = 'on';

cfg.cache.enable  = true;
cfg.cache.rebuild = false;

% Grid (must match the cache that was built)
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
% Disable all bound plots (we only need B.imin from visualize_bounds)
cfg.plot.rs4.bounds_lb    = false;
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = false;

% ---- Differential Correction (rs4_diffcorr knobs) ----
cfg.diffcorr.tol_patch     = 1e-4;   % normalized convergence threshold
cfg.diffcorr.tol_converged = 1e-4;   % report CONVERGED if ||r_sc|| <= this (can be >= tol_patch)
cfg.diffcorr.display       = 'iter'; % 'off' | 'iter' | 'final'  (keep iter for single-pair debugging)
cfg.diffcorr.MaxIterations = 300;    % fmincon iteration budget (raise to 600 for hard cases)
cfg.diffcorr.MaxFunEvals   = 8000;   % fmincon function eval budget
cfg.diffcorr.N_po_dt       = 0.003;  % PO knot spacing [ND] — auto-scales to orbit period
cfg.diffcorr.N_po_min      = 1001;   % minimum knot count (floor)

rs3_cfg_validate(cfg);

outdir = fullfile(cfg.io.out_root, cfg.io.tag);
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ===================== GRID + ATLASES =====================
grid3 = rs3_grid_make(cfg);

[SA, ~] = rs3_prepare_or_load_family(famA, cfg, grid3);
[SB, ~] = rs3_prepare_or_load_family(famB, cfg, grid3);

% ===================== OVERLAP + VOXEL METADATA =====================
O = rs4_overlap_pair(SA, SB, cfg);
V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg);

% Get B struct (best voxel index) — suppress figure output via cfg flags
B = rs4_overlap_visualize_bounds(V, SA, SB, O, cfg, outdir, tag);

fprintf('\n[traj] Best voxel: DVproxy=%.3f m/s at (%.4f, %.4f), voxel #%d\n', ...
    B.min_dvproxy, B.x_at_min, B.y_at_min, B.imin);

% ===================== EXTRACT TRUE TRAJECTORY =====================
T = rs4_voxel_traj_extract(SA, SB, V, B, cfg);

% ===================== PRINT SUMMARY =====================
fprintf('\n========== DV Summary (%s -> %s) ==========\n', famA, famB);
fprintf('  DV_turn_A        = %8.3f m/s\n', T.DV_turn_A_mps);
fprintf('  DV_patch_true    = %8.3f m/s   (delta_th = %.4f deg)\n', ...
    T.DV_patch_mps, rad2deg(T.delta_th_rad));
fprintf('  DV_turn_B        = %8.3f m/s\n', T.DV_turn_B_mps);
fprintf('  ------------------------------------------------\n');
fprintf('  DV_total_true    = %8.3f m/s\n', T.DV_total_true_mps);
fprintf('  DV_proxy (bound) = %8.3f m/s   (tightening: %.3f m/s)\n', ...
    T.DV_proxy_mps, T.DV_proxy_mps - T.DV_total_true_mps);
fprintf('\n  TOF_A = %.2f days,  TOF_B = %.2f days\n', T.tof_A_days, T.tof_B_days);
fprintf('  miss dA = %.5f nd,  dB = %.5f nd\n', T.dA_nd, T.dB_nd);
fprintf('=================================================\n\n');

% ===================== PLOT (original arcs) =====================
rs4_voxel_traj_visualize_single(T, SA, SB, cfg, outdir, tag);

% ===================== DIFFERENTIAL CORRECTION =====================
fprintf('\n[traj] Running differential correction ...\n');
Tc = rs4_diffcorr(T, SA, SB, cfg);

% ===================== PRINT DC SUMMARY =====================
fprintf('\n========== Differential Correction Summary (%s -> %s) ==========\n', famA, famB);
fprintf('  BEFORE:  DV_turn_A=%8.3f  +  DV_patch=%8.3f  +  DV_turn_B=%8.3f  =  %8.3f m/s\n', ...
    T.DV_turn_A_mps,  T.DV_patch_mps,  T.DV_turn_B_mps,  T.DV_total_true_mps);
fprintf('  AFTER:   DV_turn_A=%8.3f  +  DV_patch=%8.3f  +  DV_turn_B=%8.3f  =  %8.3f m/s\n', ...
    Tc.DV_turn_A_mps, Tc.DV_patch_mps, Tc.DV_turn_B_mps, Tc.DV_total_mps);
fprintf('  Reduction: %.3f m/s  |  fmincon exitflag=%d  iterations=%d\n', ...
    T.DV_total_true_mps - Tc.DV_total_mps, Tc.exitflag, Tc.iterations);
fprintf('  TOF_A: %.2f -> %.2f days  |  TOF_B: %.2f -> %.2f days\n', ...
    T.tof_A_days, Tc.tof_A_days, T.tof_B_days, Tc.tof_B_days);
fprintf('=================================================================\n\n');

% ===================== COMPARISON PLOT =====================
rs4_diffcorr_visualize(T, Tc, SA, SB, cfg, outdir, tag);

fprintf('[traj] Done. Output: %s\n', outdir);
