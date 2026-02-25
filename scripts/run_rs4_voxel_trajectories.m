%% RUN_RS4_VOXEL_TRAJECTORIES
% Re-integrate the argmin-DV trajectory pair for the best overlap voxel of
% one family pair and compute the true DV_total.
%
% Workflow:
%   1. Load (or rebuild) both family atlases from cache
%   2. Compute overlap O and voxel metadata V (same as run_rs4_overlap_and_visuals)
%   3. Run rs4_overlap_visualize_bounds to get B struct (best voxel index B.imin)
%   4. Call rs4_voxel_traj_extract  -> T struct with arcs + true DV values
%   5. Call rs4_voxel_traj_visualize_single -> cislunar background plot
%
% The DVtotal printed here replaces the proxy upper-bound with:
%   DV_total_true = DV_turn_A_min + DV_patch_true + DV_turn_B_min
% where DV_patch_true = 2*v_box_center*sin(|theta_A - theta_B|/2) measured
% at the closest points of each re-integrated arc to the voxel center.

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
cfg.grid.dx     = 0.01;
cfg.grid.dy     = 0.01;
cfg.grid.dtheta = deg2rad(4);

cfg.seed.ds_seed   = 0.02;

% Propagation / fan
cfg.propag.Tmax    = pi/2;
cfg.fan.DV_cap_nd  = 0.1;
cfg.fan.dtheta_fan = deg2rad(2.0);
cfg.propag.absTol  = 1e-9;
cfg.propag.relTol  = 1e-9;
cfg.propag.v2tol   = 1e-8;
cfg.log.step_len_factor = 0.5;
cfg.log.maxstep_factor  = 0.5;

% Differential-correction (DC) options for voxel patch refinement.
% Keep disabled by default so current DV_proxy/DV_patch_true outputs and
% figure behavior remain unchanged unless explicitly enabled.
cfg.rs4.dc.enable             = false;
cfg.rs4.dc.max_iter           = 25;
cfg.rs4.dc.tol_residual       = 1e-10;
cfg.rs4.dc.tol_step           = 1e-10;
cfg.rs4.dc.tol_dv_improve     = 1e-6;
cfg.rs4.dc.lambda_init        = 1e-2;   % Levenberg-Marquardt damping
cfg.rs4.dc.line_search_enable = true;
cfg.rs4.dc.phase_enable       = true;
cfg.rs4.dc.weights.position   = 1.0;
cfg.rs4.dc.weights.heading    = 1.0;

% Disable all bound plots (we only need B.imin from visualize_bounds)
cfg.plot.rs4.bounds_lb    = false;
cfg.plot.rs4.bounds_ub    = false;
cfg.plot.rs4.bounds_proxy = false;

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

% ===================== PLOT =====================
rs4_voxel_traj_visualize_single(T, SA, SB, cfg, outdir, tag);

fprintf('[traj] Done. Output: %s\n', outdir);
