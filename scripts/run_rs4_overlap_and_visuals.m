%% RUN_RS4_OVERLAP_AND_VISUALS
clear; clc;

% ===================== USER KNOBS =====================
famA = 'Lyapunov L1';
famB = 'Resonant 2to1 Unstable';   % <-- change this
tag  = [famA '_TO_' famB];

cfg = rs3_cfg_defaults();   % you can rename later to rs4_cfg_defaults

cfg.io.save_figs   = true;
cfg.io.save_fig    = false;
cfg.io.fig_visible = 'on';

cfg.cache.enable  = true;
cfg.cache.rebuild = false;

% grid settings (must match for both atlases)
cfg.grid.dx     = 0.01;
cfg.grid.dy     = 0.01;
cfg.grid.dtheta = deg2rad(4);

cfg.seed.ds_seed   = 0.02;

% propagation/fan
cfg.propag.Tmax   = pi/2;
cfg.fan.DV_cap_nd = 0.10;
cfg.fan.dtheta_fan = deg2rad(2);

% zoom optional
cfg.diag.zoom.enable = false;
cfg.diag.zoom.xlim = [0.70 1.25];
cfg.diag.zoom.ylim = [-0.45 0.45];


% Figure toggles (optional)
cfg.plot.rs4.overlap_xy = true;
cfg.plot.rs4.overlap_xyz = false;
cfg.plot.rs4.combo_xy = true;
cfg.plot.rs4.combo_xyz = false;
cfg.plot.rs4.bounds_lb = true;
cfg.plot.rs4.bounds_ub = true;
cfg.plot.rs4.bounds_proxy = true;

rs3_cfg_validate(cfg);

outdir = fullfile(cfg.io.out_root, cfg.io.tag);
if ~exist(outdir,'dir'), mkdir(outdir); end

% ===================== BUILD GRID =====================
grid3 = rs3_grid_make(cfg);
if exist('rs3_grid_validate','file')==2
    rs3_grid_validate(grid3, cfg);
end

% ===================== LOAD/BUILD ATLASES =====================
[SA, infoA] = rs3_prepare_or_load_family(famA, cfg, grid3);
[SB, infoB] = rs3_prepare_or_load_family(famB, cfg, grid3);

save(fullfile(outdir, ['rs4_' rs3_sanitize_fname(tag) '_atlases.mat']), 'SA','SB','infoA','infoB','-v7.3');

% ===================== OVERLAP + VISUALS =====================
O = rs4_overlap_pair(SA, SB, cfg);
save(fullfile(outdir, ['rs4_' rs3_sanitize_fname(tag) '_overlap.mat']), 'O', '-v7.3');

% Extract voxel-wise candidate metadata for downstream ranking (no ranking yet)
V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg);
save(fullfile(outdir, ['rs4_' rs3_sanitize_fname(tag) '_overlap_voxel_info.mat']), 'V', '-v7.3');

rs4_overlap_visualize(O, SA, SB, cfg, outdir, tag);
rs4_overlap_visualize_combo(SA, SB, O, cfg, outdir, tag);
B = rs4_overlap_visualize_bounds(V, SA, SB, cfg, outdir, tag);

if isstruct(B) && isfield(B,'min_dvproxy') && isscalar(B.min_dvproxy) && ~isempty(B.min_dvproxy) && isfinite(double(B.min_dvproxy))
    fprintf('[rs4] min DVproxy: %.3f m/s at (x,y)=(%.4f, %.4f), voxel #%d\n', ...
        B.min_dvproxy, B.x_at_min, B.y_at_min, B.imin);
end

fprintf('[rs4] Done.\n  outdir: %s\n', outdir);
