%% RUN_OVERLAP_PAIR_VISUALS
clear; clc;
% opengl software   % force software renderer — avoids JOGL/GPU deadlock with large scatter plots

% ===================== USER KNOBS =====================
famA = 'Lyapunov L1';
famB = 'Resonant 2to1 Unstable';   % <-- change this
tag  = [famA '_TO_' famB];

cfg = atlas_cfg_defaults();

cfg.io.save_figs   = true;
cfg.io.save_fig    = true;
cfg.io.fig_visible = 'on';

cfg.cache.enable  = true;
cfg.cache.rebuild = false;

% grid settings (must match for both atlases)
cfg.grid.dx               = 0.01;
cfg.grid.dy               = 0.01;
cfg.grid.dtheta           = deg2rad(5);
cfg.seed.ds_seed          = 0.05;
cfg.propag.Tmax           = pi/2;
cfg.fan.DV_cap_nd         = 0.2/1.5;
cfg.fan.dtheta_fan        = deg2rad(1);
cfg.propag.absTol         = 1e-8;
cfg.propag.relTol         = 1e-8;
cfg.propag.v2tol          = 1e-8;
cfg.log.step_len_factor   = 0.75;
cfg.log.maxstep_factor    = 2;

% zoom optional
cfg.diag.zoom.enable = false;
cfg.diag.zoom.xlim = [0.70 1.25];
cfg.diag.zoom.ylim = [-0.45 0.45];


% Figure toggles (optional)
cfg.plot.overlap.overlap_xy = true;
cfg.plot.overlap.overlap_xyz = false;
cfg.plot.overlap.combo_xy = true;
cfg.plot.overlap.combo_xyz = false;
cfg.plot.overlap.bounds_lb = true;
cfg.plot.overlap.bounds_ub = false;
cfg.plot.overlap.bounds_proxy = true;

atlas_cfg_validate(cfg);

outdir = fullfile(cfg.io.out_root, cfg.io.tag);
if ~exist(outdir,'dir'), mkdir(outdir); end

% ===================== BUILD GRID =====================
grid3 = atlas_grid_make(cfg);
if exist('atlas_grid_validate','file')==2
    atlas_grid_validate(grid3, cfg);
end

% ===================== LOAD/BUILD ATLASES =====================
[SA, infoA] = atlas_prepare_or_load(famA, cfg, grid3);
[SB, infoB] = atlas_prepare_or_load(famB, cfg, grid3);

save(fullfile(outdir, ['overlap_' sanitize_fname(tag) '_atlases.mat']), 'SA','SB','infoA','infoB','-v7.3');

% ===================== OVERLAP + VISUALS =====================
O = overlap_pair(SA, SB, cfg);
save(fullfile(outdir, ['overlap_' sanitize_fname(tag) '_overlap.mat']), 'O', '-v7.3');

% Extract voxel-wise candidate metadata for downstream ranking (no ranking yet)
V = overlap_extract_voxel_info(SA, SB, O, cfg);
save(fullfile(outdir, ['overlap_' sanitize_fname(tag) '_overlap_voxel_info.mat']), 'V', '-v7.3');

overlap_visualize(O, SA, SB, cfg, outdir, tag);
overlap_visualize_combo(SA, SB, O, cfg, outdir, tag);
B = overlap_visualize_bounds(V, SA, SB, O, cfg, outdir, tag);

if isstruct(B) && isfield(B,'min_dvproxy') && isscalar(B.min_dvproxy) && ~isempty(B.min_dvproxy) && isfinite(double(B.min_dvproxy))
    fprintf('[overlap] min DVproxy: %.3f m/s at (x,y)=(%.4f, %.4f), voxel #%d\n', ...
        B.min_dvproxy, B.x_at_min, B.y_at_min, B.imin);
end

fprintf('[overlap] Done.\n  outdir: %s\n', outdir);
