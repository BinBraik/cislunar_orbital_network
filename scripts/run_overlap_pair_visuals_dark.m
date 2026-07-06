%% RUN_OVERLAP_PAIR_VISUALS_DARK
%
% Dark-theme presentation outputs for an overlap pair, kept separate from
% RUN_OVERLAP_PAIR_VISUALS so the standard (light) figures are untouched.
%
% Produces, alongside the standard atlas/overlap cache files:
%   - overlap_<tag>_overlap_xy_dark.fig / .pdf / .eps   (paper/poster figure)
%   - overlap_<tag>_process.gif                          (background+POs -> FRS -> BRS -> overlap)

clear; clc;

% ===================== USER KNOBS =====================
famA = 'Lyapunov L1';
famB = 'Resonant 2to1 Unstable';   % <-- change this
tag  = [famA '_TO_' famB];

cfg = atlas_cfg_defaults();

cfg.io.save_figs   = true;
cfg.io.save_fig    = true;
cfg.io.fig_visible = 'off';   % rendering is headless; dark figure/GIF are saved to disk

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

% zoom optional (applied to the dark figure and the GIF alike)
cfg.diag.zoom.enable = false;
cfg.diag.zoom.xlim = [0.70 1.25];
cfg.diag.zoom.ylim = [-0.45 0.45];

% GIF pacing knobs
GIF_FPS             = 20;
GIF_HOLD_START_SEC  = 1.5;
GIF_HOLD_FRS_SEC    = 1.0;
GIF_HOLD_BRS_SEC    = 1.0;
GIF_HOLD_END_SEC    = 3.0;
GIF_N_BUILD_FRS     = 12;
GIF_N_BUILD_BRS     = 12;
GIF_N_BUILD_OVERLAP = 10;
GIF_FIG_SIZE        = [900 700];
GIF_MAX_PTS_PER_SET = 20000;  % voxel-cloud cap for the animated scatter; lower this
                              % (e.g. 4000) only if a given family pair renders too slowly

% Dark figure knobs
FIG_MAX_PTS_PER_SET = Inf;    % plot every voxel (PDF/EPS is rasterized, so this is free)
FIG_PDF_RESOLUTION  = 300;    % DPI for the rasterized PDF/EPS embed (keeps files small)

atlas_cfg_validate(cfg);

outdir = fullfile(cfg.io.out_root, cfg.io.tag);
if ~exist(outdir,'dir'), mkdir(outdir); end

% ===================== BUILD GRID =====================
grid3 = atlas_grid_make(cfg);
if exist('atlas_grid_validate','file')==2
    atlas_grid_validate(grid3, cfg);
end

% ===================== LOAD/BUILD ATLASES =====================
[SA, infoA] = atlas_prepare_or_load(famA, cfg, grid3); %#ok<ASGLU>
[SB, infoB] = atlas_prepare_or_load(famB, cfg, grid3); %#ok<ASGLU>

% ===================== OVERLAP =====================
O = overlap_pair(SA, SB, cfg);

% ===================== DARK PRESENTATION FIGURE =====================
overlap_visualize_dark(O, SA, SB, cfg, outdir, tag, ...
    'MaxPointsPerSet', FIG_MAX_PTS_PER_SET, ...
    'Resolution', FIG_PDF_RESOLUTION);

% ===================== PROCESS GIF =====================
gifPath = overlap_visualize_gif(O, SA, SB, cfg, outdir, tag, ...
    'FPS', GIF_FPS, ...
    'HoldStart', GIF_HOLD_START_SEC, ...
    'HoldFRS', GIF_HOLD_FRS_SEC, ...
    'HoldBRS', GIF_HOLD_BRS_SEC, ...
    'HoldEnd', GIF_HOLD_END_SEC, ...
    'NBuildFRS', GIF_N_BUILD_FRS, ...
    'NBuildBRS', GIF_N_BUILD_BRS, ...
    'NBuildOverlap', GIF_N_BUILD_OVERLAP, ...
    'FigSize', GIF_FIG_SIZE, ...
    'MaxPointsPerSet', GIF_MAX_PTS_PER_SET);

fprintf('[overlap-dark] Done.\n  outdir: %s\n  gif:    %s\n', outdir, gifPath);
