%% RUN_RS3_ONE_FAMILY_ATLAS_AND_PLOTS
% Build/load ONE family atlas (Steps 2–5) then plot FRS and BRS sets
% (Upper generated, Lower added via symmetry at plot time).

clear; clc;

% --- Repo setup (safe to call multiple times) ---
if exist('setup','file') ~= 2
    try
        repoRoot = fileparts(fileparts(which(mfilename)));
        addpath(repoRoot);
    catch
    end
end
if exist('setup','file') == 2
    setup();
end

% ====================== USER KNOBS ======================
familyName = 'Resonant 2to1 Unstable';      

cfg = atlas_cfg_defaults();

% Make it a "single family atlas" run
cfg.families.list = {familyName};
cfg.families.test_only = true;

% Grid / seeds / fan / propagation
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

% Output + plots
cfg.io.save_figs   = true;
cfg.io.save_fig    = true;      % also save .fig for interactive debugging
cfg.io.fig_visible = 'on';      % 'on' to see figures; 'off' for batch

% Optional zoom for the plots
cfg.diag.zoom.enable = false;
cfg.diag.zoom.xlim = [0.70 1.25];
cfg.diag.zoom.ylim = [-0.30 0.30];

% Cache (recommended)
cfg.cache.enable  = true;
cfg.cache.rebuild = false;

% Parallelism (optional)
cfg.par.enable = true;
cfg.par.mode   = 'jobs';

% ====================== VALIDATE ======================
atlas_cfg_validate(cfg);

% ====================== OUTDIR ======================
outdir = fullfile(cfg.io.out_root, cfg.io.tag);
if ~exist(outdir, 'dir'), mkdir(outdir); end



% ====================== STEP 2: GRID ======================
grid3 = atlas_grid_make(cfg);
atlas_grid_validate(grid3, cfg);

if cfg.io.save_figs
    atlas_grid_visual_validate(grid3, cfg, outdir);
end

% ====================== STEP 5: BUILD/LOAD ATLAS ======================
fprintf('[atlas] Building/loading atlas for "%s"...\n', familyName);
[S, cacheInfo] = atlas_prepare_or_load(familyName, cfg, grid3);

safeName = sanitize_fname(S.name);
save(fullfile(outdir, sprintf('step5_%s_atlas.mat', safeName)), 'S', 'cacheInfo', '-v7.3');

% ====================== PLOTS ======================
if cfg.io.save_figs
    atlas_cache_visual_validate(S, cfg, outdir, cacheInfo);
    atlas_family_visual_validate(S, cfg, outdir);
    atlas_hits_visual_validate(S, cfg, outdir);   % <-- now produces FRS and BRS plots w/ lower
end

fprintf('[atlas] Done.\n  outdir: %s\n', outdir);
