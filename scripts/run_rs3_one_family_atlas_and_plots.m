%% RUN_RS3_ONE_FAMILY_ATLAS_AND_PLOTS  —  build / load one family atlas and plot
%
% Builds (or loads from cache) a single-family reachable-set atlas and produces
% diagnostic figures for the Forward Reachable Set (FRS), Backward Reachable Set
% (BRS), seed distribution, and voxel grid.  This is the entry point for Steps
% 2–5 of the pipeline.  Run it first to verify that a family builds and plots
% correctly before committing to the larger batch sweeps.
%
% Prerequisites:
%   - Call rs3_setup (or run from the repository root) to put src/ on the path.
%   - No upstream scripts are required; atlases are built and cached here.
%
% User knobs (edit the USER KNOBS block near the top of the file):
%   familyName            — periodic orbit family to process (see list below)
%   cfg.grid.dx/dy        — spatial voxel width (nd); finer = more memory/time
%   cfg.grid.dtheta       — heading-angle voxel width (rad)
%   cfg.seed.ds_seed      — arclength spacing between seed points on the PO (nd)
%   cfg.propag.Tmax       — max integration time per arc, forward and backward (nd)
%   cfg.fan.DV_cap_nd     — max heading-offset DV budget at each seed (nd)
%   cfg.fan.dtheta_fan    — angular resolution of the heading fan (rad)
%   cfg.propag.absTol/relTol — ODE tolerances; tighten for finer spatial grids
%   cfg.io.save_figs      — true to write PNG/FIG files; false for display-only
%   cfg.io.fig_visible    — 'on' to see figures during the run; 'off' for batch
%   cfg.diag.zoom.*       — optional zoom window applied to trajectory figures
%   cfg.cache.rebuild     — true to recompute even if a valid cache already exists
%   cfg.par.enable        — true to parallelise over (seed × heading) jobs
%
% Available families:
%   'Lyapunov L1'            'Lyapunov L2'
%   'Cycler 21'              'Cycler 11a'    'Cycler 11b'    'Cycler 32'
%   'Resonant 2to1 Stable'   'Resonant 2to1 Unstable'
%   'Resonant 3to1 Stable'   'Resonant 3to1 Unstable'
%   'Resonant 5to2 Stable'   'Resonant 5to2 Unstable'
%   'Distant Prograde Orbit'
%
% Outputs written to rs3_results/<timestamp>/:
%   step5_<family>_atlas.mat   — atlas struct S + cache metadata (cacheInfo)
%   grid_validate_*.png/fig    — voxel grid diagnostic figure
%   family_validate_*.png/fig  — seed distribution along the periodic orbit
%   hits_validate_frs_*.png    — Forward Reachable Set voxel occupancy
%   hits_validate_brs_*.png    — Backward Reachable Set voxel occupancy

clear; clc;

% --- Repo setup (safe to call multiple times) ---
if exist('rs3_setup','file') ~= 2
    try
        repoRoot = fileparts(fileparts(which(mfilename)));
        addpath(repoRoot);
    catch
    end
end
if exist('rs3_setup','file') == 2
    rs3_setup();
end

% ====================== USER KNOBS ======================
familyName = 'Resonant 2to1 Unstable';      

cfg = rs3_cfg_defaults();

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
rs3_cfg_validate(cfg);

% ====================== OUTDIR ======================
outdir = fullfile(cfg.io.out_root, cfg.io.tag);
if ~exist(outdir, 'dir'), mkdir(outdir); end

% Optional log (doesn't break if logging files are removed)
if exist('rs3_log_init','file') == 2
    rs3_log_init(outdir, 'level', 'info');
    cleanupObj = onCleanup(@() rs3_log_close()); %#ok<NASGU>
else
    fprintf('[rs3] rs3_log_init not found -> continuing without logging.\n');
end


% ====================== STEP 2: GRID ======================
grid3 = rs3_grid_make(cfg);
rs3_grid_validate(grid3, cfg);

if cfg.io.save_figs
    rs3_grid_visual_validate(grid3, cfg, outdir);
end

% ====================== STEP 5: BUILD/LOAD ATLAS ======================
fprintf('[rs3] Building/loading atlas for "%s"...\n', familyName);
[S, cacheInfo] = rs3_prepare_or_load_family(familyName, cfg, grid3);

safeName = rs3_sanitize_fname(S.name);
save(fullfile(outdir, sprintf('step5_%s_atlas.mat', safeName)), 'S', 'cacheInfo', '-v7.3');

% ====================== PLOTS ======================
if cfg.io.save_figs
    rs3_cache_visual_validate(S, cfg, outdir, cacheInfo);
    rs3_family_visual_validate(S, cfg, outdir);
    rs3_hits_visual_validate(S, cfg, outdir);   % <-- now produces FRS and BRS plots w/ lower
end

fprintf('[rs3] Done.\n  outdir: %s\n', outdir);
