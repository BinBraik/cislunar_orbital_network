function cfg = rs3_cfg_defaults()
%RS3_CFG_DEFAULTS  Default configuration schema for the rs3_* pipeline.
%
% Returns
%   cfg : struct with all knobs in one place (runner may override fields).

cfg = struct();

% ---------------- System / problem ----------------
cfg.sys = struct();
cfg.sys.model = 'planar_earth_moon_cr3bp';
cfg.sys.mu = [];                    % if empty: obtained from family_ic
cfg.sys.CJ_policy = 'from_family';  % 'from_family' (recommended), 'fixed' (optional later)
cfg.sys.CJ_fixed = NaN;

% Nondimensional radii (Earth/Moon) used for stop events (baseline values)
cfg.sys.RE_nd = 6378/384400;
cfg.sys.RM_nd = 1737/384400;

% Units (for readable reporting only; dynamics remain nondimensional)
cfg.units = struct();
cfg.units.LU_m = 384400e3;              % Earth-Moon distance (m)
cfg.units.period_days = 27.321661;      % lunar sidereal period (days)
cfg.units.period_s = cfg.units.period_days * 86400;
cfg.units.TU_s = cfg.units.period_s/(2*pi);
cfg.units.VU_mps = cfg.units.LU_m / cfg.units.TU_s;
cfg.units.TU_days = cfg.units.period_days / (2*pi);   % nondimensional time unit in days

% ---------------- Family selection ----------------
cfg.families = struct();
cfg.families.list = { ...
    'Lyapunov L1', ...
    'Lyapunov L2'  ...
    };
cfg.families.test_only = true;      % Step 1–4: do not run full matrix by default

% ---------------- Seed / PO sampling ----------------
cfg.seed = struct();
cfg.seed.Tf_scale = 1.0;
cfg.seed.N_dense = 2001;
cfg.seed.y_eps = 0.0;
cfg.seed.ds_seed = 0.01;
cfg.seed.minSegPts = 5;

% ---------------- Grid / voxelization ----------------
cfg.grid = struct();
cfg.grid.Rdom = 1.2;
cfg.grid.dx = 0.01;
cfg.grid.dy = 0.01;
cfg.grid.dtheta = deg2rad(4);       % voxel theta bin width (rad)
cfg.grid.theta_convention = '[-pi,pi)';
cfg.grid.enforce_y0_edge = true;
cfg.grid.enforce_xy_symmetry = true;

% ---------------- Steering fan (departure/arrival) ----------------
cfg.fan = struct();
cfg.fan.dtheta_fan = deg2rad(1);    % heading fan resolution (rad)
cfg.fan.DV_cap_nd = 0.1;         % nondimensional cap for steering at seed
cfg.fan.delta_policy = 'symmetric';

% Back-compat alias (older cfg field name)
cfg.fan.dtheta = cfg.fan.dtheta_fan;

% ---------------- Propagation ----------------
cfg.propag = struct();
cfg.propag.Tmax = pi/2;
cfg.propag.absTol = 1e-11;
cfg.propag.relTol = 1e-11;
cfg.propag.v2tol = 1e-10;           % stop if v^2 drops below this tolerance

% ---------------- Logging ----------------
cfg.log = struct();
cfg.log.step_len_factor = 0.20;     % step_len = factor * min(dx,dy)
cfg.log.maxstep_factor = 0.25;      % ode MaxStep = factor * step_len (baseline used 0.25*step_len)

cfg.log.segwalk = struct();
cfg.log.segwalk.enable = true;
cfg.log.segwalk.frac = 0.25;        % substep <= frac*cell (smaller => safer, more work)

% Back-compat alias
cfg.log.maxStep_factor = cfg.log.maxstep_factor;

% ---------------- Caching ----------------
cfg.cache = struct();
cfg.cache.enable = true;
cfg.cache.dir = fullfile(rs3_repo_root(), 'rs3_cache');
cfg.cache.rebuild = false;
cfg.cache.version_tag = 'rs3_v2_keep_masked';  % v2: packed rows + center-level Keep filtering
cfg.cache.store_entry_state = false;
cfg.cache.store_dense_po = false;

% ---------------- LEGACY REFINEMENT (DISABLED BY DEFAULT) ----------------
cfg.refine = struct();
cfg.refine.enable = true;   % master switch for Step 7 (any refinement)
cfg.refine.dx_min = 0.0005;
cfg.refine.dy_min = 0.0005;
cfg.refine.dtheta_min = deg2rad(0.25);
cfg.refine.maxLevels = 1;
cfg.refine.maxRegions = 1500;
cfg.refine.maxLocalJobs = 2e5;
cfg.refine.maxJobsPerRegion = cfg.refine.maxLocalJobs; % alias used by Step7
cfg.refine.stop_on_exit = true;
cfg.refine.t_local_max = 1.5;
cfg.refine.split_policy = 'both';
cfg.refine.pad_vox = 1;
cfg.refine.nsamp_per_voxel = 2;
cfg.refine.maxVoxelsPerRegion = 800;

% ---------------- NEW REFINEMENT (Step 7, enabled by default) ----------------
cfg.refine7 = struct();
cfg.refine7.enable = true;          % master switch – ON by default
cfg.refine7.M_pairs_frac = 0.10;    % fraction of unique pairs to keep after dedup
cfg.refine7.M_pairs_min  = 50;     % floor: never fewer than this
cfg.refine7.M_pairs_max  = 2000;   % cap: never more than this
cfg.refine7.N_resample = 200;       % points per curve for proximity search
cfg.refine7.d_min_gate = 0.005;     % proximity threshold (nd) – ~1900 km
cfg.refine7.N_top_global = 50;      % final top candidates to keep
cfg.refine7.optim_method = 'fminsearch'; % 'fminsearch' or 'gridsearch'

% ---------------- Candidate retention / scoring ----------------
cfg.cand = struct();
cfg.cand.K_per_voxel = 10;
cfg.cand.K_pairs_per_voxel = 50;  % top pairs kept per voxel (pairwise scoring)
cfg.cand.bound_mode = 'conservative_lb';
cfg.cand.maxPairsPerVoxel = 200;

% ---------------- Step 8: 3-impulse scoring ----------------
cfg.score = struct();
cfg.score.enable = true;
cfg.score.objective = 'dvtotal';     % dvtotal | tof | pareto
cfg.score.maxPairsPerVoxel = 200;     % cap per overlap voxel (after trimming)
cfg.score.maxCandPerVoxel = 40;       % cap candidates kept per side per voxel before pairing
cfg.score.reintegrate_enable = true;
cfg.score.reintegrate_topN = 25;
cfg.score.brs_theta_convention = 'as_is'; % as_is (default) | flip_pi
cfg.score.save_top_global = 500;      % keep only best global candidates in memory

% ---------------- Parallelism ----------------
cfg.par = struct();
cfg.par.enable = true;
cfg.par.mode = 'jobs';              % per your requirement: (seed x heading) jobs
cfg.par.pool_size = [];             % empty => MATLAB default/existing pool
cfg.par.batch_size = 1000;          % reserved for later (chunking)
cfg.par.progress_every = 50;        % progress print cadence for PARFOR (jobs)

% Back-compat alias
cfg.par.use = cfg.par.enable;

% ---------------- Output ----------------
cfg.io = struct();
cfg.io.out_root = fullfile(rs3_repo_root(), 'rs3_results');
cfg.io.tag = datestr(now, 'yyyymmdd_HHMMSS');
cfg.io.save_mat = true;
cfg.io.save_csv = true;
cfg.io.save_figs = true;
cfg.io.save_fig = true;          % also save MATLAB .fig for interactive debugging
cfg.io.fig_subdir = '';     % subfolder under outdir for figures
cfg.io.fig_resolution = 220;     % default PNG resolution
cfg.io.fig_visible = 'off';
cfg.io.verbose = true;

% ---------------- Validation / diagnostics ----------------
cfg.diag = struct();
cfg.diag.plot_each_level = false;
cfg.diag.maxPlotFootprint = 8000;
cfg.diag.maxPlotOverlap = 15000;
cfg.diag.plot_style = 'story';
cfg.diag.story_plots = true;
cfg.diag.show_orbits = true;
cfg.diag.po_stride = 8;
cfg.diag.zoom = struct();
cfg.diag.zoom.enable = true;
cfg.diag.zoom.xlim = [0.70 1.25];
cfg.diag.zoom.ylim = [-0.30 0.30];
cfg.diag.zoom.show_gridlines = true;
cfg.diag.zoom.max_gridlines = 60;

cfg.diag.estimate_memory = true;
cfg.diag.print_cfg = true;
cfg.diag.smoke_checks = true;
cfg.diag.enable_asserts = true;
cfg.diag.progress = true;

% Step toggles (these are for the old pipeline; kept for compatibility)
cfg.diag.run_step3 = true;
cfg.diag.run_step4 = true;

end