function cfg = rs3_cfg_defaults()
%RS3_CFG_DEFAULTS  Default configuration schema for the rs3_* pipeline.
%
% Returns
%   cfg : struct with all knobs in one place (runner may override fields).

cfg = struct();

% ---------------- System / problem ----------------
cfg.sys = struct();
cfg.sys.RE_nd = 6378/384400;    % nondimensional Earth radius
cfg.sys.RM_nd = 1737/384400;    % nondimensional Moon radius

% Units (for readable reporting only; dynamics remain nondimensional)
cfg.units = struct();
cfg.units.LU_m        = 384400e3;                           % Earth-Moon distance (m)
cfg.units.period_days = 27.321661;                          % lunar sidereal period (days)
cfg.units.period_s    = cfg.units.period_days * 86400;
cfg.units.TU_s        = cfg.units.period_s / (2*pi);        % nondimensional time unit (s)
cfg.units.TU_days     = cfg.units.period_days / (2*pi);     % nondimensional time unit (days)
cfg.units.VU_mps      = cfg.units.LU_m / cfg.units.TU_s;   % nondimensional velocity unit (m/s)

% ---------------- Family selection ----------------
cfg.families = struct();
cfg.families.list = { ...
    'Lyapunov L1', ...
    'Lyapunov L2'  ...
    };
cfg.families.test_only = true;   % if true: only first family in list is processed

% ---------------- Seed / PO sampling ----------------
cfg.seed = struct();
cfg.seed.Tf_scale  = 1.0;    % PO period scale factor (1.0 = one full period)
cfg.seed.N_dense   = 2001;   % points for dense PO sampling
cfg.seed.y_eps     = 0.0;    % y-threshold for upper-half seed selection
cfg.seed.ds_seed   = 0.02;   % arclength spacing between seeds (nd)
cfg.seed.minSegPts = 5;      % minimum points per upper-half segment

% ---------------- Grid / voxelization ----------------
cfg.grid = struct();
cfg.grid.Rdom               = 1.2;          % domain radius (nd)
cfg.grid.dx                 = 0.01;         % x voxel width (nd)
cfg.grid.dy                 = 0.01;         % y voxel width (nd)
cfg.grid.dtheta             = deg2rad(2);   % theta voxel width (rad)
cfg.grid.enforce_y0_edge    = true;         % force y=0 to lie on a grid edge
cfg.grid.enforce_xy_symmetry = true;        % force grid symmetric about x- and y-axes

% ---------------- Steering fan (departure/arrival) ----------------
cfg.fan = struct();
cfg.fan.dtheta_fan = deg2rad(1);   % heading fan resolution (rad)
cfg.fan.DV_cap_nd  = 0.2;         % max DV budget for heading offset at seed (nd)

% ---------------- Propagation ----------------
cfg.propag = struct();
cfg.propag.Tmax   = pi;      % max integration time (nd, forward and backward)
cfg.propag.absTol = 1e-9;    % ODE absolute tolerance
cfg.propag.relTol = 1e-9;    % ODE relative tolerance
cfg.propag.v2tol  = 1e-8;    % reduced-model fallback threshold: if v^2 < v2tol, reintegrate with full 4D

% ---------------- Logging / segment-walk ----------------
cfg.log = struct();
cfg.log.step_len_factor = 0.5;    % step_len = factor * min(dx, dy)
cfg.log.maxstep_factor  = 0.5;    % ODE MaxStep = maxstep_factor * step_len

cfg.log.segwalk = struct();
cfg.log.segwalk.enable = true;    % enable segment-walk subsampling
cfg.log.segwalk.frac   = 0.25;   % substep <= frac * cell width (smaller = safer, more work)

% ---------------- Caching ----------------
cfg.cache = struct();
cfg.cache.enable        = true;
cfg.cache.dir           = fullfile(rs3_repo_root(), 'rs3_cache');
cfg.cache.rebuild       = false;                  % if true: ignore existing cache and recompute
cfg.cache.version_tag   = 'rs3_v2_keep_masked';  % bumped when cache schema changes
cfg.cache.store_dense_po = false;                 % if true: store full dense PO in cache file

% ---------------- Overlap ----------------
cfg.overlap = struct();
cfg.overlap.primary_buffer_frac = 0.05;   % fractional exclusion buffer around Earth/Moon (applied to RE_nd/RM_nd)

cfg.rs4 = struct();
cfg.rs4.extract = struct();
cfg.rs4.extract.parallel = false;  % parallelize overlap voxel metadata extraction (requires parpool)
cfg.rs4.dc = struct();
cfg.rs4.dc.jacobian_mode = 'auto';  % {'stm','fd','auto'}; auto prefers STM and falls back to FD
cfg.rs4.dc.fd_step = 1e-6;                % central-difference perturbation for Jacobian fallback
cfg.rs4.dc.stm_max_sensitivity = 1e6;   % if STM |dthf/dth0| exceeds this, fallback to FD in auto mode

% ---------------- Parallelism ----------------
cfg.par = struct();
cfg.par.enable         = true;   % enable PARFOR over (seed x heading) jobs
cfg.par.progress_every = 50;     % print progress every N jobs in PARFOR

% ---------------- Output ----------------
cfg.io = struct();
cfg.io.out_root      = fullfile(rs3_repo_root(), 'rs3_results');
cfg.io.tag           = datestr(now, 'yyyymmdd_HHMMSS');  %#ok<DATST>
cfg.io.save_figs     = true;    % save PNG figures
cfg.io.save_fig      = true;    % also save .fig for interactive use
cfg.io.fig_subdir    = '';      % subfolder under outdir for figures ('' = none)
cfg.io.fig_resolution = 220;    % PNG DPI
cfg.io.fig_visible   = 'off';   % figure visibility during rendering
cfg.io.verbose       = true;    % print progress messages

% ---------------- Diagnostics ----------------
cfg.diag = struct();
cfg.diag.progress    = true;    % enable per-job progress reporting in PARFOR
cfg.diag.smoke_checks = true;   % run lightweight sanity checks during seed preparation
cfg.diag.show_orbits = true;    % overlay PO arcs in visual validate figures
cfg.diag.po_stride   = 8;       % stride for PO arc subsampling in plots

cfg.diag.zoom = struct();
cfg.diag.zoom.enable = true;           % apply zoom window to trajectory/overlap figures
cfg.diag.zoom.xlim   = [0.70 1.25];   % zoom x-limits (nd), centered near Moon
cfg.diag.zoom.ylim   = [-0.30 0.30];  % zoom y-limits (nd)

end
