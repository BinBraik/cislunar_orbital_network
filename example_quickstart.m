%% EXAMPLE_QUICKSTART
% End-to-end demonstration of the cislunar orbital network pipeline.
%
% Runs in ~2-5 minutes on a single CPU (serial, coarse grid).
% For paper-quality results use the default grid settings and run
% with cfg.par.enable = true (Parallel Computing Toolbox required).
%
% Stages demonstrated:
%   1. Configure a fast coarse grid
%   2. Build reachable-set atlases for two families (LL1 and LL2)
%   3. Compute FRS(LL1) ∩ BRS(LL2) overlap
%   4. Extract DV proxy metrics for the best overlap voxel
%   5. Print a summary
%
% See README.md for the full execution order and runner scripts.

%% 0. Add paths
setup

%% 1. Configuration (coarse for speed)
cfg = atlas_cfg_defaults();

% Coarser grid — fast enough to run without a parpool
cfg.grid.dx     = 0.04;          % default 0.01 nd
cfg.grid.dy     = 0.04;
cfg.grid.dtheta = deg2rad(8);    % default 2 deg

% Shorter integration time
cfg.propag.Tmax = pi / 2;        % default pi nd (~6.8 days each direction)

% Narrower heading fan
cfg.fan.DV_cap_nd  = 0.1;        % default 0.2 nd  (~102 m/s)
cfg.fan.dtheta_fan = deg2rad(2); % default 1 deg

% Serial execution (remove or set true if you have a parpool)
cfg.par.enable = false;

% Output folder
cfg.io.out_root    = fullfile(repo_root(), 'quickstart_output');
cfg.io.fig_visible = 'on';
cfg.io.save_figs   = false;      % set true to save PNGs

% Cache in a separate folder so it doesn't pollute the main cache
cfg.cache.dir = fullfile(repo_root(), 'atlas_cache_quickstart');

atlas_cfg_validate(cfg);

%% 2. Build voxel grid
fprintf('\n[quickstart] Building grid...\n');
grid3 = atlas_grid_make(cfg);
fprintf('[quickstart] Grid: %d x %d x %d voxels (%d total)\n', ...
    grid3.Ny, grid3.Nx, grid3.Nth, grid3.Ny * grid3.Nx * grid3.Nth);

%% 3. Build / load family atlases
famA_name = 'Lyapunov L1';    % LL1 in paper Table 2
famB_name = 'Lyapunov L2';    % LL2

fprintf('\n[quickstart] Building atlas for "%s"...\n', famA_name);
t0 = tic;
[SA, ~] = atlas_prepare_or_load(famA_name, cfg, grid3);
fprintf('[quickstart]   Done in %.1f s  (FRS rows: %d)\n', toc(t0), SA.Step4.rows_FRS_upper.n);

fprintf('\n[quickstart] Building atlas for "%s"...\n', famB_name);
t0 = tic;
[SB, ~] = atlas_prepare_or_load(famB_name, cfg, grid3);
fprintf('[quickstart]   Done in %.1f s  (FRS rows: %d)\n', toc(t0), SB.Step4.rows_FRS_upper.n);

%% 4. Compute overlap: FRS(LL1) ∩ BRS(LL2)
fprintf('\n[quickstart] Computing overlap...\n');
t0 = tic;
O = overlap_pair(SA, SB, cfg);
fprintf('[quickstart]   Done in %.1f s\n', toc(t0));
fprintf('[quickstart]   Overlap voxels: %d\n', numel(O.ids));

if numel(O.ids) == 0
    fprintf('[quickstart] No overlap found at this budget/grid. Try increasing\n');
    fprintf('             cfg.fan.DV_cap_nd or cfg.propag.Tmax.\n');
    return
end

%% 5. Extract DV proxy for each overlap voxel
fprintf('\n[quickstart] Extracting voxel metrics...\n');
V = overlap_extract_voxel_info(SA, SB, O, cfg);

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;

% Proxy DV = min turn cost from A + min turn cost from B (already in m/s)
dv_mps   = V.dv_turn_mps_min_A + V.dv_turn_mps_min_B;
tof_days = V.t_days_min_A      + V.t_days_min_B;

[best_dv, k] = min(dv_mps);
best_tof = tof_days(k);

%% 6. Summary
fprintf('\n=== QUICKSTART RESULTS ===\n');
fprintf('Transfer:  %s  →  %s\n', famA_name, famB_name);
fprintf('Grid:      dx=dy=%.3g nd, dtheta=%.1f deg\n', cfg.grid.dx, rad2deg(cfg.grid.dtheta));
fprintf('Budget:    DV_cap=%.3g nd (%.0f m/s),  Tmax=%.4g nd (%.1f days)\n', ...
    cfg.fan.DV_cap_nd, cfg.fan.DV_cap_nd * VU_mps, ...
    cfg.propag.Tmax,   cfg.propag.Tmax   * TU_days);
fprintf('Overlap:   %d voxels\n', numel(O.ids));
fprintf('Best proxy DV:   %.1f m/s\n', best_dv);
fprintf('  at voxel:  (x=%.3f, y=%.3f, th=%.1f deg)\n', ...
    O.x(k), O.y(k), rad2deg(O.th(k)));
fprintf('  proxy TOF: %.2f days\n', best_tof);
fprintf('\nFor paper-quality results, use default grid settings\n');
fprintf('and run scripts/run_overlap_all_pairs.m for all 78 pairs.\n');
fprintf('==========================\n\n');
