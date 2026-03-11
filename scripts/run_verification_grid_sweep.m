%% RUN_VERIFICATION_GRID_SWEEP
% Automated numerical robustness verification sweep.
% Runs 7 core configurations (+ optional Run 7 stress-coarse) with different
% grid / seed / heading-fan parameters to verify that network centrality results
% are stable across reasonable numerical choices.
%
% DV_cap = 0.2 nd  and  Tmax = π  are FIXED across all runs.
% Only the spatial grid (dx, dy), heading grid (dtheta), heading fan (dtheta_fan),
% and seed spacing (ds_seed) are varied.
%
% Run table  (from verification design):
%   Run 0  Baseline           dx=dy=0.0005  dtheta=0.5°  fan=0.5°  ds=0.01
%   Run 1  Coarser spatial    dx=dy=0.001   dtheta=0.5°  fan=0.5°  ds=0.01
%   Run 2  Coarser hdg grid   dx=dy=0.0005  dtheta=1.0°  fan=0.5°  ds=0.01
%   Run 3  Coarser hdg fan    dx=dy=0.0005  dtheta=0.5°  fan=1.0°  ds=0.01
%   Run 4  Coarser seed       dx=dy=0.0005  dtheta=0.5°  fan=0.5°  ds=0.02
%   Run 5  Combined coarse    dx=dy=0.001   dtheta=1.0°  fan=1.0°  ds=0.02
%   Run 6  Combined fine      dx=dy=0.00025 dtheta=0.25° fan=0.25° ds=0.005
%   Run 7* Stress coarse      dx=dy=0.0015  dtheta=1.5°  fan=1.5°  ds=0.03
%            (* only if INCLUDE_RUN7 = true)
%
% Strategy — identical to run_rs4_dv_tmax_sweep.m + run_rs4_all_pairs_summary.m:
%   1. For each run config: build cfg, load 13 atlases, build footprints.
%   2. Slice footprints into FA_arr / FB_arr (one per pair).
%   3. Run 78-pair loop with parfor (sliced footprints — OOM-safe).
%   4. Assemble N×N minDVproxyMat + TOFmat.
%   5. Save to rs4_verification_runs/<RunID>/result.mat  (skip-if-exists).
%
% Output per run:  rs4_verification_runs/<RunID>/result.mat
%   Variables: minDVproxyMat [13×13], TOFmat [13×13], T (winners table),
%              families {13×1}, cfg (run config)
%
% After all runs complete, run run_robustness_analysis.m to compare results.

clear; clc;

% ── repo paths ───────────────────────────────────────────────────────────────
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% ════════════════════════════════════════════════════════════════════════════
%  USER KNOBS
% ════════════════════════════════════════════════════════════════════════════
N_WORKERS    = 4;      % 0 = fully serial; N ≥ 1 = parfor atlas building
INCLUDE_RUN7 = false;  % set true to include the stress-coarse Run 7
INCLUDE_TEST = false;  % set true + edit RunTest_Custom below to add a test run

% ── Fixed parameters (same across all verification runs) ────────────────────
DV_CAP_ND = 0.2;
TMAX_ND   = pi;

% ── Family list (must match your cached atlases) ─────────────────────────────
families = { ...
    'Lyapunov L1', ...
    'Lyapunov L2', ...
    'Cycler 21', ...
    'Cycler 11a', ...
    'Cycler 11b', ...
    'Cycler 32', ...
    'Resonant 2to1 Stable', ...
    'Resonant 2to1 Unstable', ...
    'Resonant 3to1 Stable', ...
    'Resonant 3to1 Unstable', ...
    'Resonant 5to2 Stable', ...
    'Resonant 5to2 Unstable', ...
    'Distant Prograde Orbit' ...
    };

% ════════════════════════════════════════════════════════════════════════════
%  RUN CONFIGURATION TABLE
%  Fields: id, dx, dy, dtheta_deg (heading grid), dtheta_fan_deg (heading fan),
%          ds_seed (seed arc spacing along PO)
% ════════════════════════════════════════════════════════════════════════════
runs(1).id           = 'Run0_Baseline';
runs(1).dx           = 0.0005;
runs(1).dy           = 0.0005;
runs(1).dtheta_deg   = 0.5;
runs(1).dtheta_fan_deg = 0.5;
runs(1).ds_seed      = 0.01;

runs(2).id           = 'Run1_CoarserSpatial';
runs(2).dx           = 0.001;
runs(2).dy           = 0.001;
runs(2).dtheta_deg   = 0.5;
runs(2).dtheta_fan_deg = 0.5;
runs(2).ds_seed      = 0.01;

runs(3).id           = 'Run2_CoarserHdgGrid';
runs(3).dx           = 0.0005;
runs(3).dy           = 0.0005;
runs(3).dtheta_deg   = 1.0;
runs(3).dtheta_fan_deg = 0.5;
runs(3).ds_seed      = 0.01;

runs(4).id           = 'Run3_CoarserHdgFan';
runs(4).dx           = 0.0005;
runs(4).dy           = 0.0005;
runs(4).dtheta_deg   = 0.5;
runs(4).dtheta_fan_deg = 1.0;
runs(4).ds_seed      = 0.01;

runs(5).id           = 'Run4_CoarserSeed';
runs(5).dx           = 0.0005;
runs(5).dy           = 0.0005;
runs(5).dtheta_deg   = 0.5;
runs(5).dtheta_fan_deg = 0.5;
runs(5).ds_seed      = 0.02;

runs(6).id           = 'Run5_CombinedCoarse';
runs(6).dx           = 0.001;
runs(6).dy           = 0.001;
runs(6).dtheta_deg   = 1.0;
runs(6).dtheta_fan_deg = 1.0;
runs(6).ds_seed      = 0.02;

runs(7).id           = 'Run6_CombinedFine';
runs(7).dx           = 0.00025;
runs(7).dy           = 0.00025;
runs(7).dtheta_deg   = 0.25;
runs(7).dtheta_fan_deg = 0.25;
runs(7).ds_seed      = 0.005;

runs(8).id           = 'Run7_StressCoarse';
runs(8).dx           = 0.0015;
runs(8).dy           = 0.0015;
runs(8).dtheta_deg   = 1.5;
runs(8).dtheta_fan_deg = 1.5;
runs(8).ds_seed      = 0.03;

% ── Optional test runs — add/remove entries freely ───────────────────────────
% Each entry needs a unique .id (= output folder name).
% All entries are ignored unless INCLUDE_TEST = true.
test_runs(1).id            = 'RunTest_A';
test_runs(1).dx            = 0.0005;
test_runs(1).dy            = 0.0005;
test_runs(1).dtheta_deg    = 0.5;
test_runs(1).dtheta_fan_deg = 0.5;
test_runs(1).ds_seed       = 0.01;
% Duplicate and extend for more test runs:
% test_runs(2).id            = 'RunTest_B';
% test_runs(2).dx            = ...;  etc.

nRunsTotal = ternary(INCLUDE_RUN7, 8, 7);
if INCLUDE_TEST
    for k = 1:numel(test_runs)
        runs(end+1) = test_runs(k);  %#ok<AGROW>
    end
    nRunsTotal = numel(runs);
end

% ── Output root ──────────────────────────────────────────────────────────────
verifyRoot = fullfile(repoRoot, 'rs4_verification_runs');
if ~exist(verifyRoot, 'dir'), mkdir(verifyRoot); end

% ── Enumerate all 78 pairs (index vectors, computed once) ────────────────────
N      = numel(families);
nPairs = N * (N-1) / 2;
pairI  = zeros(nPairs, 1);
pairJ  = zeros(nPairs, 1);
p = 0;
for ii = 1:N
    for jj = ii+1:N
        p = p + 1;
        pairI(p) = ii;
        pairJ(p) = jj;
    end
end

fprintf('[verify] Verification sweep: %d runs  |  families: %d  |  pairs: %d\n', ...
    nRunsTotal, N, nPairs);
fprintf('[verify] Output root: %s\n\n', verifyRoot);

% ════════════════════════════════════════════════════════════════════════════
%  MAIN LOOP OVER RUN CONFIGURATIONS
% ════════════════════════════════════════════════════════════════════════════
for r = 1:nRunsTotal
    runDef = runs(r);

    outDir     = fullfile(verifyRoot, runDef.id);
    resultFile = fullfile(outDir, 'result.mat');

    if ~exist(outDir, 'dir'), mkdir(outDir); end

    % Skip-if-exists guard (safe to re-run the script)
    if exist(resultFile, 'file')
        fprintf('[verify] Skipping %s — result.mat already exists.\n', runDef.id);
        continue;
    end

    fprintf('\n[verify] ══════ %d/%d  %s ══════\n', r, nRunsTotal, runDef.id);
    fprintf('[verify]   dx=%.5f  dy=%.5f  dtheta=%.3f°  fan=%.3f°  ds=%.4f\n', ...
        runDef.dx, runDef.dy, runDef.dtheta_deg, runDef.dtheta_fan_deg, runDef.ds_seed);
    tRun = tic;

    % ── Build cfg for this run ───────────────────────────────────────────────
    cfg = rs3_cfg_defaults();
    cfg.families.list      = families;
    cfg.families.test_only = false;

    cfg.grid.dx             = runDef.dx;
    cfg.grid.dy             = runDef.dy;
    cfg.grid.dtheta         = deg2rad(runDef.dtheta_deg);
    cfg.seed.ds_seed        = runDef.ds_seed;
    cfg.propag.Tmax         = TMAX_ND;
    cfg.fan.DV_cap_nd       = DV_CAP_ND;
    cfg.fan.dtheta_fan      = deg2rad(runDef.dtheta_fan_deg);
    cfg.propag.absTol       = 1e-8;
    cfg.propag.relTol       = 1e-8;
    cfg.propag.v2tol        = 1e-8;
    cfg.log.step_len_factor = 0.75;
    cfg.log.maxstep_factor  = 2;

    % Cache: enabled by default so atlas builds can be resumed if interrupted.
    % Set cfg.cache.enable = false to skip saving atlas cache files to disk.
    cfg.cache.enable  = true;
    cfg.cache.rebuild = false;

    cfg.io.save_figs   = false;
    cfg.io.save_fig    = false;
    cfg.io.fig_visible = 'off';

    cfg.plot.rs4.overlap_xy   = false;
    cfg.plot.rs4.overlap_xyz  = false;
    cfg.plot.rs4.combo_xy     = false;
    cfg.plot.rs4.combo_xyz    = false;
    cfg.plot.rs4.bounds_lb    = false;
    cfg.plot.rs4.bounds_ub    = false;
    cfg.plot.rs4.bounds_proxy = false;

    if exist('rs3_cfg_validate', 'file') == 2
        rs3_cfg_validate(cfg);
    end

    VU_mps  = cfg.units.VU_mps;
    TU_days = cfg.units.TU_days;

    % ── Build grid ───────────────────────────────────────────────────────────
    grid3 = rs3_grid_make(cfg);

    % ── Create/reuse pool right before atlas parfor ──────────────────────────
    % Pool is started here (fresh per run) so workers are never left idle during
    % long atlas builds on an HPC scheduler that reclaims idle workers.
    if N_WORKERS > 0
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('local', N_WORKERS);
            fprintf('[verify] Started parpool with %d workers.\n', N_WORKERS);
        else
            fprintf('[verify] Reusing parpool (%d workers).\n', pool.NumWorkers);
        end
    end

    % ── Load atlas + build footprint one family at a time ───────────────────
    % Atlas and footprint loops are merged so only 1 atlas is in RAM at a time.
    % For fine grids (32x more voxels) holding all 13 atlases simultaneously
    % would exceed available RAM — clear S immediately after footprint is built.
    cfg.par.enable = (N_WORKERS > 0);
    fprintf('[verify] Loading/building %d atlases → footprints (1 atlas in RAM at a time)...\n', N);
    Fall = cell(N, 1);
    for i = 1:N
        fprintf('[verify]   %2d/%d  %s\n', i, N, families{i});
        [S, ~] = rs3_prepare_or_load_family(families{i}, cfg, grid3);
        Fall{i} = local_compute_footprint(S, grid3, VU_mps, TU_days);
        clear S;
    end
    fprintf('[verify] All footprints built.\n');

    % ── Pair intersection loop (no FA_arr/FB_arr copy — access Fall directly) ──
    % Precompute the last pair index each family appears in, so footprints can
    % be freed the moment their last pair is processed.  For N=13 this means
    % family 1 is freed after pair (1,13), family 2 after (2,13), etc., so
    % peak RAM during the pair loop is at most 13→1 footprints (monotonically
    % decreasing) rather than a flat 13 throughout.
    lastUse = zeros(N, 1);
    for p = 1:nPairs
        lastUse(pairI(p)) = p;
        lastUse(pairJ(p)) = p;
    end

    pair_minDV   = nan(nPairs, 1);
    pair_DVlb    = nan(nPairs, 1);
    pair_DVpatch = nan(nPairs, 1);
    pair_TOF     = nan(nPairs, 1);
    pair_voxelId = nan(nPairs, 1);

    cfg_p   = cfg;
    grid3_p = grid3;
    VU_p    = VU_mps;

    fprintf('[verify] Running %d pair intersections (serial)...\n', nPairs);
    tPairs = tic;

    for p = 1:nPairs
        [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
            pair_TOF(p), pair_voxelId(p)] = ...
            local_run_pair(Fall{pairI(p)}, Fall{pairJ(p)}, grid3_p, cfg_p, VU_p);
        % Free footprints as soon as all their pairs have been processed
        if lastUse(pairI(p)) == p, Fall{pairI(p)} = []; end
        if lastUse(pairJ(p)) == p, Fall{pairJ(p)} = []; end
    end
    clear Fall;

    n_overlap = sum(isfinite(pair_minDV));
    fprintf('[verify] Pairs done in %.1f s — %d/%d have overlap.\n', ...
        toc(tPairs), n_overlap, nPairs);

    % ── Assemble N×N matrices ─────────────────────────────────────────────────
    minDVproxyMat = nan(N, N);
    TOFmat        = nan(N, N);
    for p = 1:nPairs
        i = pairI(p);  j = pairJ(p);
        minDVproxyMat(i, j) = pair_minDV(p);
        minDVproxyMat(j, i) = pair_minDV(p);
        TOFmat(i, j)        = pair_TOF(p);
        TOFmat(j, i)        = pair_TOF(p);
    end

    % ── Winners table ─────────────────────────────────────────────────────────
    pair_famA = families(pairI);
    pair_famB = families(pairJ);
    T = table(pair_famA(:), pair_famB(:), ...
        pair_minDV, pair_DVlb, pair_DVpatch, pair_TOF, pair_voxelId, ...
        'VariableNames', { ...
            'FamilyA', 'FamilyB', ...
            'minDVproxy_mps', 'DVlb_mps', 'DVpatch_ub_mps', ...
            'EstimatedTOF_days', 'VoxelId'});
    T = sortrows(T, 'minDVproxy_mps', 'MissingPlacement', 'last');

    % ── Save result ───────────────────────────────────────────────────────────
    save(resultFile, ...
        'minDVproxyMat', 'TOFmat', 'T', 'families', 'cfg', '-v7');

    fprintf('[verify] %s  saved in %.1f s  (%.0f/%.0f pairs with overlap)\n', ...
        runDef.id, toc(tRun), n_overlap, nPairs);
    fprintf('[verify] Output: %s\n', resultFile);

    % Free memory before next run
    clear grid3 minDVproxyMat TOFmat T lastUse;
    clear pair_minDV pair_DVlb pair_DVpatch pair_TOF pair_voxelId;

end  % for r

fprintf('\n[verify] ══════════ ALL RUNS COMPLETE ══════════\n');
fprintf('[verify] Output root: %s\n', verifyRoot);
fprintf('[verify] Next step  : run run_robustness_analysis.m\n');

% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  (identical to run_rs4_dv_tmax_sweep.m / run_rs4_all_pairs_summary.m)
% ════════════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────────────
function [minDV, dvlb, dvpatch, tof, voxelId] = ...
        local_run_pair(FA, FB, grid3, cfg, VU_mps)
minDV = NaN;  dvlb = NaN;  dvpatch = NaN;  tof = NaN;  voxelId = NaN;
try
    idsO = intersect(FA.uid_frs, FB.uid_brs);
    if isempty(idsO), return; end

    Ny = numel(grid3.y_centers);
    Nx = numel(grid3.x_centers);
    Nt = numel(grid3.th_centers);
    [iy, ix, ~] = ind2sub([Ny, Nx, Nt], idsO);

    bufFrac = 0.05;
    if isfield(cfg,'overlap') && isfield(cfg.overlap,'primary_buffer_frac') ...
            && ~isempty(cfg.overlap.primary_buffer_frac)
        bufFrac = cfg.overlap.primary_buffer_frac;
    end
    if ~(isfield(cfg,'sys') && isfield(cfg.sys,'RE_nd') && isfield(cfg.sys,'RM_nd'))
        error('cfg.sys.RE_nd and cfg.sys.RM_nd required for primary buffer filter.');
    end
    RE = cfg.sys.RE_nd;
    RM = cfg.sys.RM_nd;
    mu = FA.mu;

    if isfield(grid3,'Keep') && ~isempty(grid3.Keep)
        keepXY = logical(grid3.Keep);
        if ~isequal(size(keepXY), [Ny, Nx]), keepXY = keepXY.'; end
        okKeep = keepXY(sub2ind([Ny, Nx], iy, ix));
    else
        okKeep = true(numel(idsO), 1);
    end

    x = grid3.x_centers(ix);
    y = grid3.y_centers(iy);
    okEarth = hypot(x + mu, y)     > (1 + bufFrac) * RE;
    okMoon  = hypot(x - (1-mu), y) > (1 + bufFrac) * RM;
    ok = okKeep(:) & okEarth(:) & okMoon(:);

    idsO = idsO(ok);
    if isempty(idsO), return; end
    ix = ix(ok);  iy = iy(ok);

    [~, locA] = ismember(idsO, FA.uid_frs);
    [~, locB] = ismember(idsO, FB.uid_brs);
    dv_min_A = FA.dv_min_frs(locA);
    dv_min_B = FB.dv_min_brs(locB);
    t_mean_A = FA.t_mean_frs(locA);
    t_mean_B = FB.t_mean_brs(locB);

    x_ok = grid3.x_centers(ix);
    y_ok = grid3.y_centers(iy);
    CJstar = min(FA.CJ, FB.CJ);
    pot    = rs3_core_cr3bp_U_and_derivs(x_ok(:), y_ok(:), mu);
    v_box  = sqrt(max(2*pot.U - CJstar, 0));
    dv_patch_vec = 2 * v_box .* sin(abs(grid3.dtheta) / 2) * VU_mps;
    dv_lb_vec    = dv_min_A(:) + dv_min_B(:);
    dv_proxy     = dv_lb_vec + dv_patch_vec;

    valid = isfinite(dv_proxy);
    if ~any(valid), return; end

    idxValid  = find(valid);
    [~, iLoc] = min(dv_proxy(idxValid));
    iWin      = idxValid(iLoc);

    minDV   = dv_proxy(iWin);
    dvlb    = dv_lb_vec(iWin);
    dvpatch = dv_patch_vec(iWin);
    voxelId = idsO(iWin);
    tof     = t_mean_A(iWin) + t_mean_B(iWin);
catch ME
    warning('[verify:pair] %s', ME.message);
end
end

% ─────────────────────────────────────────────────────────────────────────────
function F = local_compute_footprint(S, grid3, VU_mps, TU_days)
Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

thm    = rs3_wrapToPi(pi - grid3.th_centers(:));
lut    = discretize(thm, grid3.th_edges);
lut(isnan(lut)) = 0;
it_lut = uint16(lut);

dlists = S.Step4.delta_lists;
Ns     = numel(dlists);
max_h  = max(1, max(cellfun(@numel, dlists)));
delta_mat = zeros(Ns, max_h);
for s = 1:Ns
    v = double(dlists{s});
    delta_mat(s, 1:numel(v)) = v;
end

pot_u    = rs3_core_cr3bp_U_and_derivs(S.SeedsUpper(:,1), S.SeedsUpper(:,2), S.mu);
v0_upper = sqrt(max(2*pot_u.U(:) - S.CJ, 0));
if isfield(S,'SeedsLower') && ~isempty(S.SeedsLower)
    pot_l    = rs3_core_cr3bp_U_and_derivs(S.SeedsLower(:,1), S.SeedsLower(:,2), S.mu);
    v0_lower = sqrt(max(2*pot_l.U(:) - S.CJ, 0));
else
    v0_lower = v0_upper;
end

% ── Upper FRS rows ──────────────────────────────────────────────────────────
nu = double(S.Step4.rows_FRS_upper.n);
if nu > 0
    [ids_u, dv_u, t_u, ix_u, iy_u, it_u] = local_fp_rows( ...
        S.Step4.rows_FRS_upper, nu, v0_upper, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
else
    ids_u=zeros(0,1); dv_u=zeros(0,1); t_u=zeros(0,1);
    ix_u=zeros(0,1);           iy_u=zeros(0,1); it_u=zeros(0,1);
end

% BRS upper — derived immediately so ix_u/iy_u/it_u can be freed before lower rows
if ~isempty(ix_u)
    biy_u = Ny - iy_u + 1;
    bit_u = double(it_lut(it_u));
    ok_u  = bit_u > 0;
    ids_brs_u = sub2ind([Ny,Nx,Nt], biy_u(ok_u), ix_u(ok_u), max(1,min(Nt,bit_u(ok_u))));
    dv_bu = dv_u(ok_u);  t_bu = t_u(ok_u);
else
    ids_brs_u = zeros(0,1);  dv_bu = zeros(0,1);  t_bu = zeros(0,1);
end
clear ix_u iy_u it_u;

% ── Lower FRS rows ──────────────────────────────────────────────────────────
nl = double(S.Step4.rows_FRS_lower.n);
if nl > 0
    [ids_l, dv_l, t_l, ix_l, iy_l, it_l] = local_fp_rows( ...
        S.Step4.rows_FRS_lower, nl, v0_lower, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
else
    ids_l=zeros(0,1); dv_l=zeros(0,1); t_l=zeros(0,1);
    ix_l=zeros(0,1);           iy_l=zeros(0,1); it_l=zeros(0,1);
end

% BRS lower — same pattern: derive then free spatial index arrays
if ~isempty(ix_l)
    biy_l = Ny - iy_l + 1;
    bit_l = double(it_lut(it_l));
    ok_l  = bit_l > 0;
    ids_brs_l = sub2ind([Ny,Nx,Nt], biy_l(ok_l), ix_l(ok_l), max(1,min(Nt,bit_l(ok_l))));
    dv_bl = dv_l(ok_l);  t_bl = t_l(ok_l);
else
    ids_brs_l = zeros(0,1);  dv_bl = zeros(0,1);  t_bl = zeros(0,1);
end
clear ix_l iy_l it_l;

% ── Aggregate FRS, then free raw rows ──────────────────────────────────────
[F.uid_frs, F.dv_min_frs, F.t_mean_frs] = local_fp_agg( ...
    [ids_u; ids_l], [dv_u; dv_l], [t_u; t_l]);
clear ids_u dv_u t_u ids_l dv_l t_l;

% ── Aggregate BRS ──────────────────────────────────────────────────────────
[F.uid_brs, F.dv_min_brs, F.t_mean_brs] = local_fp_agg( ...
    [ids_brs_u; ids_brs_l], [dv_bu; dv_bl], [t_bu; t_bl]);

F.CJ   = S.CJ;
F.mu   = S.mu;
F.name = S.name;
end

% ─────────────────────────────────────────────────────────────────────────────
function [ids, dv_mps, t_days, ix_out, iy_out, it_out] = local_fp_rows( ...
        rows, n, v0_per_seed, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days)
ix_out = double(rows.ix(1:n));
iy_out = double(rows.iy(1:n));
it_out = double(rows.it(1:n));
ids    = sub2ind([Ny, Nx, Nt], iy_out, ix_out, it_out);
iSeed  = double(rows.iSeed(1:n));
iHead  = double(rows.iHead(1:n));
t_nd   = double(rows.t(1:n));
lin    = sub2ind([Ns, max_h], iSeed, iHead);
delta  = delta_mat(lin);
v0     = v0_per_seed(iSeed);
dv_mps = 2 * v0(:) .* sin(abs(delta(:)) / 2) * VU_mps;
t_days = abs(t_nd(:)) * TU_days;
end

% ─────────────────────────────────────────────────────────────────────────────
function [uid, dv_min, t_mean] = local_fp_agg(ids, dv, t)
if isempty(ids)
    uid = zeros(0,1);  dv_min = zeros(0,1);  t_mean = zeros(0,1);
    return;
end
[uid, ~, ic] = unique(ids(:));
dv_min = accumarray(ic, dv(:), [], @min);
t_mean = accumarray(ic, t(:)) ./ accumarray(ic, ones(numel(ic), 1));
end

% ─────────────────────────────────────────────────────────────────────────────
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
