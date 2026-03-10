%% RUN_RS4_ALL_PAIRS_SUMMARY
% Single-configuration DV-proxy batch runner over all 13×12/2 = 78 family pairs.
%
% REFACTORED to use the same footprint + parfor strategy as
% run_rs4_dv_tmax_sweep.m (OOM-safe, no full atlas broadcast).
%
% Strategy:
%   1. Load all 13 base atlases once (serial, main process — Java MD5 safe).
%   2. Build compact per-family footprints (~5-25 MB each).
%   3. Slice footprints into FA_arr / FB_arr (one entry per pair).
%   4. Run pair loop: parfor over 78 pairs using sliced footprints.
%   5. Assemble N×N minDVproxyMat + TOFmat and save outputs.
%
% Outputs written to <out_root>/<tag>/rs4_pairs_13fam/Summary/:
%   minDVproxy_matrix.csv           — N×N DV matrix (m/s)
%   pair_winners_top1.csv           — long-form winners table
%   batch_summary_workspace.mat     — minDVproxyMat, T, cfg, families
%
% For per-pair figures and detailed overlap visualisation use
%   run_rs4_overlap_and_visuals.m  (separate heavy pipeline).

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
N_WORKERS = 4;   % 0 = fully serial; N ≥ 1 = parfor with N workers

% POOL_AFTER_ATLAS — when to start the parallel pool:
%   true  (default/safe): pool starts AFTER atlas loading, right before
%         footprint parfor.  Safe for fresh atlas builds (can take hours).
%   false : pool starts before atlas loading (only safe when all atlases
%         are cached and idle-worker timeout is not a concern).
POOL_AFTER_ATLAS = true;

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
%  CFG  (must match the cached atlases)
% ════════════════════════════════════════════════════════════════════════════
cfg = rs3_cfg_defaults();

cfg.families.list      = families;
cfg.families.test_only = false;

cfg.grid.dx               = 0.0005;
cfg.grid.dy               = 0.0005;
cfg.grid.dtheta           = deg2rad(0.5);
cfg.seed.ds_seed          = 0.01;
cfg.propag.Tmax           = pi;
cfg.fan.DV_cap_nd         = 0.2;
cfg.fan.dtheta_fan        = deg2rad(0.5);
cfg.propag.absTol         = 1e-8;
cfg.propag.relTol         = 1e-8;
cfg.propag.v2tol          = 1e-8;
cfg.log.step_len_factor   = 0.75;
cfg.log.maxstep_factor    = 2;

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

% ── Derived constants ────────────────────────────────────────────────────────
VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
N       = numel(families);

% ── Enumerate pairs ──────────────────────────────────────────────────────────
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

% ── Output directories ───────────────────────────────────────────────────────
batchTag   = sprintf('rs4_pairs_%dfam', N);
outRoot    = fullfile(cfg.io.out_root, cfg.io.tag, batchTag);
summaryDir = fullfile(outRoot, 'Summary');
if ~exist(summaryDir, 'dir'), mkdir(summaryDir); end

fprintf('[rs4-batch] Output: %s\n', summaryDir);
fprintf('[rs4-batch] Mode  : %s\n\n', ...
    ternary(N_WORKERS > 0, sprintf('PARALLEL (%d workers)', N_WORKERS), 'SERIAL'));

% ════════════════════════════════════════════════════════════════════════════
%  1. LOAD ALL FAMILY ATLASES  (serial, main process — Java MD5 required)
% ════════════════════════════════════════════════════════════════════════════
if ~POOL_AFTER_ATLAS && N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[rs4-batch] Started parpool with %d workers.\n', N_WORKERS);
    else
        fprintf('[rs4-batch] Using existing parpool (%d workers).\n', pool.NumWorkers);
    end
end

fprintf('[rs4-batch] Loading %d family atlases...\n', N);
grid3 = rs3_grid_make(cfg);
Sall  = cell(N, 1);
for i = 1:N
    fprintf('[rs4-batch]   %d/%d  %s\n', i, N, families{i});
    [Sall{i}, ~] = rs3_prepare_or_load_family(families{i}, cfg, grid3);
end
fprintf('[rs4-batch] Atlases loaded.\n\n');

% ════════════════════════════════════════════════════════════════════════════
%  2. BUILD COMPACT FOOTPRINTS  (~5-25 MB each vs ~0.5-2 GB full atlas)
% ════════════════════════════════════════════════════════════════════════════
if POOL_AFTER_ATLAS && N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[rs4-batch] Started parpool with %d workers.\n', N_WORKERS);
    else
        fprintf('[rs4-batch] Using existing parpool (%d workers).\n', pool.NumWorkers);
    end
end
fprintf('[rs4-batch] Building voxel footprints...\n');
Fall = cell(N, 1);
if N_WORKERS > 0
    parfor i = 1:N
        Fall{i} = local_compute_footprint(Sall{i}, grid3, VU_mps, TU_days);
    end
else
    for i = 1:N
        Fall{i} = local_compute_footprint(Sall{i}, grid3, VU_mps, TU_days);
    end
end
clear Sall;
fprintf('[rs4-batch] Footprints built. Full atlases released.\n\n');

% ════════════════════════════════════════════════════════════════════════════
%  3. SLICE INTO PER-PAIR ARRAYS  (required for parfor slicing)
% ════════════════════════════════════════════════════════════════════════════
FA_arr = cell(nPairs, 1);
FB_arr = cell(nPairs, 1);
for p = 1:nPairs
    FA_arr{p} = Fall{pairI(p)};
    FB_arr{p} = Fall{pairJ(p)};
end
clear Fall;

% ════════════════════════════════════════════════════════════════════════════
%  4. PAIR LOOP  (parfor over 78 pairs — sliced footprints only)
% ════════════════════════════════════════════════════════════════════════════
fprintf('[rs4-batch] Running %d pair intersections...\n', nPairs);
tPairs = tic;

pair_minDV   = nan(nPairs, 1);
pair_DVlb    = nan(nPairs, 1);
pair_DVpatch = nan(nPairs, 1);
pair_TOF     = nan(nPairs, 1);
pair_voxelId = nan(nPairs, 1);

cfg_p   = cfg;
grid3_p = grid3;
VU_p    = VU_mps;

if N_WORKERS > 0
    parfor p = 1:nPairs
        [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
            pair_TOF(p), pair_voxelId(p)] = ...
            local_run_pair(FA_arr{p}, FB_arr{p}, grid3_p, cfg_p, VU_p);
    end
else
    for p = 1:nPairs
        fprintf('[rs4-batch]   pair %d/%d: %-28s → %s\n', ...
            p, nPairs, families{pairI(p)}, families{pairJ(p)});
        [pair_minDV(p), pair_DVlb(p), pair_DVpatch(p), ...
            pair_TOF(p), pair_voxelId(p)] = ...
            local_run_pair(FA_arr{p}, FB_arr{p}, grid3_p, cfg_p, VU_p);
    end
end
clear FA_arr FB_arr;

n_overlap = sum(isfinite(pair_minDV));
fprintf('[rs4-batch] Pair loop done in %.1f s — %d/%d pairs have overlap.\n\n', ...
    toc(tPairs), n_overlap, nPairs);

% ════════════════════════════════════════════════════════════════════════════
%  5. ASSEMBLE N×N MATRICES
% ════════════════════════════════════════════════════════════════════════════
minDVproxyMat = nan(N, N);
TOFmat        = nan(N, N);
for p = 1:nPairs
    i = pairI(p);  j = pairJ(p);
    minDVproxyMat(i, j) = pair_minDV(p);
    minDVproxyMat(j, i) = pair_minDV(p);
    TOFmat(i, j)        = pair_TOF(p);
    TOFmat(j, i)        = pair_TOF(p);
end

% ── Winners table (same format as before — compatible with tmax sweep scanner) ──
pair_famA = families(pairI);
pair_famB = families(pairJ);
T = table(pair_famA(:), pair_famB(:), ...
    pair_minDV, pair_DVlb, pair_DVpatch, pair_TOF, pair_voxelId, ...
    'VariableNames', { ...
        'FamilyA', 'FamilyB', ...
        'minDVproxy_mps', 'DVlb_mps', 'DVpatch_ub_mps', ...
        'EstimatedTOF_days', 'VoxelId'});
T = sortrows(T, 'minDVproxy_mps', 'MissingPlacement', 'last');

% ════════════════════════════════════════════════════════════════════════════
%  6. PRINT SUMMARY
% ════════════════════════════════════════════════════════════════════════════
fprintf('[rs4-batch] ===== DV PROXY SUMMARY =====\n');
fprintf('  Total pairs        : %d\n', nPairs);
fprintf('  Pairs with overlap : %d\n', n_overlap);
fprintf('  Pairs no overlap   : %d\n', nPairs - n_overlap);
if n_overlap > 0
    valid_dv = pair_minDV(isfinite(pair_minDV));
    fprintf('  Min proxy DV       : %.1f m/s\n', min(valid_dv));
    fprintf('  Max proxy DV       : %.1f m/s\n', max(valid_dv));
    fprintf('  Mean proxy DV      : %.1f m/s\n', mean(valid_dv));
    fprintf('\n  Top-5 pairs by proxy DV:\n');
    top5 = min(5, n_overlap);
    for k = 1:top5
        fprintf('    %2d. %-26s → %-26s  %.1f m/s\n', ...
            k, T.FamilyA{k}, T.FamilyB{k}, T.minDVproxy_mps(k));
    end
end

% ════════════════════════════════════════════════════════════════════════════
%  7. SAVE OUTPUTS
% ════════════════════════════════════════════════════════════════════════════
% CSV 1: N×N matrix
matrixCell = cell(N+1, N+1);
matrixCell{1,1} = 'Family';
for i = 1:N
    matrixCell{1, i+1} = families{i};
    matrixCell{i+1, 1} = families{i};
    for j = 1:N
        matrixCell{i+1, j+1} = minDVproxyMat(i,j);
    end
end
writecell(matrixCell, fullfile(summaryDir, 'minDVproxy_matrix.csv'));

% CSV 2: long-form winners
writetable(T, fullfile(summaryDir, 'pair_winners_top1.csv'));

% MAT workspace (compatible with run_rs4_dv_tmax_sweep.m scanner)
save(fullfile(summaryDir, 'batch_summary_workspace.mat'), ...
    'families', 'minDVproxyMat', 'TOFmat', 'T', 'cfg', '-v7.3');

fprintf('\n[rs4-batch] Done.\n');
fprintf('  Summary dir: %s\n', summaryDir);

% ════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  (identical to run_rs4_dv_tmax_sweep.m)
% ════════════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────────────
function [minDV, dvlb, dvpatch, tof, voxelId] = ...
        local_run_pair(FA, FB, grid3, cfg, VU_mps)
%LOCAL_RUN_PAIR  Compute DV proxy and TOF for one pair using footprints.
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
    dv_patch_vec = 2*v_box .* sin(abs(grid3.dtheta)/2) * VU_mps;
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
    warning('[rs4-batch:pair] %s', ME.message);
end
end

% ─────────────────────────────────────────────────────────────────────────────
function F = local_compute_footprint(S, grid3, VU_mps, TU_days)
%LOCAL_COMPUTE_FOOTPRINT  Build compact per-voxel summary for one atlas family.
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

nu = double(S.Step4.rows_FRS_upper.n);
if nu > 0
    [ids_u, dv_u, t_u, ix_u, iy_u, it_u] = local_fp_rows( ...
        S.Step4.rows_FRS_upper, nu, v0_upper, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
else
    ids_u=zeros(0,1); dv_u=zeros(0,1); t_u=zeros(0,1);
    ix_u=zeros(0,1);  iy_u=zeros(0,1); it_u=zeros(0,1);
end

nl = double(S.Step4.rows_FRS_lower.n);
if nl > 0
    [ids_l, dv_l, t_l, ix_l, iy_l, it_l] = local_fp_rows( ...
        S.Step4.rows_FRS_lower, nl, v0_lower, delta_mat, Ns, max_h, Ny, Nx, Nt, VU_mps, TU_days);
else
    ids_l=zeros(0,1); dv_l=zeros(0,1); t_l=zeros(0,1);
    ix_l=zeros(0,1);  iy_l=zeros(0,1); it_l=zeros(0,1);
end

[F.uid_frs, F.dv_min_frs, F.t_mean_frs] = local_fp_agg( ...
    [ids_u; ids_l], [dv_u; dv_l], [t_u; t_l]);

if ~isempty(ix_u)
    biy_u = Ny - iy_u + 1;
    bit_u = double(it_lut(it_u));
    ok_u  = bit_u > 0;
    ids_brs_u = sub2ind([Ny,Nx,Nt], biy_u(ok_u), ix_u(ok_u), max(1,min(Nt,bit_u(ok_u))));
    dv_bu = dv_u(ok_u);  t_bu = t_u(ok_u);
else
    ids_brs_u = zeros(0,1);  dv_bu = zeros(0,1);  t_bu = zeros(0,1);
end

if ~isempty(ix_l)
    biy_l = Ny - iy_l + 1;
    bit_l = double(it_lut(it_l));
    ok_l  = bit_l > 0;
    ids_brs_l = sub2ind([Ny,Nx,Nt], biy_l(ok_l), ix_l(ok_l), max(1,min(Nt,bit_l(ok_l))));
    dv_bl = dv_l(ok_l);  t_bl = t_l(ok_l);
else
    ids_brs_l = zeros(0,1);  dv_bl = zeros(0,1);  t_bl = zeros(0,1);
end

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
