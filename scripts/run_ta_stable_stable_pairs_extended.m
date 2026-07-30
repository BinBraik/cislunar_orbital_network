%% RUN_TA_STABLE_STABLE_PAIRS_EXTENDED
% Extends the stable-stable pair test (R21-S, R31-S, R52-S vs each
% other) from 32pi/24 rungs to 64pi/32 rungs. All three pairs came in
% notably ABOVE the naive "sum of individual floors" prediction on the
% first pass, and none had clearly finished descending (R21-S<->R31-S
% and R21-S<->R52-S were still monotonically decreasing at the last
% rung; only R31-S<->R52-S looked like a firm plateau). The 8 new rungs
% are packed densely into the NEW stretch (32pi-64pi, at double the
% density of the original ladder) rather than just continuing the same
% spacing, specifically to catch lobe-switch voxel jumps -- the winning
% connection point can jump discontinuously to a cheaper voxel as Ta
% grows, and a coarser ladder can miss exactly where that happens.
%
% WHY A FRESH ATLAS BUILD IS UNAVOIDABLE: atlas_derive_subset only
% shrinks Tmax from a cached atlas, never expands it -- so the existing
% 32pi atlases for all three families can't be reused for this, even
% though 24 of the 32 new rungs represent identical Ta values as before.
% Each family needs a genuinely longer (Tmax=64pi) integration. The
% idempotent per-family cache check also validates against the WHOLE
% ladder array, so all 32 rungs get their (cheap) derive+save redone
% even for the previously-known Ta values -- that's a small, acceptable
% cost next to the one unavoidable expensive part (the fresh 64pi
% integration itself).
%
% dv_A/dv_B SPLIT NOW NATIVE TO PHASE B: local_run_pair already computed
% dv_min_A and dv_min_B internally before summing them into DVlb_mps --
% it just never returned them. Now both are written directly into the
% Phase B results table for every rung, not just the largest-Ta rung via
% a separate reconstruction pass.
%
% Same proven shape otherwise: single MATLAB process, single parpool,
% serial family loop, save-per-rung-immediately, incremental writes for
% both Phase B results and the trajectory reconstruction that follows.

clear; clc;

% ── repo paths ──────────────────────────────────────────────────────────────
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% ══════════════════════════════════════════════════════════════════════════
%  USER KNOBS
% ══════════════════════════════════════════════════════════════════════════

FAMILIES3 = {'Resonant 2to1 Stable', 'Resonant 3to1 Stable', 'Resonant 5to2 Stable'};
PAIRS = {1, 2; 1, 3; 2, 3};   % (iA, iB) index pairs into FAMILIES3

% Original 24 rungs (0.354pi..32pi) UNCHANGED, plus 8 new rungs packed
% densely into the new 32pi..64pi stretch (double the original density).
Ta_old    = 2.^linspace(-1.5, 5, 24);
Ta_newtail = 2.^(5 + (1:8)/8);
Ta_multiples_of_pi = sort([Ta_old, Ta_newtail], 'descend');
DV_cap_nd = 0.15;

N_WORKERS = 60;
envWorkers = getenv('STABLE_STABLE_EXT_NWORKERS');
if ~isempty(envWorkers)
    N_WORKERS = str2double(envWorkers);
end

OUTPUT_DIR    = fullfile(repoRoot, 'ta_stable_stable_extended_results');
BYRUNG_DIR    = fullfile(OUTPUT_DIR, 'footprints_by_rung');
RESULTS_DIR   = fullfile(repoRoot, 'ta_stable_stable_extended_final');
TRAJ_OUT_DIR  = fullfile(repoRoot, 'ta_stable_stable_extended_trajectories');
if ~exist(BYRUNG_DIR, 'dir'), mkdir(BYRUNG_DIR); end
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end
if ~exist(TRAJ_OUT_DIR, 'dir'), mkdir(TRAJ_OUT_DIR); end

% ══════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION
% ══════════════════════════════════════════════════════════════════════════
cfg = atlas_cfg_defaults();
cfg.grid.dx               = 0.001;
cfg.grid.dy               = 0.001;
cfg.grid.dtheta           = deg2rad(1);
cfg.seed.ds_seed          = 0.01;
cfg.fan.DV_cap_nd         = DV_cap_nd;
cfg.fan.dtheta_fan        = deg2rad(0.5);
cfg.propag.absTol         = 1e-8;
cfg.propag.relTol         = 1e-8;
cfg.propag.v2tol          = 1e-8;
cfg.log.step_len_factor   = 0.75;
cfg.log.maxstep_factor    = 2;
cfg.par.progress_every    = 1000;

cfg.cache.enable  = true;
cfg.cache.dir     = fullfile(repoRoot, 'atlas_cache');
cfg.cache.rebuild = false;

cfg.par.enable = (N_WORKERS > 0);

cfg.io.save_figs   = false;
cfg.io.save_fig    = false;
cfg.io.fig_visible = 'off';

if exist('atlas_cfg_validate', 'file') == 2
    atlas_cfg_validate(cfg);
end

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;

grid3 = atlas_grid_make(cfg);
nFam  = numel(FAMILIES3);
nRung = numel(Ta_multiples_of_pi);
nPairs = size(PAIRS, 1);

fprintf('[stable_stable_ext] Families: %d   Ta rungs: %d (%.1f .. %.1f days)   DV_cap_nd: %.3g\n', ...
    nFam, nRung, Ta_multiples_of_pi(end)*pi*TU_days, Ta_multiples_of_pi(1)*pi*TU_days, DV_cap_nd);

if N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[stable_stable_ext] Started parpool with %d workers.\n', N_WORKERS);
    end
end

% ══════════════════════════════════════════════════════════════════════════
%  PHASE A -- per-family footprints, serial, save-per-rung-immediately
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[stable_stable_ext] ══════════ PHASE A: per-family footprints ══════════\n');

for i = 1:nFam
    familyName = FAMILIES3{i};
    key    = local_fieldkey(familyName);
    famDir = fullfile(BYRUNG_DIR, key);
    if ~exist(famDir, 'dir'), mkdir(famDir); end

    allDone = true;
    for r = 1:nRung
        rpath = fullfile(famDir, sprintf('rung%02d.mat', r));
        if ~isfile(rpath), allDone = false; break; end
        ck = load(rpath, 'Ta_multiples_of_pi', 'DV_cap_nd');
        if ~isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) || ck.DV_cap_nd ~= DV_cap_nd
            allDone = false; break;
        end
    end
    if allDone
        fprintf('[stable_stable_ext] (%d/%d) %-25s -- all %d rungs already cached, skipping.\n', ...
            i, nFam, familyName, nRung);
        continue;
    end

    tFam = tic;
    Tmax_fat = Ta_multiples_of_pi(1) * pi;
    fprintf('[stable_stable_ext] (%d/%d) %-25s -- building fat atlas at Ta=%.4gpi (%.2f days)...\n', ...
        i, nFam, familyName, Ta_multiples_of_pi(1), Tmax_fat * TU_days);

    cfg_fat = cfg;
    cfg_fat.propag.Tmax = Tmax_fat;
    [S_fat, cacheInfo] = atlas_prepare_or_load(familyName, cfg_fat, grid3);
    fprintf('[stable_stable_ext]   fat atlas %s in %.1fs.\n', local_hitstr(cacheInfo), toc(tFam));

    for r = 1:nRung
        rpath = fullfile(famDir, sprintf('rung%02d.mat', r));
        if isfile(rpath)
            ck = load(rpath, 'Ta_multiples_of_pi', 'DV_cap_nd');
            if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) && ck.DV_cap_nd == DV_cap_nd
                fprintf('[stable_stable_ext]   %-25s rung %d/%d already cached, skipping.\n', familyName, r, nRung);
                continue;
            end
        end

        Tmax_r = Ta_multiples_of_pi(r) * pi;
        if r == 1
            Ssub = S_fat;   % top rung IS the fat atlas, no filtering needed
        else
            cfg_r = cfg;
            cfg_r.propag.Tmax = Tmax_r;
            Ssub = atlas_derive_subset(S_fat, cfg_r);
        end
        F = local_compute_footprint(Ssub, grid3, VU_mps, TU_days); %#ok<NASGU>
        clear Ssub
        if ~exist(famDir, 'dir'), mkdir(famDir); end   % defensive re-check before every save
        save(rpath, 'F', 'Ta_multiples_of_pi', 'DV_cap_nd', '-v7.3');
        clear F
        fprintf('[stable_stable_ext]   %-25s rung %d/%d (Ta=%.2fd) done and saved (%.1fs elapsed).\n', ...
            familyName, r, nRung, Tmax_r * TU_days, toc(tFam));
    end

    clear S_fat
    fprintf('[stable_stable_ext] (%d/%d) %-25s -- done in %.1fs total.\n', i, nFam, familyName, toc(tFam));
end

fprintf('\n[stable_stable_ext] ══════════ Phase A complete. ══════════\n');

% ══════════════════════════════════════════════════════════════════════════
%  PHASE B -- pairwise computation, dv_A/dv_B split now native
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[stable_stable_ext] ══════════ PHASE B: %d pairs x %d rungs ══════════\n', nPairs, nRung);

results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVA_mps', {}, 'DVB_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});

resultsCsv = fullfile(RESULTS_DIR, 'ta_stable_stable_extended_results.csv');

for r = 1:nRung
    Tmax_r = Ta_multiples_of_pi(r) * pi;
    tRung = tic;

    Frung = struct();
    for f = 1:numel(FAMILIES3)
        Frung.(local_fieldkey(FAMILIES3{f})) = local_load_rung(FAMILIES3{f}, r, BYRUNG_DIR);
    end

    for p = 1:nPairs
        iA = PAIRS{p, 1}; iB = PAIRS{p, 2};
        famA = FAMILIES3{iA}; famB = FAMILIES3{iB};
        FA = Frung.(local_fieldkey(famA));
        FB = Frung.(local_fieldkey(famB));
        m = local_run_pair(FA, FB, grid3, cfg, VU_mps);

        results(end+1) = struct( ...
            'pairA', famA, 'pairB', famB, ...
            'Ta_nd', Tmax_r, 'Ta_days', Tmax_r * TU_days, ...
            'DVproxy_mps', m.dv_proxy_mps, 'DVlb_mps', m.dvlb_mps, ...
            'DVA_mps', m.dv_A_mps, 'DVB_mps', m.dv_B_mps, ...
            'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
            'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, 'thc_deg', m.thc_deg); %#ok<AGROW>
    end
    clear Frung

    writetable(struct2table(results), resultsCsv);
    fprintf('[stable_stable_ext] Rung %d/%d (Ta=%.2fd) done in %.1fs.\n', r, nRung, Tmax_r * TU_days, toc(tRung));
end

T = struct2table(results);
T = sortrows(T, {'pairA', 'pairB', 'Ta_nd'});
writetable(T, resultsCsv);
save(fullfile(RESULTS_DIR, 'ta_stable_stable_extended_results.mat'), 'T', 'Ta_multiples_of_pi', 'FAMILIES3', 'PAIRS', 'cfg');

fprintf('\n[stable_stable_ext] ══════════ Phase B complete. ══════════\n');
fprintf('  Results table : %s\n', resultsCsv);
disp(T(T.Ta_nd == max(T.Ta_nd), :));

% ══════════════════════════════════════════════════════════════════════════
%  TRAJECTORY RECONSTRUCTION -- at each pair's largest-Ta rung
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[stable_stable_ext] ══════════ TRAJECTORY RECONSTRUCTION ══════════\n');

targetVoxel = containers.Map('KeyType', 'char', 'ValueType', 'double');
pairKey = @(a, b) sprintf('%s|%s', a, b);
for p = 1:nPairs
    iA = PAIRS{p, 1}; iB = PAIRS{p, 2};
    famA = FAMILIES3{iA}; famB = FAMILIES3{iB};
    rowsThisPair = T(strcmp(T.pairA, famA) & strcmp(T.pairB, famB), :);
    [~, iMax] = max(rowsThisPair.Ta_nd);
    targetVoxel(pairKey(famA, famB)) = rowsThisPair.VoxelId(iMax);
end

familyQueries = containers.Map('KeyType', 'char', 'ValueType', 'any');
for f = 1:numel(FAMILIES3)
    familyQueries(FAMILIES3{f}) = {};
end
for p = 1:nPairs
    iA = PAIRS{p, 1}; iB = PAIRS{p, 2};
    famA = FAMILIES3{iA}; famB = FAMILIES3{iB};
    vid = targetVoxel(pairKey(famA, famB));
    q = familyQueries(famA); q{end+1} = struct('partner', famB, 'side', 'FRS', 'voxelId', vid); familyQueries(famA) = q; %#ok<AGROW>
    q = familyQueries(famB); q{end+1} = struct('partner', famA, 'side', 'BRS', 'voxelId', vid); familyQueries(famB) = q; %#ok<AGROW>
end

dvParts = containers.Map('KeyType', 'char', 'ValueType', 'double');
contribCsv = fullfile(RESULTS_DIR, 'ta_stable_stable_extended_contribution_split.csv');

for f = 1:numel(FAMILIES3)
    familyName = FAMILIES3{f};
    queries = familyQueries(familyName);
    if isempty(queries), continue; end

    fprintf('[stable_stable_ext] Loading %s fat atlas (Tmax=64pi, cap=0.15)...\n', familyName);
    tF = tic;
    cfg_full = cfg;
    cfg_full.propag.Tmax = Ta_multiples_of_pi(1) * pi;
    [S_f, cacheInfo] = atlas_prepare_or_load(familyName, cfg_full, grid3);
    fprintf('[stable_stable_ext]   %s in %.1fs.\n', local_hitstr(cacheInfo), toc(tF));

    for qi = 1:numel(queries)
        q = queries{qi};
        [found, dv_mps, ~, iSeed_w, iHead_w, half_w, traj] = ...
            local_find_winner(S_f, grid3, VU_mps, TU_days, q.voxelId, q.side);

        if strcmp(q.side, 'FRS')
            famA = familyName; famB = q.partner;
        else
            famA = q.partner; famB = familyName;
        end

        if found
            fprintf('[stable_stable_ext]   %s (%s side, vs %s): dv=%.4f m/s (seed=%d head=%d half=%s)\n', ...
                familyName, q.side, q.partner, dv_mps, iSeed_w, iHead_w, half_w);
            key = local_fieldkey(familyName);
            sideTag = 'A'; if strcmp(q.side, 'BRS'), sideTag = 'B'; end
            fname = sprintf('%s_vs_%s_%sside.csv', key, matlab.lang.makeValidName(q.partner), sideTag);
            writetable(traj, fullfile(TRAJ_OUT_DIR, fname));
        else
            warning('[stable_stable_ext] %s (%s side, vs %s): NO MATCH FOUND for voxel %d.', ...
                familyName, q.side, q.partner, q.voxelId);
            dv_mps = NaN;
        end

        dvParts(sprintf('%s|%s|%s', famA, famB, q.side)) = dv_mps;

        rows2 = {};
        for p2 = 1:nPairs
            iA2 = PAIRS{p2, 1}; iB2 = PAIRS{p2, 2};
            fA2 = FAMILIES3{iA2}; fB2 = FAMILIES3{iB2};
            keyA = sprintf('%s|%s|FRS', fA2, fB2);
            keyB = sprintf('%s|%s|BRS', fA2, fB2);
            dvA = NaN; dvB = NaN;
            if isKey(dvParts, keyA), dvA = dvParts(keyA); end
            if isKey(dvParts, keyB), dvB = dvParts(keyB); end
            rows2(end+1, :) = {fA2, fB2, targetVoxel(pairKey(fA2, fB2)), dvA, dvB, dvA + dvB}; %#ok<AGROW>
        end
        Tc = cell2table(rows2, 'VariableNames', {'pairA', 'pairB', 'voxelId', 'dv_A_mps', 'dv_B_mps', 'DVlb_check_mps'});
        writetable(Tc, contribCsv);
    end

    clear S_f
end

fprintf('\n[stable_stable_ext] ══════════ DONE ══════════\n');
fprintf('  Contribution split : %s\n', contribCsv);
fprintf('  Trajectory CSVs    : %s\n', TRAJ_OUT_DIR);

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function k = local_fieldkey(familyName)
k = matlab.lang.makeValidName(familyName);
end

function s = local_hitstr(cacheInfo)
if cacheInfo.hit, s = 'loaded from cache'; else, s = 'built'; end
end

% ─────────────────────────────────────────────────────────────────────────────
function F = local_compute_footprint(S, grid3, VU_mps, TU_days)
Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

thm = wrap_to_pi(pi - grid3.th_centers(:));
lut = discretize(thm, grid3.th_edges);
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

pot_u = cr3bp_potential(S.SeedsUpper(:,1), S.SeedsUpper(:,2), S.mu);
v0_upper = sqrt(max(2*pot_u.U(:) - S.CJ, 0));

if isfield(S,'SeedsLower') && ~isempty(S.SeedsLower)
    pot_l = cr3bp_potential(S.SeedsLower(:,1), S.SeedsLower(:,2), S.mu);
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
    ids_brs_u = sub2ind([Ny,Nx,Nt], biy_u(ok_u), ix_u(ok_u), max(1, min(Nt, bit_u(ok_u))));
    dv_bu = dv_u(ok_u);  t_bu = t_u(ok_u);
else
    ids_brs_u = zeros(0,1);  dv_bu = zeros(0,1);  t_bu = zeros(0,1);
end

if ~isempty(ix_l)
    biy_l = Ny - iy_l + 1;
    bit_l = double(it_lut(it_l));
    ok_l  = bit_l > 0;
    ids_brs_l = sub2ind([Ny,Nx,Nt], biy_l(ok_l), ix_l(ok_l), max(1, min(Nt, bit_l(ok_l))));
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

iSeed = double(rows.iSeed(1:n));
iHead = double(rows.iHead(1:n));
t_nd  = double(rows.t(1:n));

lin   = sub2ind([Ns, max_h], iSeed, iHead);
delta = delta_mat(lin);

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
function F = local_load_rung(fam, r, byRungDir)
key = local_fieldkey(fam);
rpath = fullfile(byRungDir, key, sprintf('rung%02d.mat', r));
ck = load(rpath, 'F');
F = ck.F;
end

% ─────────────────────────────────────────────────────────────────────────────
function m = local_run_pair(FA, FB, grid3, cfg, VU_mps)
m = struct('voxelId', NaN, 'dv_proxy_mps', NaN, 'dvlb_mps', NaN, ...
    'dv_A_mps', NaN, 'dv_B_mps', NaN, ...
    'dvpatch_mps', NaN, 'tof_days', NaN, 'xc', NaN, 'yc', NaN, 'thc_deg', NaN);

if isempty(FA) || isempty(FB) || ~isfield(FA, 'uid_frs') || ~isfield(FB, 'uid_brs')
    return;
end

idsO = intersect(FA.uid_frs, FB.uid_brs);
if isempty(idsO), return; end

Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);
[iy, ix, it] = ind2sub([Ny, Nx, Nt], idsO);

bufFrac = 0.05;
if isfield(cfg,'overlap') && isfield(cfg.overlap,'primary_buffer_frac') && ...
        ~isempty(cfg.overlap.primary_buffer_frac)
    bufFrac = cfg.overlap.primary_buffer_frac;
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
ix = ix(ok); iy = iy(ok); it = it(ok);

[~, locA] = ismember(idsO, FA.uid_frs);
[~, locB] = ismember(idsO, FB.uid_brs);
dv_min_A = FA.dv_min_frs(locA);
dv_min_B = FB.dv_min_brs(locB);
t_mean_A = FA.t_mean_frs(locA);
t_mean_B = FB.t_mean_brs(locB);

x_ok = grid3.x_centers(ix);
y_ok = grid3.y_centers(iy);
CJstar = min(FA.CJ, FB.CJ);
pot = cr3bp_potential(x_ok(:), y_ok(:), mu);
v_box = sqrt(max(2 * pot.U - CJstar, 0));
dv_patch_vec = 2 * v_box .* sin(abs(grid3.dtheta) / 2) * VU_mps;
dv_lb_vec    = dv_min_A(:) + dv_min_B(:);
dv_proxy     = dv_lb_vec + dv_patch_vec;

valid = isfinite(dv_proxy);
if ~any(valid), return; end

idxValid  = find(valid);
[~, iLoc] = min(dv_proxy(idxValid));
iWin      = idxValid(iLoc);

m.voxelId      = idsO(iWin);
m.dv_proxy_mps = dv_proxy(iWin);
m.dvlb_mps     = dv_lb_vec(iWin);
m.dv_A_mps     = dv_min_A(iWin);
m.dv_B_mps     = dv_min_B(iWin);
m.dvpatch_mps  = dv_patch_vec(iWin);
m.tof_days     = t_mean_A(iWin) + t_mean_B(iWin);
m.xc           = grid3.x_centers(ix(iWin));
m.yc           = grid3.y_centers(iy(iWin));
m.thc_deg      = rad2deg(grid3.th_centers(it(iWin)));
end

% ─────────────────────────────────────────────────────────────────────────────
function [found, dv_mps, t_days, iSeed_w, iHead_w, half_w, traj] = ...
    local_find_winner(S, grid3, VU_mps, TU_days, targetVoxelId, side)

found = false; dv_mps = NaN; t_days = NaN; iSeed_w = NaN; iHead_w = NaN; half_w = ''; traj = [];

Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

dlists = S.Step4.delta_lists;
Ns     = numel(dlists);
max_h  = max(1, max(cellfun(@numel, dlists)));
delta_mat = zeros(Ns, max_h);
for s = 1:Ns
    v = double(dlists{s});
    delta_mat(s, 1:numel(v)) = v;
end

pot_u = cr3bp_potential(S.SeedsUpper(:,1), S.SeedsUpper(:,2), S.mu);
v0_upper = sqrt(max(2*pot_u.U(:) - S.CJ, 0));
if isfield(S, 'SeedsLower') && ~isempty(S.SeedsLower)
    pot_l = cr3bp_potential(S.SeedsLower(:,1), S.SeedsLower(:,2), S.mu);
    v0_lower = sqrt(max(2*pot_l.U(:) - S.CJ, 0));
else
    v0_lower = v0_upper;
end

if strcmp(side, 'BRS')
    thm = wrap_to_pi(pi - grid3.th_centers(:));
    lut = discretize(thm, grid3.th_edges);
    lut(isnan(lut)) = 0;
    it_lut = uint32(lut);
end

targetVoxelId = uint32(targetVoxelId);
Ny32 = uint32(Ny); Nx32 = uint32(Nx); Nt32 = uint32(Nt); %#ok<NASGU>

bestDv = inf;
halves = {'upper', 'lower'};
for hh = 1:2
    half = halves{hh};
    if strcmp(half, 'upper')
        rows = S.Step4.rows_FRS_upper; v0_per_seed = v0_upper;
    else
        rows = S.Step4.rows_FRS_lower; v0_per_seed = v0_lower;
    end
    n = double(rows.n);
    if n == 0, continue; end

    ix_i = uint32(rows.ix(1:n));
    iy_i = uint32(rows.iy(1:n));
    it_i = uint32(rows.it(1:n));

    if strcmp(side, 'FRS')
        ids = sub2ind([Ny, Nx, Nt], iy_i, ix_i, it_i);
        mask = (ids == targetVoxelId);
    else
        biy_i = Ny32 - iy_i + 1;
        bit_i = it_lut(it_i);
        okMirror = bit_i > 0;
        ids = zeros(n, 1, 'uint32');
        ids(okMirror) = sub2ind([Ny, Nx, Nt], biy_i(okMirror), ix_i(okMirror), ...
            max(uint32(1), min(Nt32, bit_i(okMirror))));
        mask = okMirror & (ids == targetVoxelId);
        clear biy_i bit_i okMirror
    end
    clear ix_i iy_i it_i ids

    if ~any(mask), continue; end
    idxCand = find(mask);
    clear mask
    iSeedCand = double(rows.iSeed(idxCand));
    iHeadCand = double(rows.iHead(idxCand));
    tCand     = double(rows.t(idxCand));

    lin   = sub2ind([Ns, max_h], iSeedCand, iHeadCand);
    delta = delta_mat(lin);
    v0    = v0_per_seed(iSeedCand);
    dv    = 2 * v0(:) .* sin(abs(delta(:)) / 2) * VU_mps;

    [thisMin, iLoc] = min(dv);
    if thisMin < bestDv
        bestDv     = thisMin;
        found      = true;
        dv_mps     = thisMin;
        t_days     = abs(tCand(iLoc)) * TU_days;
        t_match_nd = tCand(iLoc);
        iSeed_w    = iSeedCand(iLoc);
        iHead_w    = iHeadCand(iLoc);
        half_w     = half;
    end
end

if found
    if strcmp(half_w, 'upper')
        rows = S.Step4.rows_FRS_upper;
    else
        rows = S.Step4.rows_FRS_lower;
    end
    n = double(rows.n);
    trajMask = (double(rows.iSeed(1:n)) == iSeed_w) & (double(rows.iHead(1:n)) == iHead_w);
    idxT = find(trajMask);
    tt = double(rows.t(idxT));
    [tt, order] = sort(tt);
    idxT = idxT(order);

    [~, matchIdxInSorted] = min(abs(tt - t_match_nd));
    tt   = tt(1:matchIdxInSorted);
    idxT = idxT(1:matchIdxInSorted);

    ixT = double(rows.ix(idxT));
    iyT = double(rows.iy(idxT));
    itT = double(rows.it(idxT));

    xVals  = reshape(grid3.x_centers(ixT), [], 1);
    yVals  = reshape(grid3.y_centers(iyT), [], 1);
    thVals = reshape(rad2deg(grid3.th_centers(itT)), [], 1);
    traj = table(tt(:), tt(:) * TU_days, xVals, yVals, thVals, ...
        'VariableNames', {'t_nd', 't_days', 'x_nd', 'y_nd', 'theta_deg'});
end
end
