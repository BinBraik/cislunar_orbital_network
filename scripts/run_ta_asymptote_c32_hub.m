%% RUN_TA_ASYMPTOTE_C32_HUB
% Shane's economical test of the zero-cost unstable-unstable subgraph
% hypothesis: rather than computing all 45 unstable-unstable pairs, use
% (3,2)-cycler (C32) -- the paper's established hub/gateway/relay center
% -- as a single reference family. Compute the direct ΔV_lb(Ta) curve
% from C32 to each of the other 9 unstable families and each of the 3
% stable families (12 pairs total).
%
% If all 9 unstable-partner curves -> 0 as Ta -> infinity, then for any
% two unstable families F_i, F_j, the relay path F_i -> C32 -> F_j has
% asymptotic cost 0 + 0 = 0, establishing that all 10 unstable families
% form one zero-cost connected component WITHOUT computing all 45 direct
% pairs. The 3 stable-partner curves are predicted to approach nonzero
% nonzero lower bounds (the paper's DV_min,s) instead.
%
% This is a companion to run_ta_asymptote_sweep_full.m (the full 78-pair
% sweep), NOT a replacement -- both are intended to run at the same time.
% Same physics/config: DV_cap_nd = 0.2 fixed, Table 3 grid unchanged,
% same 12-rung Ta ladder, Ta=16pi fat build derived downward.
%
% CRITICAL: this script points FOOTPRINT_DIR at the SAME directory the
% full-sweep script writes to (ta_asymptote_full_results/footprints/),
% so any family the other job has already finished is reused instantly
% here (and vice versa -- whatever this script builds is picked up by
% the other job automatically too, since both use the identical cache
% file naming/schema and Ta ladder / DV_cap fingerprint check).
%
% atlas_cache_save.m does a plain (non-atomic) save(), so two jobs must
% never build the SAME family at the same time -- that's a real risk of
% a corrupted cache file, not just wasted compute. SKIP_FAMILIES below is
% the manual guard for that: list any family currently being built by
% the other job's log so this script will not touch it. Remove an entry
% once the other job's log confirms that family is done (its footprint
% cache file will exist and get reused automatically -- no code change
% needed beyond editing this list, just rerun the script).

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

HUB_FAMILY = 'Cycler 32';

PARTNER_FAMILIES = { ...
    'Lyapunov L1', 'Lyapunov L2', 'Cycler 21', 'Cycler 11a', 'Cycler 11b', ...
    'Resonant 2to1 Unstable', 'Resonant 3to1 Unstable', 'Resonant 5to2 Unstable', ...
    'Distant Prograde Orbit', ...
    'Resonant 2to1 Stable', 'Resonant 3to1 Stable', 'Resonant 5to2 Stable' ...
};

% Families currently being built by the OTHER (full 78-pair) job -- do
% NOT attempt to build these here, to avoid a concurrent-write race on
% the same atlas_cache/ file. Per the running job's log at the time this
% script was written, family 11/13 (Resonant 5to2 Stable) was mid-build
% (Step 4 done, still deriving per-rung footprints). Edit/clear this list
% once that job's log shows it finished -- the footprint cache file will
% then simply be found and reused, no other change needed.
SKIP_FAMILIES = {'Resonant 5to2 Stable'};

% Same ladder as run_ta_asymptote_sweep_full.m -- MUST match bit-for-bit
% for the shared footprint cache's Ta_multiples_of_pi equality check to
% pass (both compute it with the identical deterministic expression, so
% this holds automatically as long as neither script edits it).
Ta_multiples_of_pi = sort(2.^((-3:8)/2), 'descend');

N_WORKERS = 60;

% Points at the SAME footprint cache the full-sweep job uses -- this is
% the whole point (mutual reuse, not a separate cache).
FULL_SWEEP_OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_full_results');
FOOTPRINT_DIR = fullfile(FULL_SWEEP_OUTPUT_DIR, 'footprints');

OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_c32_hub_results');

% ══════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION  (Table 3 / paper reference values, unchanged)
% ══════════════════════════════════════════════════════════════════════════
cfg = atlas_cfg_defaults();

cfg.grid.dx               = 0.001;
cfg.grid.dy               = 0.001;
cfg.grid.dtheta           = deg2rad(1);
cfg.seed.ds_seed          = 0.01;
cfg.fan.DV_cap_nd         = 0.2;
cfg.fan.dtheta_fan        = deg2rad(0.5);
cfg.propag.absTol         = 1e-8;
cfg.propag.relTol         = 1e-8;
cfg.propag.v2tol          = 1e-8;
cfg.log.step_len_factor   = 0.75;
cfg.log.maxstep_factor    = 2;

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

if ~exist(OUTPUT_DIR, 'dir'),    mkdir(OUTPUT_DIR);    end
if ~exist(FOOTPRINT_DIR, 'dir'), mkdir(FOOTPRINT_DIR); end

nRung = numel(Ta_multiples_of_pi);
grid3 = atlas_grid_make(cfg);
Tmax_fat = Ta_multiples_of_pi(1) * pi;

fprintf('[c32_hub] Hub family: %s\n', HUB_FAMILY);
fprintf('[c32_hub] Partner families: %d   Skipped this run: %s\n', ...
    numel(PARTNER_FAMILIES), strjoin(SKIP_FAMILIES, ', '));
fprintf('[c32_hub] Shared footprint cache: %s\n', FOOTPRINT_DIR);

% ══════════════════════════════════════════════════════════════════════════
%  PARALLEL POOL
% ══════════════════════════════════════════════════════════════════════════
if N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[c32_hub] Started parpool with %d workers.\n', N_WORKERS);
    end
end

% ══════════════════════════════════════════════════════════════════════════
%  PER-FAMILY FOOTPRINTS  (shared cache: reuse if present, build if not)
% ══════════════════════════════════════════════════════════════════════════
ALL_NEEDED = [{HUB_FAMILY}, PARTNER_FAMILIES];
Fp = containers.Map();   % family key -> 1xnRung cell of footprints

for i = 1:numel(ALL_NEEDED)
    fam = ALL_NEEDED{i};
    key = local_fieldkey(fam);
    fpath = fullfile(FOOTPRINT_DIR, [key '.mat']);

    if isfile(fpath)
        ck = load(fpath, 'Fcell', 'Ta_multiples_of_pi', 'DV_cap_nd');
        if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) && ...
                ck.DV_cap_nd == cfg.fan.DV_cap_nd
            Fp(fam) = ck.Fcell;
            fprintf('[c32_hub] %-25s -- footprints found in shared cache, reused instantly.\n', fam);
            continue;
        end
    end

    if any(strcmp(SKIP_FAMILIES, fam))
        fprintf('[c32_hub] %-25s -- SKIPPED (in SKIP_FAMILIES; other job may be building it now).\n', fam);
        continue;
    end

    tFam = tic;
    fprintf('[c32_hub] %-25s -- not cached, building fat atlas at Ta=%.4gpi (%.2f days)...\n', ...
        fam, Ta_multiples_of_pi(1), Tmax_fat * TU_days);

    cfg_fat = cfg;
    cfg_fat.propag.Tmax = Tmax_fat;
    [S_fat_i, cacheInfo] = atlas_prepare_or_load(fam, cfg_fat, grid3);
    fprintf('[c32_hub]   fat atlas %s in %.1fs.\n', local_hitstr(cacheInfo), toc(tFam));

    Fcell = cell(1, nRung);
    for r = 1:nRung
        Tmax_r = Ta_multiples_of_pi(r) * pi;
        if r == 1
            Ssub = S_fat_i;
        else
            cfg_r = cfg;
            cfg_r.propag.Tmax = Tmax_r;
            Ssub = atlas_derive_subset(S_fat_i, cfg_r);
        end
        Fcell{r} = local_compute_footprint(Ssub, grid3, VU_mps, TU_days);
        clear Ssub
    end
    clear S_fat_i

    DV_cap_nd = cfg.fan.DV_cap_nd; %#ok<NASGU>
    save(fpath, 'Fcell', 'Ta_multiples_of_pi', 'DV_cap_nd', '-v7.3');
    Fp(fam) = Fcell;

    fprintf('[c32_hub] %-25s -- done in %.1fs total (now cached for both jobs).\n', fam, toc(tFam));
end

% ══════════════════════════════════════════════════════════════════════════
%  PAIRWISE OVERLAP: hub vs. every AVAILABLE partner, all rungs
% ══════════════════════════════════════════════════════════════════════════
availablePartners = PARTNER_FAMILIES(cellfun(@(f) isKey(Fp, f), PARTNER_FAMILIES));
skippedPartners   = setdiff(PARTNER_FAMILIES, availablePartners);

fprintf('\n[c32_hub] Computing %d/%d pairs (missing: %s)\n', ...
    numel(availablePartners), numel(PARTNER_FAMILIES), strjoin(skippedPartners, ', '));

results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});

FHub = Fp(HUB_FAMILY);

for pIdx = 1:numel(availablePartners)
    fam = availablePartners{pIdx};
    Fpartner = Fp(fam);

    for r = 1:nRung
        Tmax_r = Ta_multiples_of_pi(r) * pi;
        m = local_run_pair(FHub{r}, Fpartner{r}, grid3, cfg, VU_mps);

        row = struct( ...
            'pairA', HUB_FAMILY, 'pairB', fam, ...
            'Ta_nd', Tmax_r, 'Ta_days', Tmax_r * TU_days, ...
            'DVproxy_mps', m.dv_proxy_mps, 'DVlb_mps', m.dvlb_mps, ...
            'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
            'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, 'thc_deg', m.thc_deg);
        results(end+1) = row; %#ok<AGROW>
    end
    fprintf('[c32_hub]   %s -> %s : %d rungs done.\n', HUB_FAMILY, fam, nRung);
end

% ══════════════════════════════════════════════════════════════════════════
%  SAVE RESULTS
% ══════════════════════════════════════════════════════════════════════════
T = struct2table(results);
T = sortrows(T, {'pairB', 'Ta_nd'});

outMat = fullfile(OUTPUT_DIR, 'ta_asymptote_c32_hub_results.mat');
outCsv = fullfile(OUTPUT_DIR, 'ta_asymptote_c32_hub_results.csv');
save(outMat, 'T', 'Ta_multiples_of_pi', 'HUB_FAMILY', 'PARTNER_FAMILIES', 'skippedPartners', 'cfg');
writetable(T, outCsv);

fprintf('\n[c32_hub] ══════════ DONE ══════════\n');
fprintf('  Computed pairs : %d/%d\n', numel(availablePartners), numel(PARTNER_FAMILIES));
if ~isempty(skippedPartners)
    fprintf('  Still missing  : %s\n', strjoin(skippedPartners, ', '));
    fprintf('  -> rerun this script once those families finish (in this job''s SKIP_FAMILIES\n');
    fprintf('     or the other job''s log) -- already-done pairs are skipped/reused, only the\n');
    fprintf('     missing ones get computed.\n');
end
fprintf('  Results table  : %s\n', outMat);
fprintf('  CSV            : %s\n', outCsv);

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  (identical to run_ta_asymptote_sweep_full.m)
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
function m = local_run_pair(FA, FB, grid3, cfg, VU_mps)
m = struct('voxelId', NaN, 'dv_proxy_mps', NaN, 'dvlb_mps', NaN, ...
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
m.dvpatch_mps  = dv_patch_vec(iWin);
m.tof_days     = t_mean_A(iWin) + t_mean_B(iWin);
m.xc           = grid3.x_centers(ix(iWin));
m.yc           = grid3.y_centers(iy(iWin));
m.thc_deg      = rad2deg(grid3.th_centers(it(iWin)));
end
