%% RUN_TA_ASYMPTOTE_SWEEP_FULL
% Full 13-family / 78-pair extension of run_ta_asymptote_sweep.m, testing
% the same chaotic-sea asymptotic hypothesis across the whole family set
% instead of just DPO<->R31-U.
%
% Same physics/config as the single-pair run: DV_cap_nd = 0.2 (paper
% nominal) held FIXED, Table 3 grid/discretization unchanged, Ta ladder
% built largest-first with every smaller rung derived from that one fat
% atlas via atlas_derive_subset (shrink-only, never expansion).
%
% Two changes from the single-pair script:
%
%   1. PAIRS is now all 78 unordered pairs among the 13 families,
%      auto-generated -- not hand-listed.
%
%   2. Memory-bounded two-phase architecture. Holding all 13 families'
%      fat atlases (Ta=16pi) in memory at once -- what the single-pair
%      script does with its 2 families -- does not scale to 13. Instead:
%
%        PHASE A (family-by-family): build the fat atlas for family i,
%          immediately derive every rung's compact footprint (~5-25MB,
%          vs. ~0.5-2GB+ for a full atlas) from it, cache those footprints
%          to disk, then DISCARD the full atlas before moving to family
%          i+1. Peak memory is therefore ~one family's fat atlas plus the
%          (small) footprint cache -- not 13 fat atlases simultaneously.
%          Checkpointed per family: a family whose footprint cache file
%          already exists (and matches the current Ta ladder / DV_cap)
%          is skipped entirely on resubmission, no re-integration.
%
%        PHASE B (rung-by-rung, pair-by-pair): using ONLY the cached
%          footprints from Phase A, compute the 78-pair overlap/proxy
%          decomposition per rung. Pair loop is parfor, matching
%          run_overlap_dv_tmax_sweep.m's pattern -- footprints for the
%          pairs in play are pre-sliced into flat per-pair arrays BEFORE
%          the parfor so each worker only receives the ~2 footprints it
%          needs, not the whole (family x rung) cache broadcast to every
%          worker. Checkpointed per rung.
%
%   Ta ladder reverted to the original 12-rung ladder (2.^((-3:8)/2)*pi,
%   i.e. 0.354pi ... 16pi) -- the denser 25-rung ladder was for the single
%   DPO->R31-U deep-dive and isn't needed for a first full-network pass.

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

FAMILIES = { ...
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

% Ta ladder, expressed as multiples of pi (nd). Back to the original
% 12-rung geometric (sqrt(2)-ratio) ladder: 0.354pi ... 16pi.
Ta_multiples_of_pi = sort(2.^((-3:8)/2), 'descend');

% Parallel workers -- reused for both the (seed x heading) propagation in
% Phase A and the pair loop in Phase B. Match --cpus-per-task in SLURM
% (N_WORKERS + 1 for the MATLAB driver).
N_WORKERS = 60;

OUTPUT_DIR    = fullfile(repoRoot, 'ta_asymptote_full_results');
FOOTPRINT_DIR = fullfile(OUTPUT_DIR, 'footprints');
CHECKPOINT_B  = fullfile(OUTPUT_DIR, 'checkpoint_phaseB.mat');

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

nFam  = numel(FAMILIES);
nRung = numel(Ta_multiples_of_pi);

% ── All 78 unordered pairs, auto-generated ─────────────────────────────────
nPairs = nFam * (nFam - 1) / 2;
pairIdxA = zeros(nPairs, 1);
pairIdxB = zeros(nPairs, 1);
p = 0;
for i = 1:nFam-1
    for j = i+1:nFam
        p = p + 1;
        pairIdxA(p) = i;
        pairIdxB(p) = j;
    end
end

fprintf('[full_sweep] Families: %d   Pairs: %d   Ta rungs: %d\n', nFam, nPairs, nRung);
fprintf('[full_sweep] Ta ladder (xpi): %s\n', mat2str(Ta_multiples_of_pi, 4));
fprintf('[full_sweep] DV_cap_nd = %.4g (%.1f m/s one-sided)\n', ...
    cfg.fan.DV_cap_nd, cfg.fan.DV_cap_nd * VU_mps);

% ══════════════════════════════════════════════════════════════════════════
%  PARALLEL POOL
% ══════════════════════════════════════════════════════════════════════════
if N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[full_sweep] Started parpool with %d workers.\n', N_WORKERS);
    end
end

grid3 = atlas_grid_make(cfg);
Tmax_fat = Ta_multiples_of_pi(1) * pi;

% ══════════════════════════════════════════════════════════════════════════
%  PHASE A -- family-by-family: fat atlas -> per-rung footprints -> discard
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[full_sweep] ══════════ PHASE A: per-family footprints ══════════\n');

Fp = cell(nFam, nRung);   % in-memory footprint cache, reused in Phase B

for i = 1:nFam
    key   = local_fieldkey(FAMILIES{i});
    fpath = fullfile(FOOTPRINT_DIR, [key '.mat']);

    if isfile(fpath)
        ck = load(fpath, 'Fcell', 'Ta_multiples_of_pi', 'DV_cap_nd');
        if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) && ...
                ck.DV_cap_nd == cfg.fan.DV_cap_nd
            Fp(i, :) = ck.Fcell;
            fprintf('[full_sweep] (%d/%d) %-25s -- footprints cached, skipping.\n', ...
                i, nFam, FAMILIES{i});
            continue;
        else
            fprintf('[full_sweep] (%d/%d) %-25s -- cached footprints stale (config changed), rebuilding.\n', ...
                i, nFam, FAMILIES{i});
        end
    end

    tFam = tic;
    fprintf('[full_sweep] (%d/%d) %-25s -- building fat atlas at Ta=%.4gpi (%.2f days)...\n', ...
        i, nFam, FAMILIES{i}, Ta_multiples_of_pi(1), Tmax_fat * TU_days);

    cfg_fat = cfg;
    cfg_fat.propag.Tmax = Tmax_fat;
    [S_fat_i, cacheInfo] = atlas_prepare_or_load(FAMILIES{i}, cfg_fat, grid3);
    fprintf('[full_sweep]   fat atlas %s in %.1fs.\n', local_hitstr(cacheInfo), toc(tFam));

    Fcell = cell(1, nRung);
    for r = 1:nRung
        Tmax_r = Ta_multiples_of_pi(r) * pi;
        if r == 1
            Ssub = S_fat_i;   % top rung IS the fat atlas, no filtering needed
        else
            cfg_r = cfg;
            cfg_r.propag.Tmax = Tmax_r;
            Ssub = atlas_derive_subset(S_fat_i, cfg_r);
        end
        Fcell{r} = local_compute_footprint(Ssub, grid3, VU_mps, TU_days);
        clear Ssub
    end
    clear S_fat_i   % discard the full atlas -- only the footprints survive

    DV_cap_nd = cfg.fan.DV_cap_nd; %#ok<NASGU> -- saved for cache validation
    save(fpath, 'Fcell', 'Ta_multiples_of_pi', 'DV_cap_nd', '-v7.3');
    Fp(i, :) = Fcell;

    fprintf('[full_sweep] (%d/%d) %-25s -- done in %.1fs total.\n', ...
        i, nFam, FAMILIES{i}, toc(tFam));
end

fprintf('[full_sweep] Phase A complete. All %d families'' footprints cached.\n', nFam);

% ══════════════════════════════════════════════════════════════════════════
%  PHASE B -- rung-by-rung, all 78 pairs (parfor), footprints pre-sliced
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[full_sweep] ══════════ PHASE B: pairwise overlap, all %d pairs ══════════\n', nPairs);

results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});
rung_done = false(nRung, 1);

if isfile(CHECKPOINT_B)
    fprintf('[full_sweep] Phase B checkpoint found -- loading...\n');
    ck = load(CHECKPOINT_B, 'results', 'rung_done', 'Ta_multiples_of_pi');
    if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi)
        results   = ck.results;
        rung_done = ck.rung_done;
        fprintf('[full_sweep] %d/%d rungs already done.\n', sum(rung_done), nRung);
    else
        warning('[full_sweep] Phase B checkpoint Ta ladder mismatch -- starting fresh.');
    end
end

for r = 1:nRung

    if rung_done(r)
        fprintf('[full_sweep] Rung %d/%d already done -- skipping.\n', r, nRung);
        continue;
    end

    Tmax_r = Ta_multiples_of_pi(r) * pi;
    fprintf('\n[full_sweep] ── Rung %d/%d: Ta = %.4gpi (%.2f days) ──\n', ...
        r, nRung, Ta_multiples_of_pi(r), Tmax_r * TU_days);
    tRung = tic;

    % Pre-slice footprints for THIS rung into flat per-pair arrays so the
    % parfor below receives only the ~2 footprints each iteration needs,
    % not the whole (family x rung) cache broadcast to every worker.
    FA_arr = Fp(pairIdxA, r);
    FB_arr = Fp(pairIdxB, r);

    % NOTE: FAMILIES is a 1xN row cell, and MATLAB indexing a vector with
    % any index array preserves the ORIENTATION OF THE ARRAY BEING INDEXED
    % (not the index's shape) -- so this is already 1xnPairs. Do not
    % transpose: a column shape here would make `results` grow into a 2-D
    % matrix of structs across rungs instead of a flat, growing list, and
    % struct2table would only fail on that at the very end of the run.
    pair_famA = FAMILIES(pairIdxA);
    pair_famB = FAMILIES(pairIdxB);

    pair_DVproxy = nan(1, nPairs);
    pair_DVlb    = nan(1, nPairs);
    pair_DVpatch = nan(1, nPairs);
    pair_TOF     = nan(1, nPairs);
    pair_Voxel   = nan(1, nPairs);
    pair_xc      = nan(1, nPairs);
    pair_yc      = nan(1, nPairs);
    pair_thc     = nan(1, nPairs);

    cfg_b = cfg;   % broadcast copy for the parfor (sys/overlap fields only)

    if N_WORKERS > 0
        parfor p = 1:nPairs
            m = local_run_pair(FA_arr{p}, FB_arr{p}, grid3, cfg_b, VU_mps); %#ok<PFBNS>
            pair_DVproxy(p) = m.dv_proxy_mps;
            pair_DVlb(p)    = m.dvlb_mps;
            pair_DVpatch(p) = m.dvpatch_mps;
            pair_TOF(p)     = m.tof_days;
            pair_Voxel(p)   = m.voxelId;
            pair_xc(p)      = m.xc;
            pair_yc(p)      = m.yc;
            pair_thc(p)     = m.thc_deg;
        end
    else
        for p = 1:nPairs
            m = local_run_pair(FA_arr{p}, FB_arr{p}, grid3, cfg_b, VU_mps);
            pair_DVproxy(p) = m.dv_proxy_mps;
            pair_DVlb(p)    = m.dvlb_mps;
            pair_DVpatch(p) = m.dvpatch_mps;
            pair_TOF(p)     = m.tof_days;
            pair_Voxel(p)   = m.voxelId;
            pair_xc(p)      = m.xc;
            pair_yc(p)      = m.yc;
            pair_thc(p)     = m.thc_deg;
        end
    end

    rowsThisRung = struct( ...
        'pairA', pair_famA, 'pairB', pair_famB, ...
        'Ta_nd', num2cell(repmat(Tmax_r, 1, nPairs)), ...
        'Ta_days', num2cell(repmat(Tmax_r * TU_days, 1, nPairs)), ...
        'DVproxy_mps', num2cell(pair_DVproxy), 'DVlb_mps', num2cell(pair_DVlb), ...
        'DVpatch_mps', num2cell(pair_DVpatch), 'TOF_days', num2cell(pair_TOF), ...
        'VoxelId', num2cell(pair_Voxel), 'xc_nd', num2cell(pair_xc), ...
        'yc_nd', num2cell(pair_yc), 'thc_deg', num2cell(pair_thc));

    results = [results, rowsThisRung]; %#ok<AGROW>
    rung_done(r) = true;

    save(CHECKPOINT_B, 'results', 'rung_done', 'Ta_multiples_of_pi', '-v7.3');
    fprintf('[full_sweep] Rung %d/%d done in %.1fs -- checkpoint saved.\n', r, nRung, toc(tRung));

end

clear Fp   % done with the footprint cache

% ══════════════════════════════════════════════════════════════════════════
%  SAVE FINAL RESULTS
% ══════════════════════════════════════════════════════════════════════════
T = struct2table(results);
T = sortrows(T, {'pairA', 'pairB', 'Ta_nd'});

outMat = fullfile(OUTPUT_DIR, 'ta_asymptote_full_results.mat');
outCsv = fullfile(OUTPUT_DIR, 'ta_asymptote_full_results.csv');
save(outMat, 'T', 'Ta_multiples_of_pi', 'FAMILIES', 'cfg');
writetable(T, outCsv);

fprintf('\n[full_sweep] ══════════ DONE ══════════\n');
fprintf('  Results table : %s\n', outMat);
fprintf('  CSV           : %s\n', outCsv);

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
%LOCAL_COMPUTE_FOOTPRINT  Compact per-voxel (dv_min, t_mean) summary for FRS/BRS.
% Identical formulation to run_ta_asymptote_sweep.m / run_overlap_dv_tmax_sweep.m.

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
%LOCAL_RUN_PAIR  Overlap + proxy decomposition for one pair at one rung.
% Same eq. (48)-matching convention as run_ta_asymptote_sweep.m: the
% winning voxel minimizes the TOTAL proxy cost; DVlb/DVpatch reported here
% are the decomposition of that SAME voxel.

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
