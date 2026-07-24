%% RUN_TA_ASYMPTOTE_C32_HUB_FAST
% Faster, still-safe version of run_ta_asymptote_c32_hub.m, scoped to
% ONLY the families already fully cached (per your request: focus on C32
% pairs with the finished families now; the big 78-pair job can keep
% running separately, or not, at your leisure).
%
% WHY THE ORIGINAL HUB SCRIPT / PHASE-B-SAFE SCRIPT ARE SLOW HERE:
% the footprint cache files run 4-40GB EACH (measured on disk), bundled
% as all 12 rungs in one file per family. The memory-safe rung-at-a-time
% design re-reads each family's FULL bundle file on every one of the 12
% rungs just to extract 1/12 of it and discard the rest -- for the 10
% cached families that's up to ~2.4TB of redundant disk I/O total. That
% is the real cost, not the (fast, vectorized) intersect/argmin math.
%
% FIX: two steps.
%   STEP 1 (one-time, idempotent): read each family's big bundle file
%     ONCE, split it into 12 small per-rung files (~1-3GB each instead of
%     the full 4-40GB), save, discard. Parallelized across families with
%     a DELIBERATELY BOUNDED pool (N_WORKERS below) -- this is bounded by
%     worst-case simultaneous file size, not CPU count. Do not set
%     N_WORKERS to 60 here: that would mean up to 60 concurrent multi-GB
%     reads (~2.4TB peak) and defeat the entire point. Already-repacked
%     families are skipped instantly on rerun.
%   STEP 2: the actual C32-vs-partner computation now loads only the
%     small per-rung files -- fast, and safe to parallelize broadly since
%     each load is a few GB, not tens of GB.
%
% Only operates on families whose ORIGINAL bundle already exists in the
% shared footprints/ cache (i.e. already fully built by the other jobs).
% Families not yet built (per your last check: Resonant 5to2 Stable,
% Resonant 5to2 Unstable, Distant Prograde Orbit) are simply reported as
% missing -- rerun this script later once they're built and it picks
% them up automatically, same idempotent design as the other scripts.

clear; clc;

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

Ta_multiples_of_pi = sort(2.^((-3:8)/2), 'descend');   % must match the writing scripts
DV_cap_nd_expected  = 0.2;

% Bounded on purpose -- see header. This caps worst-case simultaneous
% large-file reads during the one-time repack step (Step 1). 8 x ~40GB
% (the biggest family seen so far, Cycler 32 itself) = ~320GB peak,
% comfortably inside a 400G allocation. Step 2 (small per-rung files) is
% safe at this same worker count too -- no need for a second pool size.
N_WORKERS = 8;

FULL_SWEEP_OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_full_results');
BUNDLE_DIR = fullfile(FULL_SWEEP_OUTPUT_DIR, 'footprints');            % existing big per-family files
BYRUNG_DIR = fullfile(FULL_SWEEP_OUTPUT_DIR, 'footprints_by_rung');    % new small per-family-per-rung files

OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_c32_hub_results');
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end
if ~exist(BYRUNG_DIR, 'dir'), mkdir(BYRUNG_DIR); end

cfg = atlas_cfg_defaults();
cfg.grid.dx      = 0.001;
cfg.grid.dy      = 0.001;
cfg.grid.dtheta  = deg2rad(1);
cfg.fan.DV_cap_nd = DV_cap_nd_expected;

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
grid3   = atlas_grid_make(cfg);
nRung   = numel(Ta_multiples_of_pi);

ALL_NEEDED = [{HUB_FAMILY}, PARTNER_FAMILIES];

% ── Which families are actually available (bundle file exists)? ────────────
famReady = false(size(ALL_NEEDED));
for i = 1:numel(ALL_NEEDED)
    famReady(i) = isfile(fullfile(BUNDLE_DIR, [local_fieldkey(ALL_NEEDED{i}) '.mat']));
end
readyFamilies = ALL_NEEDED(famReady);
missingFamilies = ALL_NEEDED(~famReady);

fprintf('[c32_fast] %d/%d families ready (bundle exists): %s\n', ...
    numel(readyFamilies), numel(ALL_NEEDED), strjoin(readyFamilies, ', '));
if ~isempty(missingFamilies)
    fprintf('[c32_fast] Not yet built, skipped this run: %s\n', strjoin(missingFamilies, ', '));
end

if ~any(strcmp(readyFamilies, HUB_FAMILY))
    error('[c32_fast] Hub family %s is not built yet -- nothing to do.', HUB_FAMILY);
end

% ══════════════════════════════════════════════════════════════════════════
%  PARALLEL POOL  (deliberately small -- see header)
% ══════════════════════════════════════════════════════════════════════════
if N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[c32_fast] Started parpool with %d workers (bounded by design).\n', N_WORKERS);
    end
end

% ══════════════════════════════════════════════════════════════════════════
%  STEP 1 -- one-time repack: big bundle -> 12 small per-rung files
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[c32_fast] ══════════ STEP 1: repack bundles into per-rung files ══════════\n');

needsRepack = false(numel(readyFamilies), 1);
for i = 1:numel(readyFamilies)
    famDir = fullfile(BYRUNG_DIR, local_fieldkey(readyFamilies{i}));
    if ~exist(famDir, 'dir'), needsRepack(i) = true; continue; end
    for r = 1:nRung
        if ~isfile(fullfile(famDir, sprintf('rung%02d.mat', r)))
            needsRepack(i) = true;
            break;
        end
    end
end

toRepack = readyFamilies(needsRepack);
fprintf('[c32_fast] %d/%d families already repacked; %d need repacking: %s\n', ...
    sum(~needsRepack), numel(readyFamilies), numel(toRepack), strjoin(toRepack, ', '));

if N_WORKERS > 0 && numel(toRepack) > 1
    parfor i = 1:numel(toRepack)
        local_repack_family(toRepack{i}, BUNDLE_DIR, BYRUNG_DIR, Ta_multiples_of_pi, cfg.fan.DV_cap_nd); %#ok<PFBNS>
    end
else
    for i = 1:numel(toRepack)
        local_repack_family(toRepack{i}, BUNDLE_DIR, BYRUNG_DIR, Ta_multiples_of_pi, cfg.fan.DV_cap_nd);
    end
end
fprintf('[c32_fast] Step 1 complete.\n');

% ══════════════════════════════════════════════════════════════════════════
%  STEP 2 -- fast pairwise computation, rung by rung, small files only
% ══════════════════════════════════════════════════════════════════════════
availablePartners = intersect(PARTNER_FAMILIES, readyFamilies, 'stable');
fprintf('\n[c32_fast] ══════════ STEP 2: %d pairs x %d rungs ══════════\n', ...
    numel(availablePartners), nRung);

results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});

for r = 1:nRung
    Tmax_r = Ta_multiples_of_pi(r) * pi;
    tRung = tic;

    FHub_r = local_load_rung(HUB_FAMILY, r, BYRUNG_DIR);

    nP = numel(availablePartners);
    rowsThisRung = cell(1, nP);
    if N_WORKERS > 0 && nP > 1
        parfor pIdx = 1:nP
            Fp_r = local_load_rung(availablePartners{pIdx}, r, BYRUNG_DIR); %#ok<PFBNS>
            m = local_run_pair(FHub_r, Fp_r, grid3, cfg, VU_mps); %#ok<PFBNS>
            rowsThisRung{pIdx} = struct( ...
                'pairA', HUB_FAMILY, 'pairB', availablePartners{pIdx}, ...
                'Ta_nd', Tmax_r, 'Ta_days', Tmax_r * TU_days, ...
                'DVproxy_mps', m.dv_proxy_mps, 'DVlb_mps', m.dvlb_mps, ...
                'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
                'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, 'thc_deg', m.thc_deg);
        end
    else
        for pIdx = 1:nP
            Fp_r = local_load_rung(availablePartners{pIdx}, r, BYRUNG_DIR);
            m = local_run_pair(FHub_r, Fp_r, grid3, cfg, VU_mps);
            rowsThisRung{pIdx} = struct( ...
                'pairA', HUB_FAMILY, 'pairB', availablePartners{pIdx}, ...
                'Ta_nd', Tmax_r, 'Ta_days', Tmax_r * TU_days, ...
                'DVproxy_mps', m.dv_proxy_mps, 'DVlb_mps', m.dvlb_mps, ...
                'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
                'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, 'thc_deg', m.thc_deg);
        end
    end

    for pIdx = 1:nP
        results(end+1) = rowsThisRung{pIdx}; %#ok<AGROW>
    end
    clear FHub_r rowsThisRung

    fprintf('[c32_fast] Rung %d/%d (Ta=%.2fd) done in %.1fs.\n', r, nRung, Tmax_r * TU_days, toc(tRung));
end

T = struct2table(results);
T = sortrows(T, {'pairB', 'Ta_nd'});

outMat = fullfile(OUTPUT_DIR, 'ta_asymptote_c32_hub_results.mat');
outCsv = fullfile(OUTPUT_DIR, 'ta_asymptote_c32_hub_results.csv');
save(outMat, 'T', 'Ta_multiples_of_pi', 'HUB_FAMILY', 'availablePartners', 'missingFamilies', 'cfg');
writetable(T, outCsv);

fprintf('\n[c32_fast] ══════════ DONE ══════════\n');
fprintf('  Computed pairs : %d/%d\n', numel(availablePartners), numel(PARTNER_FAMILIES));
if ~isempty(missingFamilies)
    fprintf('  Still missing  : %s (rerun once built)\n', strjoin(setdiff(missingFamilies, {HUB_FAMILY}), ', '));
end
fprintf('  Results table  : %s\n', outMat);
fprintf('  CSV            : %s\n', outCsv);

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function k = local_fieldkey(familyName)
k = matlab.lang.makeValidName(familyName);
end

% ─────────────────────────────────────────────────────────────────────────────
function local_repack_family(fam, bundleDir, byRungDir, Ta_multiples_of_pi, DV_cap_expected)
%LOCAL_REPACK_FAMILY  Read one family's big bundle file ONCE, write 12
% small per-rung files, discard. Safe to call from a bounded parfor: each
% call touches only its own family's files, no shared state.
key = local_fieldkey(fam);
bundlePath = fullfile(bundleDir, [key '.mat']);
famDir = fullfile(byRungDir, key);
if ~exist(famDir, 'dir'), mkdir(famDir); end

ck = load(bundlePath, 'Fcell', 'Ta_multiples_of_pi', 'DV_cap_nd');
if ~isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) || ck.DV_cap_nd ~= DV_cap_expected
    warning('[c32_fast:repack] %s bundle config mismatch -- skipping repack.', fam);
    return;
end

nRung = numel(ck.Fcell);
for r = 1:nRung
    rpath = fullfile(famDir, sprintf('rung%02d.mat', r));
    if isfile(rpath), continue; end
    F = ck.Fcell{r}; %#ok<NASGU>
    save(rpath, 'F', '-v7.3');
end
clear ck   % discard the full bundle now that all 12 small files are on disk
fprintf('[c32_fast:repack]   %s repacked into %d small files.\n', fam, nRung);
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
