%% RUN_TA_ASYMPTOTE_C32_UNSTABLE_FAST
% C32-hub test restricted to the 9 UNSTABLE partners only (LL1, LL2, C21,
% C11a, C11b, R21-U, R31-U, R52-U, DPO). The 3 stable partners (R21-S,
% R31-S, R52-S) are deliberately excluded here -- they need a much larger
% Ta ladder to actually converge (both were still decreasing at Ta=218.6
% days in the combined run), which means genuinely new, longer atlas
% builds for both the stable orbits AND C32 itself (atlas_derive_subset
% can only shrink a cached atlas, never expand it). That's real new
% integration work belonging to a separate, dedicated script -- this one
% stays scoped to the unstable pairs, which already converge cleanly
% within the existing 12-rung/16pi ladder and need nothing beyond what's
% already cached.
%
% Same repack-then-compute design as run_ta_asymptote_c32_hub_fast.m:
% each family's big bundle file (footprints/) gets split into 12 small
% per-rung files (footprints_by_rung/) ONCE, then the pairwise computation
% reads only the small files. Both directories are the SAME shared cache
% the other scripts read/write, so anything already repacked (from the
% earlier combined run) is reused instantly.
%
% Run this once run_ta_asymptote_full_v2.slurm has finished R52-U and DPO
% (the two unstable families that weren't ready yet as of the last check).
% If either is still missing, this script reports it and computes
% whatever it can -- safe and idempotent to rerun once they land.

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

% Unstable partners ONLY -- the 3 stable ones (Resonant *to1/2 Stable) are
% intentionally not listed here; see header.
PARTNER_FAMILIES = { ...
    'Lyapunov L1', 'Lyapunov L2', 'Cycler 21', 'Cycler 11a', 'Cycler 11b', ...
    'Resonant 2to1 Unstable', 'Resonant 3to1 Unstable', 'Resonant 5to2 Unstable', ...
    'Distant Prograde Orbit' ...
};

Ta_multiples_of_pi = sort(2.^((-3:8)/2), 'descend');   % must match the writing scripts
DV_cap_nd_expected  = 0.2;

% Bounded on purpose -- see run_ta_asymptote_c32_hub_fast.m's header for
% the full rationale (worst-case simultaneous large-file reads during the
% one-time repack step, not CPU count).
N_WORKERS = 8;

FULL_SWEEP_OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_full_results');
BUNDLE_DIR = fullfile(FULL_SWEEP_OUTPUT_DIR, 'footprints');            % existing big per-family files
BYRUNG_DIR = fullfile(FULL_SWEEP_OUTPUT_DIR, 'footprints_by_rung');    % shared small per-family-per-rung cache

OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_c32_unstable_results');
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

fprintf('[c32_unstable] %d/%d families ready (bundle exists): %s\n', ...
    numel(readyFamilies), numel(ALL_NEEDED), strjoin(readyFamilies, ', '));
if ~isempty(missingFamilies)
    fprintf('[c32_unstable] Not yet built, skipped this run: %s\n', strjoin(missingFamilies, ', '));
end

if ~any(strcmp(readyFamilies, HUB_FAMILY))
    error('[c32_unstable] Hub family %s is not built yet -- nothing to do.', HUB_FAMILY);
end

% ══════════════════════════════════════════════════════════════════════════
%  PARALLEL POOL  (deliberately small -- see header)
% ══════════════════════════════════════════════════════════════════════════
if N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[c32_unstable] Started parpool with %d workers (bounded by design).\n', N_WORKERS);
    end
end

% ══════════════════════════════════════════════════════════════════════════
%  STEP 1 -- one-time repack: big bundle -> 12 small per-rung files
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[c32_unstable] ══════════ STEP 1: repack bundles into per-rung files ══════════\n');

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
fprintf('[c32_unstable] %d/%d families already repacked; %d need repacking: %s\n', ...
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
fprintf('[c32_unstable] Step 1 complete.\n');

% ══════════════════════════════════════════════════════════════════════════
%  STEP 2 -- fast pairwise computation, rung by rung, small files only
% ══════════════════════════════════════════════════════════════════════════
availablePartners = intersect(PARTNER_FAMILIES, readyFamilies, 'stable');
fprintf('\n[c32_unstable] ══════════ STEP 2: %d pairs x %d rungs ══════════\n', ...
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

    fprintf('[c32_unstable] Rung %d/%d (Ta=%.2fd) done in %.1fs.\n', r, nRung, Tmax_r * TU_days, toc(tRung));
end

T = struct2table(results);
T = sortrows(T, {'pairB', 'Ta_nd'});

outMat = fullfile(OUTPUT_DIR, 'ta_asymptote_c32_unstable_results.mat');
outCsv = fullfile(OUTPUT_DIR, 'ta_asymptote_c32_unstable_results.csv');
save(outMat, 'T', 'Ta_multiples_of_pi', 'HUB_FAMILY', 'availablePartners', 'missingFamilies', 'cfg');
writetable(T, outCsv);

fprintf('\n[c32_unstable] ══════════ DONE ══════════\n');
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
    warning('[c32_unstable:repack] %s bundle config mismatch -- skipping repack.', fam);
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
fprintf('[c32_unstable:repack]   %s repacked into %d small files.\n', fam, nRung);
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
