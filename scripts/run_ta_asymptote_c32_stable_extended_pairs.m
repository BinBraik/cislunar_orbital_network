%% RUN_TA_ASYMPTOTE_C32_STABLE_EXTENDED_PAIRS
% Phase-B-equivalent for the extended-Ta (32pi, 24 rungs), reduced-cap
% (DV_cap_nd=0.15) C32-vs-stable test. Reads the per-rung footprint files
% written directly by run_ta_asymptote_c32_stable_extended_family.m --
% those files already ARE the by-rung cache (no monolithic bundle stage,
% no repack needed, unlike the original full-sweep/hub scripts).
%
% Computes both C32-vs-Resonant2to1Stable and C32-vs-Resonant3to1Stable
% curves across all 24 rungs. Run this only after all three families'
% run_ta_asymptote_c32_stable_extended_family.m calls have finished.

clear; clc;

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

HUB_FAMILY = 'Cycler 32';
PARTNER_FAMILIES = {'Resonant 2to1 Stable', 'Resonant 3to1 Stable'};

Ta_multiples_of_pi = sort(2.^linspace(-1.5, 5, 24), 'descend');   % must match the family script
DV_cap_nd_expected = 0.15;

cfg = atlas_cfg_defaults();
cfg.grid.dx       = 0.001;
cfg.grid.dy       = 0.001;
cfg.grid.dtheta   = deg2rad(1);
cfg.fan.DV_cap_nd = DV_cap_nd_expected;

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
grid3   = atlas_grid_make(cfg);
nRung   = numel(Ta_multiples_of_pi);

OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_c32_stable_extended_results');
BYRUNG_DIR = fullfile(OUTPUT_DIR, 'footprints_by_rung');

ALL_NEEDED = [{HUB_FAMILY}, PARTNER_FAMILIES];
famReady = false(size(ALL_NEEDED));
for i = 1:numel(ALL_NEEDED)
    famDir = fullfile(BYRUNG_DIR, local_fieldkey(ALL_NEEDED{i}));
    ready = true;
    for r = 1:nRung
        rpath = fullfile(famDir, sprintf('rung%02d.mat', r));
        if ~isfile(rpath), ready = false; break; end
        ck = load(rpath, 'Ta_multiples_of_pi', 'DV_cap_nd');
        if ~isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) || ck.DV_cap_nd ~= DV_cap_nd_expected
            ready = false; break;
        end
    end
    famReady(i) = ready;
end
readyFamilies   = ALL_NEEDED(famReady);
missingFamilies = ALL_NEEDED(~famReady);

fprintf('[c32_stable_ext_pairs] %d/%d families ready: %s\n', ...
    numel(readyFamilies), numel(ALL_NEEDED), strjoin(readyFamilies, ', '));
if ~isempty(missingFamilies)
    fprintf('[c32_stable_ext_pairs] Not ready, skipped: %s\n', strjoin(missingFamilies, ', '));
end
if ~any(strcmp(readyFamilies, HUB_FAMILY))
    error('[c32_stable_ext_pairs] Hub family %s is not ready -- nothing to do.', HUB_FAMILY);
end

availablePartners = intersect(PARTNER_FAMILIES, readyFamilies, 'stable');

results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});

for r = 1:nRung
    Tmax_r = Ta_multiples_of_pi(r) * pi;
    tRung = tic;
    FHub_r = local_load_rung(HUB_FAMILY, r, BYRUNG_DIR);

    for pIdx = 1:numel(availablePartners)
        Fp_r = local_load_rung(availablePartners{pIdx}, r, BYRUNG_DIR);
        m = local_run_pair(FHub_r, Fp_r, grid3, cfg, VU_mps);

        results(end+1) = struct( ...
            'pairA', HUB_FAMILY, 'pairB', availablePartners{pIdx}, ...
            'Ta_nd', Tmax_r, 'Ta_days', Tmax_r * TU_days, ...
            'DVproxy_mps', m.dv_proxy_mps, 'DVlb_mps', m.dvlb_mps, ...
            'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
            'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, 'thc_deg', m.thc_deg); %#ok<AGROW>
    end
    clear FHub_r
    fprintf('[c32_stable_ext_pairs] Rung %d/%d (Ta=%.2fd) done in %.1fs.\n', r, nRung, Tmax_r * TU_days, toc(tRung));
end

T = struct2table(results);
T = sortrows(T, {'pairB', 'Ta_nd'});

% Flag any pair/rung whose found minimum sits suspiciously close to the
% cap -- that's the signature of the cap clipping the true optimum rather
% than finding it (see the DV_cap_nd=0.15 sizing rationale in the family
% script's header).
capMps = cfg.fan.DV_cap_nd * VU_mps;
nearCap = T.DVlb_mps > 0.95 * capMps;
if any(nearCap)
    warning(['[c32_stable_ext_pairs] %d row(s) have DVlb within 5%% of the ' ...
        'DV_cap_nd=%.3g ceiling (%.1f m/s) -- the cap may be clipping the true ' ...
        'minimum there; consider rerunning with a larger cap.'], ...
        sum(nearCap), cfg.fan.DV_cap_nd, capMps);
end

if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end
outMat = fullfile(OUTPUT_DIR, 'ta_asymptote_c32_stable_extended_results.mat');
outCsv = fullfile(OUTPUT_DIR, 'ta_asymptote_c32_stable_extended_results.csv');
save(outMat, 'T', 'Ta_multiples_of_pi', 'HUB_FAMILY', 'availablePartners', 'missingFamilies', 'cfg');
writetable(T, outCsv);

fprintf('\n[c32_stable_ext_pairs] ══════════ DONE ══════════\n');
fprintf('  Computed pairs : %d/%d\n', numel(availablePartners), numel(PARTNER_FAMILIES));
fprintf('  Results table  : %s\n', outMat);
fprintf('  CSV            : %s\n', outCsv);

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function k = local_fieldkey(familyName)
k = matlab.lang.makeValidName(familyName);
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
