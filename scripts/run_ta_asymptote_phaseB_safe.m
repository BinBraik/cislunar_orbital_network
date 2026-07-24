%% RUN_TA_ASYMPTOTE_PHASEB_SAFE
% Memory-safe Phase B for the full 78-pair Ta-asymptote sweep.
%
% CONTEXT: the per-family footprint cache files written by
% run_ta_asymptote_sweep_full.m and run_ta_asymptote_c32_hub.m turned out
% to be gigabytes each at Ta=16pi (NOT the ~5-25MB the original codebase's
% comments assumed at Ta=pi -- trajectories touch far more distinct voxels
% over 218 days). That makes BOTH of the following unsafe at this Ta:
%   (a) holding all (13 families x 12 rungs) footprints in memory at once
%       (the full-sweep script's Fp cache -- could be 300-500GB+), and
%   (b) parfor-ing the 78-pair loop, since each worker is a separate OS
%       process and gigabyte-sized sliced footprints get duplicated per
%       worker, not shared.
%
% This script assumes Phase A is already complete (or complete enough for
% the pairs you care about) -- i.e. it ONLY reads the shared footprint
% cache, never builds atlases -- and processes ONE RUNG AT A TIME:
%   for each rung r:
%     for each of the 13 families: load that family's FULL .mat file,
%       keep ONLY Fcell{r} (this rung), discard the other 11 rungs
%       immediately (clear the loaded struct after extracting).
%     run all 78 pairs for this rung SERIALLY (no parfor -- Phase B's
%       per-pair computation is a fast vectorized intersect/argmin, it
%       was never the bottleneck; Phase A's atlas builds are).
%     save results for this rung, clear this rung's 13 footprints,
%       move to the next rung.
%
% Peak memory is therefore bounded by roughly (one family's full file,
% transiently, during its own load) + (13 families' worth of ONE rung),
% not all 156 (family x rung) combinations at once. Trade-off: each
% family's multi-GB file gets read from disk up to 12 times (once per
% rung) instead of once -- more disk I/O, but safe regardless of how
% large an individual footprint file is.
%
% Checkpointed per rung, same as the other scripts.

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

FAMILIES = { ...
    'Lyapunov L1', 'Lyapunov L2', 'Cycler 21', 'Cycler 11a', 'Cycler 11b', ...
    'Cycler 32', 'Resonant 2to1 Stable', 'Resonant 2to1 Unstable', ...
    'Resonant 3to1 Stable', 'Resonant 3to1 Unstable', ...
    'Resonant 5to2 Stable', 'Resonant 5to2 Unstable', 'Distant Prograde Orbit' ...
};

Ta_multiples_of_pi = sort(2.^((-3:8)/2), 'descend');   % MUST match the cache-writing scripts
DV_cap_nd_expected  = 0.2;

% Same shared cache the two running jobs read/write.
FULL_SWEEP_OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_full_results');
FOOTPRINT_DIR = fullfile(FULL_SWEEP_OUTPUT_DIR, 'footprints');

OUTPUT_DIR   = fullfile(repoRoot, 'ta_asymptote_phaseB_safe_results');
CHECKPOINT_B = fullfile(OUTPUT_DIR, 'checkpoint_phaseB_safe.mat');

if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

% ══════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION  (only cfg.sys / cfg.overlap / grid needed -- no
%  propagation happens in this script)
% ══════════════════════════════════════════════════════════════════════════
cfg = atlas_cfg_defaults();
cfg.grid.dx      = 0.001;
cfg.grid.dy      = 0.001;
cfg.grid.dtheta  = deg2rad(1);
cfg.fan.DV_cap_nd = DV_cap_nd_expected;

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
grid3   = atlas_grid_make(cfg);

nFam  = numel(FAMILIES);
nRung = numel(Ta_multiples_of_pi);

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

% ── Check which families are actually available in the shared cache ────────
famReady = false(nFam, 1);
for i = 1:nFam
    fpath = fullfile(FOOTPRINT_DIR, [local_fieldkey(FAMILIES{i}) '.mat']);
    famReady(i) = isfile(fpath);
end
fprintf('[phaseB_safe] %d/%d families available in shared cache:\n', sum(famReady), nFam);
for i = 1:nFam
    fprintf('  [%s] %s\n', local_ternary(famReady(i), 'x', ' '), FAMILIES{i});
end

readyPairMask = famReady(pairIdxA) & famReady(pairIdxB);
fprintf('[phaseB_safe] %d/%d pairs computable right now (rest need families not yet cached).\n', ...
    sum(readyPairMask), nPairs);

% ══════════════════════════════════════════════════════════════════════════
%  CHECKPOINT
% ══════════════════════════════════════════════════════════════════════════
results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});
rung_done = false(nRung, 1);

if isfile(CHECKPOINT_B)
    ck = load(CHECKPOINT_B, 'results', 'rung_done', 'Ta_multiples_of_pi');
    if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi)
        results   = ck.results;
        rung_done = ck.rung_done;
        fprintf('[phaseB_safe] Checkpoint loaded: %d/%d rungs already done.\n', sum(rung_done), nRung);
    end
end

% ══════════════════════════════════════════════════════════════════════════
%  MAIN LOOP -- one rung at a time, only THIS rung's 13 footprints resident
% ══════════════════════════════════════════════════════════════════════════
for r = 1:nRung

    if rung_done(r)
        fprintf('[phaseB_safe] Rung %d/%d already done -- skipping.\n', r, nRung);
        continue;
    end

    Tmax_r = Ta_multiples_of_pi(r) * pi;
    fprintf('\n[phaseB_safe] ── Rung %d/%d: Ta = %.4gpi (%.2f days) ──\n', ...
        r, nRung, Ta_multiples_of_pi(r), Tmax_r * TU_days);
    tRung = tic;

    % ── Load ONLY this rung's footprint for each ready family ──────────────
    Fr = cell(nFam, 1);
    for i = 1:nFam
        if ~famReady(i), continue; end
        fpath = fullfile(FOOTPRINT_DIR, [local_fieldkey(FAMILIES{i}) '.mat']);
        ck = load(fpath, 'Fcell', 'Ta_multiples_of_pi', 'DV_cap_nd');
        if ~isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) || ...
                ck.DV_cap_nd ~= cfg.fan.DV_cap_nd
            warning('[phaseB_safe] %s cache config mismatch -- skipping this family.', FAMILIES{i});
            clear ck
            continue;
        end
        Fr{i} = ck.Fcell{r};   % keep ONLY this rung
        clear ck               % discard the other 11 rungs immediately
    end

    % ── Serial pair loop (no parfor -- avoid duplicating GB footprints
    %    across 60 worker processes; this loop itself is fast) ─────────────
    nDone = 0;
    for pIdx = 1:nPairs
        iA = pairIdxA(pIdx);  iB = pairIdxB(pIdx);
        if ~famReady(iA) || ~famReady(iB), continue; end
        if isempty(Fr{iA}) || isempty(Fr{iB}), continue; end

        m = local_run_pair(Fr{iA}, Fr{iB}, grid3, cfg, VU_mps);

        row = struct( ...
            'pairA', FAMILIES{iA}, 'pairB', FAMILIES{iB}, ...
            'Ta_nd', Tmax_r, 'Ta_days', Tmax_r * TU_days, ...
            'DVproxy_mps', m.dv_proxy_mps, 'DVlb_mps', m.dvlb_mps, ...
            'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
            'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, 'thc_deg', m.thc_deg);
        results(end+1) = row; %#ok<AGROW>
        nDone = nDone + 1;
    end

    clear Fr   % discard this rung's footprints before loading the next rung

    rung_done(r) = true;
    save(CHECKPOINT_B, 'results', 'rung_done', 'Ta_multiples_of_pi', '-v7.3');
    fprintf('[phaseB_safe] Rung %d/%d: %d/%d pairs computed in %.1fs -- checkpoint saved.\n', ...
        r, nRung, nDone, nPairs, toc(tRung));
end

% ══════════════════════════════════════════════════════════════════════════
%  SAVE FINAL RESULTS
% ══════════════════════════════════════════════════════════════════════════
T = struct2table(results);
T = sortrows(T, {'pairA', 'pairB', 'Ta_nd'});

outMat = fullfile(OUTPUT_DIR, 'ta_asymptote_phaseB_safe_results.mat');
outCsv = fullfile(OUTPUT_DIR, 'ta_asymptote_phaseB_safe_results.csv');
save(outMat, 'T', 'Ta_multiples_of_pi', 'FAMILIES', 'cfg');
writetable(T, outCsv);

fprintf('\n[phaseB_safe] ══════════ DONE ══════════\n');
fprintf('  Rows written  : %d (out of %d possible = %d pairs x %d rungs)\n', ...
    height(T), nPairs*nRung, nPairs, nRung);
fprintf('  Results table : %s\n', outMat);
fprintf('  CSV           : %s\n', outCsv);
fprintf('  Rerun this script any time -- families that finish building later\n');
fprintf('  (in either of the other two jobs) will simply fill in more pairs.\n');

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function k = local_fieldkey(familyName)
k = matlab.lang.makeValidName(familyName);
end

function out = local_ternary(cond, a, b)
if cond, out = a; else, out = b; end
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
