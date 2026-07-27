%% RECONSTRUCT_R52S_ONLY
% Standalone re-run of JUST the one piece that OOM'd in
% reconstruct_c32_trajectories.m: Resonant 5to2 Stable's target-side
% (BRS) lookup and trajectory -- the largest atlas on disk (~11GB
% compressed), the last of 12 pairs processed, and the only one that
% failed. The other 11 pairs' outputs are already saved from the earlier
% successful run; this just fills in the missing one and appends it to
% the same summary CSV.

clear; clc;

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

OUTPUT_DIR = fullfile(repoRoot, 'c32_trajectory_reconstruction');
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

HUB_FAMILY = 'Cycler 32';
PARTNER = 'Resonant 5to2 Stable';
VOXEL_ID = 355961646;

% ══════════════════════════════════════════════════════════════════════════
%  CONFIG B -- stable-pairs test (Tmax=32pi, DV_cap_nd=0.15)
%  Must match run_ta_asymptote_c32_stable_extended.m exactly for a cache hit.
% ══════════════════════════════════════════════════════════════════════════
cfgB = atlas_cfg_defaults();
cfgB.grid.dx               = 0.001;
cfgB.grid.dy               = 0.001;
cfgB.grid.dtheta           = deg2rad(1);
cfgB.seed.ds_seed          = 0.01;
cfgB.fan.DV_cap_nd         = 0.15;
cfgB.fan.dtheta_fan        = deg2rad(0.5);
cfgB.propag.absTol         = 1e-8;
cfgB.propag.relTol         = 1e-8;
cfgB.propag.v2tol          = 1e-8;
cfgB.log.step_len_factor   = 0.75;
cfgB.log.maxstep_factor    = 2;
cfgB.cache.enable  = true;
cfgB.cache.dir     = fullfile(repoRoot, 'atlas_cache');
cfgB.cache.rebuild = false;
cfgB.par.enable    = false;
Ta_multiples_of_pi_B = sort(2.^linspace(-1.5, 5, 24), 'descend');
cfgB.propag.Tmax = Ta_multiples_of_pi_B(1) * pi;   % 32pi

grid3B  = atlas_grid_make(cfgB);
VU_mps  = cfgB.units.VU_mps;
TU_days = cfgB.units.TU_days;

% ── C32's side (FRS) -- recomputed here rather than trusting a value
% carried over from the earlier crashed run, since that run's dvA_map was
% only ever in memory and never written to disk. ────────────────────────
fprintf('[r52s_only] Loading Cycler 32 fat atlas (Tmax=32pi, cap=0.15)...\n');
tC = tic;
[S_c32, cacheInfoC] = atlas_prepare_or_load(HUB_FAMILY, cfgB, grid3B);
fprintf('[r52s_only]   %s in %.1fs.\n', local_hitstr(cacheInfoC), toc(tC));

[foundA, dv_A_mps, ~, iSeedA, iHeadA, halfA, trajA] = ...
    local_find_winner(S_c32, grid3B, VU_mps, TU_days, VOXEL_ID, 'FRS');
if ~foundA
    error('[r52s_only] C32 side: NO MATCH FOUND for voxel %d.', VOXEL_ID);
end
fprintf('[r52s_only]   C32 side: dv_A=%.4f m/s (seed=%d head=%d half=%s)\n', ...
    dv_A_mps, iSeedA, iHeadA, halfA);
key = matlab.lang.makeValidName(PARTNER);
writetable(trajA, fullfile(OUTPUT_DIR, sprintf('%s_C32side.csv', key)));
clear S_c32

% ── target's side (BRS) ─────────────────────────────────────────────────
fprintf('[r52s_only] Loading %s fat atlas (Tmax=32pi, cap=0.15)...\n', PARTNER);
tP = tic;
[S_p, cacheInfoP] = atlas_prepare_or_load(PARTNER, cfgB, grid3B);
fprintf('[r52s_only]   %s in %.1fs.\n', local_hitstr(cacheInfoP), toc(tP));

[found, dv_mps, ~, iSeed_w, iHead_w, half_w, traj] = ...
    local_find_winner(S_p, grid3B, VU_mps, TU_days, VOXEL_ID, 'BRS');

if found
    fprintf('[r52s_only]   %s side: dv_B=%.4f m/s (seed=%d head=%d half=%s)\n', ...
        PARTNER, dv_mps, iSeed_w, iHead_w, half_w);
    writetable(traj, fullfile(OUTPUT_DIR, sprintf('%s_targetside.csv', key)));
else
    error('[r52s_only] NO MATCH FOUND for voxel %d -- something is wrong (config mismatch? voxel ID typo?).', VOXEL_ID);
end

% ── append/update this one row in the shared summary CSV ──────────────────
summaryCsv = fullfile(OUTPUT_DIR, 'c32_trajectory_contribution_split.csv');
row = table({PARTNER}, {'stable'}, VOXEL_ID, dv_A_mps, dv_mps, dv_A_mps + dv_mps, ...
    'VariableNames', {'partner', 'config', 'voxelId', 'dv_A_mps', 'dv_B_mps', 'DVlb_check_mps'});

if isfile(summaryCsv)
    T = readtable(summaryCsv);
    T(strcmp(T.partner, PARTNER), :) = [];   % replace any stale/partial row for this partner
    T = [T; row];
else
    T = row;
end
writetable(T, summaryCsv);

fprintf('\n[r52s_only] ══════════ DONE ══════════\n');
fprintf('  Updated summary : %s\n', summaryCsv);
disp(row);

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  (identical to reconstruct_c32_trajectories.m)
% ══════════════════════════════════════════════════════════════════════════

function s = local_hitstr(cacheInfo)
if cacheInfo.hit, s = 'loaded from cache'; else, s = 'built (NOT a cache hit -- check cfg matches!)'; end
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
        bestDv  = thisMin;
        found   = true;
        dv_mps  = thisMin;
        t_days  = abs(tCand(iLoc)) * TU_days;
        iSeed_w = iSeedCand(iLoc);
        iHead_w = iHeadCand(iLoc);
        half_w  = half;
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
