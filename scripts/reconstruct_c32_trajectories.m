%% RECONSTRUCT_C32_TRAJECTORIES
% For all 12 C32 pairs (9 unstable + 3 stable), at the largest-Ta rung
% already computed, reconstructs the actual winning trajectory on BOTH
% sides and reports the individual ΔV contribution from each side
% (dv_min_A from C32's forward-reachable side, dv_min_B from the
% target's side) -- answering exactly what fraction of DVlb each side
% contributes, and giving the real (x,y,theta,t) path for plotting.
%
% WHY THIS IS CHEAP: uses ONLY already-cached, already-built fat atlases
% (via atlas_prepare_or_load, which cache-hits instantly if the config
% fingerprint matches -- see atlas_cache_fingerprint.m). No new
% integration happens here. C32 needs its atlas loaded twice (once per
% config below, since the unstable and stable tests used different
% Tmax/DV_cap_nd); each partner needs its atlas loaded once.
%
% HOW THE LOOKUP WORKS: the footprint files only keep the MINIMUM dv per
% voxel (via accumarray(...,@min)), discarding which seed/heading
% achieved it. Rather than recomputing that aggregation globally (which
% would mean sorting arrays with up to ~1.3 BILLION rows), this does a
% TARGETED lookup: for the one specific winning voxel ID we already know
% from the results CSV, filter the raw Step4 rows for exact matches
% (a single vectorized comparison, not a sort), then take whichever
% matching seed/heading has the lowest dv. Once identified, the full
% discretized trajectory for that (iSeed, iHead) is extracted directly
% from the raw voxel-crossing log (Step4 rows already record every voxel
% crossed over time for every seed/heading) -- no re-propagation needed.
%
% CONFIG NOTE: cfg overrides below are copied EXACTLY from the scripts
% that built these atlases (run_ta_asymptote_sweep_full.m for the
% unstable-pairs config, run_ta_asymptote_c32_stable_extended.m for the
% stable-pairs config) -- atlas_cache_fingerprint.m hashes these fields,
% so any mismatch means a cache MISS and a very expensive rebuild.

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

% Winning VoxelId at the largest-Ta rung, taken directly from the
% already-computed results CSVs (ta_asymptote_c32_unstable_results.csv
% and ta_asymptote_c32_stable_extended_results.csv).
unstableJobs = { ...
    'Cycler 11a',             1424675512; ...
    'Cycler 11b',             1945379346; ...
    'Cycler 21',              1573390840; ...
    'Distant Prograde Orbit', 1400645807; ...
    'Lyapunov L1',            1116060712; ...
    'Lyapunov L2',            1132506744; ...
    'Resonant 2to1 Unstable',  238114089; ...
    'Resonant 3to1 Unstable', 1155926940; ...
    'Resonant 5to2 Unstable',    9247648; ...
};

stableJobs = { ...
    'Resonant 2to1 Stable', 1190050053; ...
    'Resonant 3to1 Stable',  804330715; ...
    'Resonant 5to2 Stable',  355961646; ...
};

% ══════════════════════════════════════════════════════════════════════════
%  CONFIG A -- unstable-pairs test (Tmax=16pi, DV_cap_nd=0.2)
%  Must match run_ta_asymptote_sweep_full.m exactly for a cache hit.
% ══════════════════════════════════════════════════════════════════════════
cfgA = atlas_cfg_defaults();
cfgA.grid.dx               = 0.001;
cfgA.grid.dy               = 0.001;
cfgA.grid.dtheta           = deg2rad(1);
cfgA.seed.ds_seed          = 0.01;
cfgA.fan.DV_cap_nd         = 0.2;
cfgA.fan.dtheta_fan        = deg2rad(0.5);
cfgA.propag.absTol         = 1e-8;
cfgA.propag.relTol         = 1e-8;
cfgA.propag.v2tol          = 1e-8;
cfgA.log.step_len_factor   = 0.75;
cfgA.log.maxstep_factor    = 2;
cfgA.cache.enable  = true;
cfgA.cache.dir     = fullfile(repoRoot, 'atlas_cache');
cfgA.cache.rebuild = false;
cfgA.par.enable    = false;   % no parpool needed -- pure lookup, no fresh integration
Ta_multiples_of_pi_A = sort(2.^((-3:8)/2), 'descend');
cfgA.propag.Tmax = Ta_multiples_of_pi_A(1) * pi;   % 16pi

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

grid3A = atlas_grid_make(cfgA);
grid3B = atlas_grid_make(cfgB);
VU_mps  = cfgA.units.VU_mps;   % same unit system for both configs
TU_days = cfgA.units.TU_days;

summary = struct('partner', {}, 'config', {}, 'voxelId', {}, ...
    'dv_A_mps', {}, 'dv_B_mps', {}, 'DVlb_check_mps', {});

% ══════════════════════════════════════════════════════════════════════════
%  PART 1 -- C32's side, config A (unstable partners)
% ══════════════════════════════════════════════════════════════════════════
fprintf('[reconstruct] Loading Cycler 32 fat atlas, config A (Tmax=16pi, cap=0.2)...\n');
tA = tic;
[S_c32_A, cacheInfoA] = atlas_prepare_or_load(HUB_FAMILY, cfgA, grid3A);
fprintf('[reconstruct]   %s in %.1fs.\n', local_hitstr(cacheInfoA), toc(tA));

dvA_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:size(unstableJobs, 1)
    partner = unstableJobs{i,1};
    voxelId = unstableJobs{i,2};
    [found, dv_mps, ~, iSeed_w, iHead_w, half_w, traj] = ...
        local_find_winner(S_c32_A, grid3A, VU_mps, TU_days, voxelId, 'FRS');
    if found
        fprintf('[reconstruct]   C32 side vs %-25s: dv_A=%.4f m/s (seed=%d head=%d half=%s)\n', ...
            partner, dv_mps, iSeed_w, iHead_w, half_w);
        local_save_traj(OUTPUT_DIR, partner, 'C32side', traj);
    else
        warning('[reconstruct] C32 side vs %s: NO MATCH FOUND for voxel %d.', partner, voxelId);
        dv_mps = NaN;
    end
    dvA_map(partner) = dv_mps;
end
clear S_c32_A
fprintf('[reconstruct] Cycler 32 (config A) done, atlas discarded.\n\n');

% ══════════════════════════════════════════════════════════════════════════
%  PART 2 -- C32's side, config B (stable partners)
% ══════════════════════════════════════════════════════════════════════════
fprintf('[reconstruct] Loading Cycler 32 fat atlas, config B (Tmax=32pi, cap=0.15)...\n');
tB = tic;
[S_c32_B, cacheInfoB] = atlas_prepare_or_load(HUB_FAMILY, cfgB, grid3B);
fprintf('[reconstruct]   %s in %.1fs.\n', local_hitstr(cacheInfoB), toc(tB));

for i = 1:size(stableJobs, 1)
    partner = stableJobs{i,1};
    voxelId = stableJobs{i,2};
    [found, dv_mps, ~, iSeed_w, iHead_w, half_w, traj] = ...
        local_find_winner(S_c32_B, grid3B, VU_mps, TU_days, voxelId, 'FRS');
    if found
        fprintf('[reconstruct]   C32 side vs %-25s: dv_A=%.4f m/s (seed=%d head=%d half=%s)\n', ...
            partner, dv_mps, iSeed_w, iHead_w, half_w);
        local_save_traj(OUTPUT_DIR, partner, 'C32side', traj);
    else
        warning('[reconstruct] C32 side vs %s: NO MATCH FOUND for voxel %d.', partner, voxelId);
        dv_mps = NaN;
    end
    dvA_map(partner) = dv_mps;
end
clear S_c32_B
fprintf('[reconstruct] Cycler 32 (config B) done, atlas discarded.\n\n');

% ══════════════════════════════════════════════════════════════════════════
%  PART 3 -- each partner's own side (BRS role), one atlas load each
% ══════════════════════════════════════════════════════════════════════════
allJobs = [unstableJobs; stableJobs];
configLabel = [repmat({'unstable'}, size(unstableJobs,1), 1); repmat({'stable'}, size(stableJobs,1), 1)];

for i = 1:size(allJobs, 1)
    partner = allJobs{i,1};
    voxelId = allJobs{i,2};
    isStable = strcmp(configLabel{i}, 'stable');
    if isStable
        cfg_i = cfgB; grid3_i = grid3B;
    else
        cfg_i = cfgA; grid3_i = grid3A;
    end

    fprintf('[reconstruct] Loading %s fat atlas (%s config)...\n', partner, configLabel{i});
    tP = tic;
    [S_p, cacheInfoP] = atlas_prepare_or_load(partner, cfg_i, grid3_i);
    fprintf('[reconstruct]   %s in %.1fs.\n', local_hitstr(cacheInfoP), toc(tP));

    [found, dv_mps, ~, iSeed_w, iHead_w, half_w, traj] = ...
        local_find_winner(S_p, grid3_i, VU_mps, TU_days, voxelId, 'BRS');
    if found
        fprintf('[reconstruct]   %-25s side: dv_B=%.4f m/s (seed=%d head=%d half=%s)\n', ...
            partner, dv_mps, iSeed_w, iHead_w, half_w);
        local_save_traj(OUTPUT_DIR, partner, 'targetside', traj);
    else
        warning('[reconstruct] %s side: NO MATCH FOUND for voxel %d.', partner, voxelId);
        dv_mps = NaN;
    end

    dv_A = dvA_map(partner);
    summary(end+1) = struct('partner', partner, 'config', configLabel{i}, ...
        'voxelId', voxelId, 'dv_A_mps', dv_A, 'dv_B_mps', dv_mps, ...
        'DVlb_check_mps', dv_A + dv_mps); %#ok<AGROW>

    clear S_p
end

T = struct2table(summary);
outCsv = fullfile(OUTPUT_DIR, 'c32_trajectory_contribution_split.csv');
writetable(T, outCsv);

fprintf('\n[reconstruct] ══════════ DONE ══════════\n');
fprintf('  Contribution split table : %s\n', outCsv);
fprintf('  Trajectory CSVs          : %s\n', OUTPUT_DIR);
disp(T);

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function s = local_hitstr(cacheInfo)
if cacheInfo.hit, s = 'loaded from cache'; else, s = 'built (NOT a cache hit -- check cfg matches!)'; end
end

% ─────────────────────────────────────────────────────────────────────────────
function local_save_traj(outDir, partner, tag, traj)
key = matlab.lang.makeValidName(partner);
fname = fullfile(outDir, sprintf('%s_%s.csv', key, tag));
if isempty(traj)
    warning('[reconstruct] Empty trajectory for %s (%s) -- not saved.', partner, tag);
    return;
end
writetable(traj, fname);
end

% ─────────────────────────────────────────────────────────────────────────────
function [found, dv_mps, t_days, iSeed_w, iHead_w, half_w, traj] = ...
    local_find_winner(S, grid3, VU_mps, TU_days, targetVoxelId, side)
%LOCAL_FIND_WINNER  Targeted lookup: given a specific voxel ID already
% known to be the pairwise-optimal choice (from the results CSV), find
% which (iSeed, iHead, half) achieved the minimum dv there, on either the
% forward (FRS) or mirrored-backward (BRS) side, then extract that
% trajectory's full voxel-crossing history for plotting.
%
% Deliberately does NOT recompute the global per-voxel argmin (which
% would need sorting up to ~1.3 billion rows) -- filters for the ONE
% target voxel ID directly, a single vectorized comparison.

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

    % Stay in compact integer types as long as possible -- converting a
    % billion-row uint16 array to double quadruples its memory footprint
    % for no reason; sub2ind and equality comparisons work directly on
    % integer types. Explicit clears below release each big temporary as
    % soon as it's no longer needed rather than waiting on the next loop
    % iteration's reassignment.
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
