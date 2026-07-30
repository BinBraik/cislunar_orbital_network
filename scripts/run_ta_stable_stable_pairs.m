%% RUN_TA_STABLE_STABLE_PAIRS
% Tests the third leg of the chaotic-sea hypothesis: two STABLE families
% paired against each other should converge, as Ta grows, to
% ΔV_min,s1 + ΔV_min,s2 -- the SUM of each stable orbit's own individual
% floor (already found against multiple unstable hubs: ~125-127 m/s for
% R21-S, ~86-88 m/s for R31-S, ~0.5-3.3 m/s for R52-S), rather than
% collapsing toward zero the way unstable-unstable pairs do.
%
% WHY THIS IS CHEAP: all 3 stable families (R21-S, R31-S, R52-S) are
% already fully built and cached at this exact Ta ladder (32pi, 24
% rungs) and DV_cap_nd=0.15 from the earlier C32 and partner-independence
% runs. This is Phase-B only -- no atlas building, just reading already-
% cached per-rung footprints, same as the fast hub-comparison scripts.
%
% TRAJECTORY RECONSTRUCTION, BUILT IN THIS TIME: the C32-pairs
% reconstruction was done as an awkward separate pass after the fact,
% which cost an extra debugging round (table() shape bug, meeting-point
% truncation, an OOM on the largest atlas). With only 3 pairs and 3
% already-known-cheap atlases here, this script does the reconstruction
% immediately after Phase B, at each pair's largest-Ta rung, using the
% same fixed local_find_winner (compact uint32 arithmetic, truncates at
% the exact meeting-point crossing time) -- and writes both the
% contribution-split summary and each pair's trajectory CSVs
% incrementally, not batched at the end.
%
% PAIR/ROLE BOOKKEEPING: each of the 3 families appears in exactly 2 of
% the 3 pairs (playing "A"/FRS in one, sometimes "B"/BRS in the other),
% so each atlas only needs to be loaded ONCE and can serve both of that
% family's queries before being discarded.

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
% (iA, iB) index pairs into FAMILIES3 -- consistent A/B labeling used
% throughout for both Phase B and the trajectory reconstruction below.
PAIRS = {1, 2; 1, 3; 2, 3};

Ta_multiples_of_pi = sort(2.^linspace(-1.5, 5, 24), 'descend');   % must match the C32/partner-independence runs for cache hits
DV_cap_nd = 0.15;

BYRUNG_DIR    = fullfile(repoRoot, 'ta_asymptote_c32_stable_extended_results', 'footprints_by_rung');
RESULTS_DIR   = fullfile(repoRoot, 'ta_stable_stable_final');
TRAJ_OUT_DIR  = fullfile(repoRoot, 'ta_stable_stable_trajectories');
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end
if ~exist(TRAJ_OUT_DIR, 'dir'), mkdir(TRAJ_OUT_DIR); end

% ══════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION -- must match the stable-pairs runs for cache hits
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
cfg.cache.enable  = true;
cfg.cache.dir     = fullfile(repoRoot, 'atlas_cache');
cfg.cache.rebuild = false;
cfg.par.enable    = false;   % pure cache-read + lookup, no fresh integration
Ta_fat = Ta_multiples_of_pi(1) * pi;
cfg.propag.Tmax = Ta_fat;

grid3   = atlas_grid_make(cfg);
VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;
nRung   = numel(Ta_multiples_of_pi);
nPairs  = size(PAIRS, 1);

fprintf('[stable_stable] %d pairs among %d stable families, %d rungs, DV_cap_nd=%.3g\n', ...
    nPairs, numel(FAMILIES3), nRung, DV_cap_nd);

% ══════════════════════════════════════════════════════════════════════════
%  PHASE B -- pairwise computation, pure cache read (no atlas building)
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[stable_stable] ══════════ PHASE B: %d pairs x %d rungs ══════════\n', nPairs, nRung);

results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});

resultsCsv = fullfile(RESULTS_DIR, 'ta_stable_stable_results.csv');

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
            'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
            'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, 'thc_deg', m.thc_deg); %#ok<AGROW>
    end
    clear Frung

    % Write after every rung, not just once at the end.
    writetable(struct2table(results), resultsCsv);
    fprintf('[stable_stable] Rung %d/%d (Ta=%.2fd) done in %.1fs.\n', r, nRung, Tmax_r * TU_days, toc(tRung));
end

T = struct2table(results);
T = sortrows(T, {'pairA', 'pairB', 'Ta_nd'});
writetable(T, resultsCsv);
save(fullfile(RESULTS_DIR, 'ta_stable_stable_results.mat'), 'T', 'Ta_multiples_of_pi', 'FAMILIES3', 'PAIRS', 'cfg');

fprintf('\n[stable_stable] ══════════ Phase B complete. ══════════\n');
fprintf('  Results table : %s\n', resultsCsv);
disp(T(T.Ta_nd == max(T.Ta_nd), :));

% ══════════════════════════════════════════════════════════════════════════
%  TRAJECTORY RECONSTRUCTION -- at each pair's largest-Ta rung
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[stable_stable] ══════════ TRAJECTORY RECONSTRUCTION ══════════\n');

% For each pair, the winning VoxelId at the largest Ta (last rung in the
% descending ladder = smallest Ta_nd... careful: ladder is DESCENDING, so
% "largest Ta" is Ta_multiples_of_pi(1), i.e. the FIRST rung, not the last).
targetVoxel = containers.Map('KeyType', 'char', 'ValueType', 'double');
pairKey = @(a, b) sprintf('%s|%s', a, b);
for p = 1:nPairs
    iA = PAIRS{p, 1}; iB = PAIRS{p, 2};
    famA = FAMILIES3{iA}; famB = FAMILIES3{iB};
    rowsThisPair = T(strcmp(T.pairA, famA) & strcmp(T.pairB, famB), :);
    [~, iMax] = max(rowsThisPair.Ta_nd);
    targetVoxel(pairKey(famA, famB)) = rowsThisPair.VoxelId(iMax);
end

% Which (pair, side) queries does each family need to serve?
% side = 'FRS' when the family plays role A in that pair, 'BRS' when B.
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

dvParts = containers.Map('KeyType', 'char', 'ValueType', 'double');   % key: "famA|famB|side" -> dv contribution
contribCsv = fullfile(RESULTS_DIR, 'ta_stable_stable_contribution_split.csv');

for f = 1:numel(FAMILIES3)
    familyName = FAMILIES3{f};
    queries = familyQueries(familyName);
    if isempty(queries), continue; end

    fprintf('[stable_stable] Loading %s fat atlas (Tmax=32pi, cap=0.15)...\n', familyName);
    tF = tic;
    [S_f, cacheInfo] = atlas_prepare_or_load(familyName, cfg, grid3);
    fprintf('[stable_stable]   %s in %.1fs.\n', local_hitstr(cacheInfo), toc(tF));

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
            fprintf('[stable_stable]   %s (%s side, vs %s): dv=%.4f m/s (seed=%d head=%d half=%s)\n', ...
                familyName, q.side, q.partner, dv_mps, iSeed_w, iHead_w, half_w);
            key = local_fieldkey(familyName);
            sideTag = 'A'; if strcmp(q.side, 'BRS'), sideTag = 'B'; end
            fname = sprintf('%s_vs_%s_%sside.csv', key, matlab.lang.makeValidName(q.partner), sideTag);
            writetable(traj, fullfile(TRAJ_OUT_DIR, fname));
        else
            warning('[stable_stable] %s (%s side, vs %s): NO MATCH FOUND for voxel %d.', ...
                familyName, q.side, q.partner, q.voxelId);
            dv_mps = NaN;
        end

        dvParts(sprintf('%s|%s|%s', famA, famB, q.side)) = dv_mps;

        % Write the growing contribution table after every single lookup.
        pk = pairKey(famA, famB);
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

fprintf('\n[stable_stable] ══════════ DONE ══════════\n');
fprintf('  Contribution split : %s\n', contribCsv);
fprintf('  Trajectory CSVs    : %s\n', TRAJ_OUT_DIR);

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function k = local_fieldkey(familyName)
k = matlab.lang.makeValidName(familyName);
end

function s = local_hitstr(cacheInfo)
if cacheInfo.hit, s = 'loaded from cache'; else, s = 'built (NOT a cache hit -- check cfg matches!)'; end
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
