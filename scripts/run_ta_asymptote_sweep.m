%% RUN_TA_ASYMPTOTE_SWEEP
% Test the chaotic-sea asymptotic hypothesis (Ross, email 2026-07-23):
% for an unstable-unstable family pair, does the direct proxy connection
% cost decay to zero as the one-sided propagation horizon Ta -> infinity,
% at FIXED maneuver budget DV_a?
%
% Scope (locked in for this first pass):
%   - ONE pair only: Distant Prograde Orbit (DPO) <-> 3:1 unstable resonant
%     (R31-U). Both endpoints unstable; DPO->R31-U is the pair Shane
%     flagged (direct proxy ~174 m/s at Ta=pi in the paper's Table 7; the
%     C11a relay path already gets it to 59 m/s, so a much cheaper direct
%     connection is plausible at larger Ta).
%   - DV_cap_nd stays at the paper's nominal 0.2 (Table 3) for now. We are
%     NOT reducing it yet -- first see how memory/runtime behave on the
%     HPC for the fattest (largest-Ta) atlas build, then decide.
%   - Grid, seed spacing, fan resolution, ODE tolerances: unchanged from
%     Table 3 / the paper's reference configuration. Only Tmax (Ta) moves.
%
% Build order: LARGEST Ta first (real propagation via atlas_prepare_or_load),
% then every smaller Ta rung is obtained by atlas_derive_subset() filtering
% down from that one fat atlas (no re-integration -- atlas_derive_subset can
% only shrink Tmax/DV_cap, never expand, so the fat build must come first).
%
% Three quantities are tracked per rung, all read off the SAME winning
% voxel (the one that minimizes the total proxy cost, matching the paper's
% own eq. 48 convention):
%   DVlb_mps     -- base turning cost (dv_min_A + dv_min_B). This is the
%                   quantity the hypothesis predicts decays to ~0.
%   DVpatch_mps  -- voxel-patch correction (eq. 41). Depends only on local
%                   speed at the winning voxel and the fixed dtheta=1 deg
%                   grid, so it should NOT trend with Ta. We also record
%                   the winning voxel's (x,y,theta) directly, so "is patch
%                   roughly constant because it's literally the same
%                   voxel winning every time" can be checked exactly
%                   rather than only inferred.
%   DVproxy_mps  -- DVlb + DVpatch, the paper's official metric (eq. 44),
%                   kept for continuity with Table 7 / Fig. 8 and as the
%                   selection criterion (argmin) for which voxel we report.
%
% Expand later: add more pairs to PAIRS (unstable-unstable stress cases,
% stable-unstable pairs to estimate DV_min,s, and the 3 stable-stable
% pairs) once this single-pair run is validated.

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

% Pairs to test. Each entry is a 2-element cell of family-name strings
% (must match src/cr3bp_family_ic.m naming, same strings used throughout
% scripts/run_overlap_dv_tmax_sweep.m etc.).
PAIRS = { ...
    {'Distant Prograde Orbit', 'Resonant 3to1 Unstable'} ...
};

% Ta ladder, expressed as multiples of pi (nd). Build order is largest
% first (real propagation), then atlas_derive_subset() peels off every
% smaller rung from that one fat atlas.
Ta_multiples_of_pi = sort([1 2 4 8 16], 'descend');

% Parallel workers for the (seed x heading) propagation inside each fat
% atlas build. This is orthogonal to the rung/pair loop below, which is
% cheap and runs serially.
N_WORKERS = 4;

CHECKPOINT_FILE = fullfile(repoRoot, 'ta_asymptote_results', 'checkpoint.mat');
OUTPUT_DIR      = fullfile(repoRoot, 'ta_asymptote_results');

% ══════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION  (Table 3 / paper reference values, unchanged)
% ══════════════════════════════════════════════════════════════════════════
cfg = atlas_cfg_defaults();

cfg.grid.dx               = 0.001;
cfg.grid.dy               = 0.001;
cfg.grid.dtheta           = deg2rad(1);
cfg.seed.ds_seed          = 0.01;
cfg.fan.DV_cap_nd         = 0.2;      % TODO: revisit after HPC memory check on fat build
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

% Suppress figure/file output from sub-functions we are not using here.
cfg.io.save_figs   = false;
cfg.io.save_fig    = false;
cfg.io.fig_visible = 'off';

if exist('atlas_cfg_validate', 'file') == 2
    atlas_cfg_validate(cfg);
end

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;

if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

% ── Derive the unique family list actually needed ──────────────────────────
families = {};
for p = 1:numel(PAIRS)
    for k = 1:2
        if ~any(strcmp(families, PAIRS{p}{k}))
            families{end+1} = PAIRS{p}{k}; %#ok<AGROW>
        end
    end
end
nFam  = numel(families);
nRung = numel(Ta_multiples_of_pi);

fprintf('[ta_sweep] Families: %s\n', strjoin(families, ', '));
fprintf('[ta_sweep] Ta ladder (xpi): %s\n', mat2str(Ta_multiples_of_pi));
fprintf('[ta_sweep] DV_cap_nd = %.4g (%.1f m/s one-sided)\n', ...
    cfg.fan.DV_cap_nd, cfg.fan.DV_cap_nd * VU_mps);

% ══════════════════════════════════════════════════════════════════════════
%  LOAD / INITIALISE CHECKPOINT
% ══════════════════════════════════════════════════════════════════════════
results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'build_source', {}, 'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});
rung_done = false(nRung, 1);

if isfile(CHECKPOINT_FILE)
    fprintf('[ta_sweep] Checkpoint found -- loading...\n');
    ck = load(CHECKPOINT_FILE, 'results', 'rung_done', 'Ta_multiples_of_pi');
    if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi)
        results   = ck.results;
        rung_done = ck.rung_done;
        fprintf('[ta_sweep] %d/%d rungs already done.\n', sum(rung_done), nRung);
    else
        warning('[ta_sweep] Checkpoint Ta ladder mismatch -- starting fresh.');
    end
end

% ══════════════════════════════════════════════════════════════════════════
%  PARALLEL POOL (for the fat-atlas build's seed x heading propagation)
% ══════════════════════════════════════════════════════════════════════════
if N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[ta_sweep] Started parpool with %d workers.\n', N_WORKERS);
    end
end

% ══════════════════════════════════════════════════════════════════════════
%  BUILD THE FAT ATLAS  (largest Ta, real propagation, cached by fingerprint)
% ══════════════════════════════════════════════════════════════════════════
grid3 = atlas_grid_make(cfg);

Tmax_fat = Ta_multiples_of_pi(1) * pi;
cfg_fat  = cfg;
cfg_fat.propag.Tmax = Tmax_fat;

fprintf('\n[ta_sweep] ══ Building FAT atlas(es) at Ta = %d*pi = %.4f nd (%.2f days) ══\n', ...
    Ta_multiples_of_pi(1), Tmax_fat, Tmax_fat * TU_days);

S_fat = struct();
for i = 1:nFam
    tFam = tic;
    fprintf('[ta_sweep]   family %d/%d: %s ...\n', i, nFam, families{i});
    [S_fat.(local_fieldkey(families{i})), cacheInfo] = atlas_prepare_or_load(families{i}, cfg_fat, grid3);
    fprintf('[ta_sweep]   %s done in %.1fs (cache %s).\n', ...
        families{i}, toc(tFam), local_hitstr(cacheInfo));
end

% ══════════════════════════════════════════════════════════════════════════
%  MAIN LOOP: for each Ta rung (largest -> smallest), derive atlases and
%  compute the pairwise proxy metrics.
% ══════════════════════════════════════════════════════════════════════════
for r = 1:nRung

    if rung_done(r)
        fprintf('[ta_sweep] Rung %d/%d (Ta=%d*pi) already done -- skipping.\n', ...
            r, nRung, Ta_multiples_of_pi(r));
        continue;
    end

    Tmax_r = Ta_multiples_of_pi(r) * pi;
    fprintf('\n[ta_sweep] ── Rung %d/%d: Ta = %d*pi = %.4f nd (%.2f days) ──\n', ...
        r, nRung, Ta_multiples_of_pi(r), Tmax_r, Tmax_r * TU_days);
    tRung = tic;

    cfg_r = cfg;
    cfg_r.propag.Tmax = Tmax_r;

    % ── Get (or derive) the atlas for each needed family at this rung ──────
    Ssub = struct();
    for i = 1:nFam
        key = local_fieldkey(families{i});
        if r == 1
            % Top rung IS the fat atlas -- use directly, no filtering needed.
            Ssub.(key) = S_fat.(key);
        else
            Ssub.(key) = atlas_derive_subset(S_fat.(key), cfg_r);
        end
    end

    % ── Footprints (per family, reused across all pairs at this rung) ─────
    Fp = struct();
    for i = 1:nFam
        key = local_fieldkey(families{i});
        Fp.(key) = local_compute_footprint(Ssub.(key), grid3, VU_mps, TU_days);
    end

    % ── Pairwise proxy metrics ──────────────────────────────────────────────
    for p = 1:numel(PAIRS)
        famA = PAIRS{p}{1};
        famB = PAIRS{p}{2};
        FA = Fp.(local_fieldkey(famA));
        FB = Fp.(local_fieldkey(famB));

        m = local_run_pair(FA, FB, grid3, cfg_r, VU_mps);

        row = struct( ...
            'pairA', famA, 'pairB', famB, ...
            'Ta_nd', Tmax_r, 'Ta_days', Tmax_r * TU_days, ...
            'build_source', local_ternary(r==1, 'fat', 'derived'), ...
            'DVproxy_mps', m.dv_proxy_mps, 'DVlb_mps', m.dvlb_mps, ...
            'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
            'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, ...
            'thc_deg', m.thc_deg);

        results(end+1) = row; %#ok<AGROW>

        fprintf('[ta_sweep]   %s -> %s : DVproxy=%.3f m/s  DVlb=%.3f m/s  DVpatch=%.3f m/s  TOF=%.2f d  vox=%d\n', ...
            famA, famB, m.dv_proxy_mps, m.dvlb_mps, m.dvpatch_mps, m.tof_days, m.voxelId);
    end

    rung_done(r) = true;

    % ── Checkpoint after every rung ─────────────────────────────────────────
    local_save_checkpoint(CHECKPOINT_FILE, results, rung_done, Ta_multiples_of_pi);
    fprintf('[ta_sweep] Rung %d/%d done in %.1fs -- checkpoint saved.\n', r, nRung, toc(tRung));

end

clear S_fat Ssub Fp   % free the fat atlases now that all rungs are derived

% ══════════════════════════════════════════════════════════════════════════
%  SAVE FINAL RESULTS
% ══════════════════════════════════════════════════════════════════════════
T = struct2table(results);
T = sortrows(T, {'pairA', 'pairB', 'Ta_nd'});

outMat = fullfile(OUTPUT_DIR, 'ta_asymptote_results.mat');
outCsv = fullfile(OUTPUT_DIR, 'ta_asymptote_results.csv');
save(outMat, 'T', 'Ta_multiples_of_pi', 'PAIRS', 'cfg');
writetable(T, outCsv);

fprintf('\n[ta_sweep] ══════════ DONE ══════════\n');
fprintf('  Results table : %s\n', outMat);
fprintf('  CSV           : %s\n', outCsv);

% ── Quick diagnostic plot per pair ──────────────────────────────────────────
for p = 1:numel(PAIRS)
    famA = PAIRS{p}{1};  famB = PAIRS{p}{2};
    mask = strcmp(T.pairA, famA) & strcmp(T.pairB, famB);
    Tp = sortrows(T(mask, :), 'Ta_nd');
    if height(Tp) < 2, continue; end

    fig = figure('Visible', 'off');
    semilogx(Tp.Ta_days, Tp.DVlb_mps, '-o', 'DisplayName', 'DV_{lb} (base)'); hold on;
    semilogx(Tp.Ta_days, Tp.DVpatch_mps, '-s', 'DisplayName', 'DV_{patch}');
    semilogx(Tp.Ta_days, Tp.DVproxy_mps, '-^', 'DisplayName', 'DV_{proxy} (total)');
    xlabel('T_a [days]'); ylabel('\DeltaV [m/s]');
    title(sprintf('%s \\rightarrow %s', famA, famB), 'Interpreter', 'tex');
    legend('Location', 'best'); grid on;

    figPath = fullfile(OUTPUT_DIR, sprintf('ta_asymptote_%s_%s.png', ...
        local_fieldkey(famA), local_fieldkey(famB)));
    print(fig, figPath, '-dpng', '-r150');
    close(fig);
    fprintf('  Plot          : %s\n', figPath);
end

% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function k = local_fieldkey(familyName)
%LOCAL_FIELDKEY  Turn a family-name string into a valid struct field name.
k = matlab.lang.makeValidName(familyName);
end

function s = local_hitstr(cacheInfo)
if cacheInfo.hit, s = 'HIT'; else, s = 'BUILT'; end
end

function out = local_ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function local_save_checkpoint(fpath, results, rung_done, Ta_multiples_of_pi)
tmp = [fpath '.tmp'];
save(tmp, 'results', 'rung_done', 'Ta_multiples_of_pi', '-v7.3');
if isfile(fpath), delete(fpath); end
movefile(tmp, fpath);
end

% ─────────────────────────────────────────────────────────────────────────────
function F = local_compute_footprint(S, grid3, VU_mps, TU_days)
%LOCAL_COMPUTE_FOOTPRINT  Compact per-voxel (dv_min, t_mean) summary for FRS/BRS.
% Identical formulation to run_overlap_dv_tmax_sweep.m's local_compute_footprint.

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
%
% Selection convention matches the paper's eq. (48): the winning voxel is
% the one that minimizes the TOTAL proxy cost (DVlb + DVpatch), not DVlb
% alone. DVlb and DVpatch reported here are therefore the decomposition of
% that SAME voxel -- by construction DVproxy == DVlb + DVpatch exactly.
% The voxel center (xc,yc,thc) is also returned so "is this literally the
% same voxel winning at every Ta" can be checked directly rather than only
% inferred from DVpatch staying flat.

m = struct('voxelId', NaN, 'dv_proxy_mps', NaN, 'dvlb_mps', NaN, ...
    'dvpatch_mps', NaN, 'tof_days', NaN, 'xc', NaN, 'yc', NaN, 'thc_deg', NaN);

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

m.voxelId     = idsO(iWin);
m.dv_proxy_mps = dv_proxy(iWin);
m.dvlb_mps     = dv_lb_vec(iWin);
m.dvpatch_mps  = dv_patch_vec(iWin);
m.tof_days     = t_mean_A(iWin) + t_mean_B(iWin);
m.xc           = grid3.x_centers(ix(iWin));
m.yc           = grid3.y_centers(iy(iWin));
m.thc_deg      = rad2deg(grid3.th_centers(it(iWin)));
end
