%% RUN_TA_ASYMPTOTE_C32_STABLE_EXTENDED
% Extended-Ta (32pi, 24 rungs), reduced-cap (DV_cap_nd=0.15) C32-vs-stable
% test. Single MATLAB process, single parpool, plain serial loop over the
% 4 families -- same proven shape as run_ta_asymptote_sweep_full.m, which
% already completed successfully twice on this project. No multi-process
% splitting: that design (one matlab -batch per family, launched
% concurrently from bash) failed silently -- `wait $PID1 $PID2 $PID3
% $PID4` doesn't check exit codes, so when all four background processes
% died, the script sailed past it and only surfaced the failure at the
% very end, after a 16+ hour queue wait for nothing. A single serial
% process either works or throws a real MATLAB error immediately -- no
% hidden failure path.
%
% WHY THIS Ta/CAP COMBINATION: doubles Tmax to 32pi (~437 days) to chase
% the stable-pair asymptote further than the 16pi/218-day C32-hub run
% got (C32-R21S and C32-R31S were still slowly decreasing at Ta=218.6d,
% ~5.4 m/s/rung, not yet plateaued). DV_cap_nd=0.15 (153.5 m/s) was never
% the binding constraint on those curves (130.5 and 91.8 m/s at
% Ta=218.6d) -- lowering it from 0.2 cuts fan search space/memory by 25%
% while keeping ~18%+ margin above the harder (R21S) curve's last known
% value. R52-S has no prior data to check the cap against; the near-cap
% warning at the end flags it if that assumption doesn't hold.
%
% Per-rung save: each rung's footprint is computed and written to its
% own small file immediately, not batched into one big per-family bundle
% -- so nothing is lost if the job is interrupted, and reruns skip
% whatever's already done (checked per rung, not just per family).

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
PARTNER_FAMILIES = {'Resonant 2to1 Stable', 'Resonant 3to1 Stable', 'Resonant 5to2 Stable'};
FAMILIES = [{HUB_FAMILY}, PARTNER_FAMILIES];

Ta_multiples_of_pi = sort(2.^linspace(-1.5, 5, 24), 'descend');   % 0.354pi ... 32pi, 24 rungs
DV_cap_nd = 0.15;

% Single shared pool for the whole serial run -- same worker count that
% already worked cleanly on this exact hardware (700G/61 cores, Owl base
% nodes) for a fresh full-atlas build (job 640852). Overridable via env
% var so a .slurm file can size it to whatever cpus-per-task it actually
% requested, without needing a second copy of this script.
N_WORKERS = 60;
envWorkers = getenv('C32_STABLE_EXT_NWORKERS');
if ~isempty(envWorkers)
    N_WORKERS = str2double(envWorkers);
end

OUTPUT_DIR  = fullfile(repoRoot, 'ta_asymptote_c32_stable_extended_results');
BYRUNG_DIR  = fullfile(OUTPUT_DIR, 'footprints_by_rung');
RESULTS_DIR = fullfile(repoRoot, 'ta_asymptote_c32_stable_extended_final');

% ══════════════════════════════════════════════════════════════════════════
%  BASE CONFIGURATION
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
cfg.par.progress_every    = 1000;   % was the 50-job default -- far too chatty at this scale

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

grid3 = atlas_grid_make(cfg);
nFam  = numel(FAMILIES);
nRung = numel(Ta_multiples_of_pi);

if ~exist(BYRUNG_DIR, 'dir'), mkdir(BYRUNG_DIR); end

fprintf('[c32_stable_ext] Families: %d   Ta rungs: %d   DV_cap_nd: %.3g (%.1f m/s)\n', ...
    nFam, nRung, DV_cap_nd, DV_cap_nd * VU_mps);
fprintf('[c32_stable_ext] Ta ladder (xpi): %s\n', mat2str(Ta_multiples_of_pi, 4));

if N_WORKERS > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', N_WORKERS);
        fprintf('[c32_stable_ext] Started parpool with %d workers.\n', N_WORKERS);
    end
end

% ══════════════════════════════════════════════════════════════════════════
%  PHASE A -- per-family footprints, serial, save-per-rung-immediately
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[c32_stable_ext] ══════════ PHASE A: per-family footprints ══════════\n');

for i = 1:nFam
    familyName = FAMILIES{i};
    key    = local_fieldkey(familyName);
    famDir = fullfile(BYRUNG_DIR, key);
    if ~exist(famDir, 'dir'), mkdir(famDir); end

    allDone = true;
    for r = 1:nRung
        rpath = fullfile(famDir, sprintf('rung%02d.mat', r));
        if ~isfile(rpath), allDone = false; break; end
        ck = load(rpath, 'Ta_multiples_of_pi', 'DV_cap_nd');
        if ~isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) || ck.DV_cap_nd ~= DV_cap_nd
            allDone = false; break;
        end
    end
    if allDone
        fprintf('[c32_stable_ext] (%d/%d) %-25s -- all %d rungs already cached, skipping.\n', ...
            i, nFam, familyName, nRung);
        continue;
    end

    tFam = tic;
    Tmax_fat = Ta_multiples_of_pi(1) * pi;
    fprintf('[c32_stable_ext] (%d/%d) %-25s -- building fat atlas at Ta=%.4gpi (%.2f days)...\n', ...
        i, nFam, familyName, Ta_multiples_of_pi(1), Tmax_fat * TU_days);

    cfg_fat = cfg;
    cfg_fat.propag.Tmax = Tmax_fat;
    [S_fat, cacheInfo] = atlas_prepare_or_load(familyName, cfg_fat, grid3);
    fprintf('[c32_stable_ext]   fat atlas %s in %.1fs.\n', local_hitstr(cacheInfo), toc(tFam));

    for r = 1:nRung
        rpath = fullfile(famDir, sprintf('rung%02d.mat', r));
        if isfile(rpath)
            ck = load(rpath, 'Ta_multiples_of_pi', 'DV_cap_nd');
            if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) && ck.DV_cap_nd == DV_cap_nd
                fprintf('[c32_stable_ext]   %-25s rung %d/%d already cached, skipping.\n', familyName, r, nRung);
                continue;
            end
        end

        Tmax_r = Ta_multiples_of_pi(r) * pi;
        if r == 1
            Ssub = S_fat;   % top rung IS the fat atlas, no filtering needed
        else
            cfg_r = cfg;
            cfg_r.propag.Tmax = Tmax_r;
            Ssub = atlas_derive_subset(S_fat, cfg_r);
        end
        F = local_compute_footprint(Ssub, grid3, VU_mps, TU_days); %#ok<NASGU>
        clear Ssub
        % Re-check/recreate famDir immediately before every save rather
        % than trusting the one-time mkdir at the top of this family's
        % loop -- a companion run (partner-independence) hit "Cannot
        % create rung01.mat because ... does not exist" after an
        % hour-plus atlas build sat between the directory being created
        % and the first save. Cheap insurance against that recurring.
        if ~exist(famDir, 'dir'), mkdir(famDir); end
        save(rpath, 'F', 'Ta_multiples_of_pi', 'DV_cap_nd', '-v7.3');
        clear F
        fprintf('[c32_stable_ext]   %-25s rung %d/%d (Ta=%.2fd) done and saved (%.1fs elapsed).\n', ...
            familyName, r, nRung, Tmax_r * TU_days, toc(tFam));
    end

    clear S_fat   % discard the full atlas -- only the per-rung footprints on disk survive
    fprintf('[c32_stable_ext] (%d/%d) %-25s -- done in %.1fs total.\n', i, nFam, familyName, toc(tFam));
end

fprintf('\n[c32_stable_ext] ══════════ Phase A complete. ══════════\n');

% ══════════════════════════════════════════════════════════════════════════
%  PHASE B -- pairwise computation (C32 vs each of the 3 stable partners)
% ══════════════════════════════════════════════════════════════════════════
fprintf('\n[c32_stable_ext] ══════════ PHASE B: %d pairs x %d rungs ══════════\n', ...
    numel(PARTNER_FAMILIES), nRung);

results = struct('pairA', {}, 'pairB', {}, 'Ta_nd', {}, 'Ta_days', {}, ...
    'DVproxy_mps', {}, 'DVlb_mps', {}, 'DVpatch_mps', {}, ...
    'TOF_days', {}, 'VoxelId', {}, 'xc_nd', {}, 'yc_nd', {}, 'thc_deg', {});

for r = 1:nRung
    Tmax_r = Ta_multiples_of_pi(r) * pi;
    tRung = tic;
    FHub_r = local_load_rung(HUB_FAMILY, r, BYRUNG_DIR);

    for pIdx = 1:numel(PARTNER_FAMILIES)
        Fp_r = local_load_rung(PARTNER_FAMILIES{pIdx}, r, BYRUNG_DIR);
        m = local_run_pair(FHub_r, Fp_r, grid3, cfg, VU_mps);

        results(end+1) = struct( ...
            'pairA', HUB_FAMILY, 'pairB', PARTNER_FAMILIES{pIdx}, ...
            'Ta_nd', Tmax_r, 'Ta_days', Tmax_r * TU_days, ...
            'DVproxy_mps', m.dv_proxy_mps, 'DVlb_mps', m.dvlb_mps, ...
            'DVpatch_mps', m.dvpatch_mps, 'TOF_days', m.tof_days, ...
            'VoxelId', m.voxelId, 'xc_nd', m.xc, 'yc_nd', m.yc, 'thc_deg', m.thc_deg); %#ok<AGROW>
    end
    clear FHub_r
    fprintf('[c32_stable_ext] Rung %d/%d (Ta=%.2fd) done in %.1fs.\n', r, nRung, Tmax_r * TU_days, toc(tRung));
end

T = struct2table(results);
T = sortrows(T, {'pairB', 'Ta_nd'});

% Flag any pair/rung whose found minimum sits suspiciously close to the
% cap -- that's the signature of the cap clipping the true optimum rather
% than finding it. R52-S in particular has no prior data this cap was
% checked against.
capMps = cfg.fan.DV_cap_nd * VU_mps;
nearCap = T.DVlb_mps > 0.95 * capMps;
if any(nearCap)
    warning(['[c32_stable_ext] %d row(s) have DVlb within 5%% of the ' ...
        'DV_cap_nd=%.3g ceiling (%.1f m/s) -- the cap may be clipping the true ' ...
        'minimum there; consider rerunning with a larger cap.'], ...
        sum(nearCap), cfg.fan.DV_cap_nd, capMps);
end

if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end
outMat = fullfile(RESULTS_DIR, 'ta_asymptote_c32_stable_extended_results.mat');
outCsv = fullfile(RESULTS_DIR, 'ta_asymptote_c32_stable_extended_results.csv');
save(outMat, 'T', 'Ta_multiples_of_pi', 'HUB_FAMILY', 'PARTNER_FAMILIES', 'cfg');
writetable(T, outCsv);

fprintf('\n[c32_stable_ext] ══════════ DONE ══════════\n');
fprintf('  Results table  : %s\n', outMat);
fprintf('  CSV            : %s\n', outCsv);

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
