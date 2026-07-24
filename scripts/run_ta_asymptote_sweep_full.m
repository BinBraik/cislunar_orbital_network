%% RUN_TA_ASYMPTOTE_SWEEP_FULL
% Atlas-building (Phase A) half of the full 13-family / 78-pair extension
% of run_ta_asymptote_sweep.m. Builds and caches every family's per-rung
% footprint; does NOT compute any pairwise overlaps itself -- that is
% run_ta_asymptote_phaseB_safe.m's job, run as a separate step/job once
% however many families you need are cached here.
%
% Same physics/config as the single-pair run: DV_cap_nd = 0.2 (paper
% nominal) held FIXED, Table 3 grid/discretization unchanged, Ta ladder
% built largest-first with every smaller rung derived from that one fat
% atlas via atlas_derive_subset (shrink-only, never expansion).
%
% HISTORY / WHY PHASE B WAS REMOVED FROM THIS SCRIPT: an earlier revision
% also computed all 78 pairs here, keeping every family's footprint
% resident in memory across the whole run on the assumption (from the
% original codebase's own comments, calibrated at the paper's Ta=pi) that
% footprints were ~5-25MB each. At Ta=16pi they turned out to be 4-40GB
% EACH (measured directly on disk) -- trajectories touch far more
% distinct voxels over 218 days. Accumulating 10+ of those pushed a 900GB
% allocation to 94%+ full while STILL in the middle of Phase A, before
% ever reaching the also-unsafe parfor-broadcast Phase B. Now: build,
% save to disk, discard from memory, never accumulate. Nothing in this
% script needs a family's footprint again once it's saved.
%
% Family-by-family: build the fat atlas for family i, immediately derive
% every rung's footprint from it, save to disk, discard the full atlas
% AND the footprint before moving to family i+1. Checkpointed per family
% via the footprint cache file itself: a family whose file already exists
% (and matches the current Ta ladder / DV_cap) is skipped entirely on
% resubmission, no re-integration.
%
% PAIRS is all 78 unordered pairs among the 13 families -- computed here
% only as a count for the log message; the actual pairwise work happens
% in run_ta_asymptote_phaseB_safe.m, which streams one rung's worth of
% footprints at a time regardless of how large the cached files are.
%
% Ta ladder is the 12-rung ladder (2.^((-3:8)/2)*pi, i.e. 0.354pi ...
% 16pi) -- the denser 25-rung ladder was for the single DPO->R31-U
% deep-dive and isn't needed for a first full-network pass.

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

% Informational only -- pairwise computation itself lives in the separate
% run_ta_asymptote_phaseB_safe.m (this script only builds/caches atlases).
nPairs = nFam * (nFam - 1) / 2;

fprintf('[full_sweep] Families: %d   Pairs (for later Phase B): %d   Ta rungs: %d\n', nFam, nPairs, nRung);
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
fprintf('\n[full_sweep] ══════════ PHASE A: per-family footprints (build/save/discard) ══════════\n');

% NOTE: this loop does NOT keep a persistent Fp(family, rung) cache across
% iterations. Earlier revisions did (`Fp(i,:) = Fcell` retained for the
% whole script), which turned out to be a serious bug once real footprint
% sizes were measured: at Ta=16pi they run 4-40GB EACH (not the ~5-25MB
% the original codebase assumed at Ta=pi), so accumulating 10+ of them in
% memory pushed a 900GB allocation to 94%+ full while still only in the
% MIDDLE of Phase A -- before ever reaching the (also-unsafe) old Phase B.
% Each family's footprint is saved to disk and then dropped from memory
% immediately; nothing here needs it again. Pairwise computation is a
% separate, deliberately disk-streaming step -- see
% run_ta_asymptote_phaseB_safe.m, which loads one rung's worth of
% footprints at a time instead of holding everything at once.

for i = 1:nFam
    key   = local_fieldkey(FAMILIES{i});
    fpath = fullfile(FOOTPRINT_DIR, [key '.mat']);

    if isfile(fpath)
        ck = load(fpath, 'Ta_multiples_of_pi', 'DV_cap_nd');   % header only, not Fcell
        if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) && ...
                ck.DV_cap_nd == cfg.fan.DV_cap_nd
            fprintf('[full_sweep] (%d/%d) %-25s -- footprints cached, skipping.\n', ...
                i, nFam, FAMILIES{i});
            clear ck
            continue;
        else
            fprintf('[full_sweep] (%d/%d) %-25s -- cached footprints stale (config changed), rebuilding.\n', ...
                i, nFam, FAMILIES{i});
        end
        clear ck
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
        fprintf('[full_sweep]   %-25s rung %d/%d footprint done (%.1fs elapsed).\n', ...
            FAMILIES{i}, r, nRung, toc(tFam));
    end
    clear S_fat_i   % discard the full atlas -- only the footprints survive

    DV_cap_nd = cfg.fan.DV_cap_nd; %#ok<NASGU> -- saved for cache validation
    save(fpath, 'Fcell', 'Ta_multiples_of_pi', 'DV_cap_nd', '-v7.3');
    clear Fcell   % do NOT accumulate -- saved to disk is enough, Phase B reads it fresh

    fprintf('[full_sweep] (%d/%d) %-25s -- done in %.1fs total.\n', ...
        i, nFam, FAMILIES{i}, toc(tFam));
end

fprintf('\n[full_sweep] ══════════ Phase A complete. All %d families'' footprints cached to:\n', nFam);
fprintf('  %s\n', FOOTPRINT_DIR);
fprintf('Run scripts/run_ta_asymptote_phaseB_safe.m to compute the 78 pairwise curves\n');
fprintf('from this cache (safe regardless of footprint file size -- streams one rung at a time).\n');

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

