function run_ta_asymptote_c32_stable_extended_family(familyName, nWorkers)
%RUN_TA_ASYMPTOTE_C32_STABLE_EXTENDED_FAMILY  Build the extended-Ta,
% reduced-cap atlas for ONE family and save each rung's footprint to its
% own file immediately as it's computed -- no monolithic per-family
% bundle stage. Unlike run_ta_asymptote_sweep_full.m (which builds a
% single Fcell for all rungs, saves it once, and needs a later repack
% pass to split it into per-rung files), this script IS the per-rung
% cache: nothing to repack afterwards.
%
% WHY A SEPARATE, REDUCED CAP: this test doubles Tmax to 32pi (~437
% days) to chase the stable-pair asymptote further than the original
% 16pi/218-day run got. DV_cap_nd=0.2 (204.6 m/s) was never the binding
% constraint on the C32-hub run's stable pairs -- both C32-R21S and
% C32-R31S landed well under it (130.5 and 91.8 m/s at Ta=218.6d, still
% slowly decreasing, ~5.4 m/s/rung). DV_cap_nd=0.15 (153.5 m/s) keeps
% ~18%+ margin above the harder (R21S) curve's last known value while
% cutting the fan search space -- and therefore build cost/memory -- by
% 25% relative to 0.2, on top of the already much longer Tmax.
%
% MEANT TO RUN CONCURRENTLY: launched as one of three simultaneous
% `matlab -batch` processes (one per family: Cycler 32, Resonant 2to1
% Stable, Resonant 3to1 Stable) sharing one big-memory node -- see
% run_ta_asymptote_c32_stable_extended.slurm. Each instance gets its own
% parpool of nWorkers, sized by the launcher so all three pools fit the
% node's core budget together. The three families are otherwise fully
% independent (different orbits, no shared state), so this is safe --
% unlike parallelizing the derive-subset/footprint loop WITHIN a family,
% which would mean broadcasting one (multi-GB-to-tens-of-GB) fat atlas to
% every worker and multiplying memory by worker count -- exactly the bug
% that caused the OOM on run_ta_asymptote_sweep_full.m's idle 60-worker
% pool. Parallelizing across families instead avoids that entirely: each
% process holds only its own atlas.
%
% Idempotent per rung: rerunning skips any rung whose file already
% exists and matches this exact ladder/cap.

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

Ta_multiples_of_pi = sort(2.^linspace(-1.5, 5, 24), 'descend');   % 0.354pi ... 32pi, 24 rungs
DV_cap_nd = 0.15;

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

cfg.par.enable = (nWorkers > 0);

cfg.io.save_figs   = false;
cfg.io.save_fig    = false;
cfg.io.fig_visible = 'off';

if exist('atlas_cfg_validate', 'file') == 2
    atlas_cfg_validate(cfg);
end

VU_mps  = cfg.units.VU_mps;
TU_days = cfg.units.TU_days;

grid3 = atlas_grid_make(cfg);
nRung = numel(Ta_multiples_of_pi);

OUTPUT_DIR = fullfile(repoRoot, 'ta_asymptote_c32_stable_extended_results');
BYRUNG_DIR = fullfile(OUTPUT_DIR, 'footprints_by_rung');
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
    fprintf('[c32_stable_ext:%s] all %d rungs already cached at this ladder/cap -- nothing to do.\n', ...
        familyName, nRung);
    return;
end

if nWorkers > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', nWorkers);
        fprintf('[c32_stable_ext:%s] Started parpool with %d workers.\n', familyName, nWorkers);
    end
end

tFam = tic;
Tmax_fat = Ta_multiples_of_pi(1) * pi;
fprintf('[c32_stable_ext:%s] building fat atlas at Ta=%.4gpi (%.2f days), DV_cap_nd=%.3g...\n', ...
    familyName, Ta_multiples_of_pi(1), Tmax_fat * TU_days, DV_cap_nd);

cfg_fat = cfg;
cfg_fat.propag.Tmax = Tmax_fat;
[S_fat, cacheInfo] = atlas_prepare_or_load(familyName, cfg_fat, grid3);
fprintf('[c32_stable_ext:%s]   fat atlas %s in %.1fs.\n', familyName, local_hitstr(cacheInfo), toc(tFam));

for r = 1:nRung
    rpath = fullfile(famDir, sprintf('rung%02d.mat', r));
    if isfile(rpath)
        ck = load(rpath, 'Ta_multiples_of_pi', 'DV_cap_nd');
        if isequal(ck.Ta_multiples_of_pi, Ta_multiples_of_pi) && ck.DV_cap_nd == DV_cap_nd
            fprintf('[c32_stable_ext:%s] rung %d/%d already cached, skipping.\n', familyName, r, nRung);
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
    save(rpath, 'F', 'Ta_multiples_of_pi', 'DV_cap_nd', '-v7.3');
    clear F
    fprintf('[c32_stable_ext:%s] rung %d/%d (Ta=%.2fd) done and saved (%.1fs elapsed).\n', ...
        familyName, r, nRung, Tmax_r * TU_days, toc(tFam));
end

clear S_fat   % discard the full atlas -- every rung's footprint is already on disk
fprintf('[c32_stable_ext:%s] ALL %d rungs done in %.1fs total.\n', familyName, nRung, toc(tFam));

end

% ─────────────────────────────────────────────────────────────────────────────
%  LOCAL FUNCTIONS  (verbatim from run_ta_asymptote_sweep_full.m)
% ─────────────────────────────────────────────────────────────────────────────

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
