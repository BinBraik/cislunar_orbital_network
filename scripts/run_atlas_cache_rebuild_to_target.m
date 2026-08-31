%% RUN_ATLAS_CACHE_REBUILD_TO_TARGET
% Reconcile an existing (possibly old / mixed-version / unknown-contents)
% atlas cache directory against the TARGET config below (the current
% codebase's paper-matching config), and produce a fresh, clean cache
% directory that the current, unmodified pipeline can use directly.
%
% For each of the 13 families:
%   1. Scan OLD_CACHE_DIR for any cached atlas belonging to that family
%      (any fingerprint scheme / version_tag / filename convention — the
%      scan is done by VALUE, via atlas_cache_fingerprint_parse.m, not by
%      filename or literal fingerprint string, so it doesn't matter which
%      pipeline version built it).
%   2. If a cached atlas's grid/seed/fan/propag/log parameters match
%      TARGET_CFG within tolerance -> REUSE it: load the atlas struct as-is
%      (no recomputation) and re-save it into NEW_CACHE_DIR under the
%      CURRENT fingerprint scheme, so plain atlas_prepare_or_load() will
%      find it going forward.
%   3. Otherwise -> REBUILD it from scratch under TARGET_CFG (the same
%      atlas_family_prepare_seeds + atlas_family_build_hits pipeline
%      atlas_prepare_or_load.m uses internally on a cache miss) and save
%      it into NEW_CACHE_DIR.
%
% OLD_CACHE_DIR is only ever READ. NEW_CACHE_DIR is where output is
% written; make it the (emptied) atlas_cache/ you intend to keep using, or
% a scratch directory you move into place afterward.
%
% This script does NOT touch atlas_cache_fingerprint.m or the current
% cache scheme — it only orchestrates reuse-vs-rebuild decisions and calls
% the same build/save functions atlas_prepare_or_load.m already uses.

clear; clc;

% ── repo paths ────────────────────────────────────────────────────────────────
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'src', 'network'));
addpath(fullfile(repoRoot, 'scripts'));
rehash;

% ══════════════════════════════════════════════════════════════════════════════
%  USER KNOBS
% ══════════════════════════════════════════════════════════════════════════════

% Directory to SCAN for reusable atlases (read-only — never written to).
% Point this at your larger/uncertain existing cache.
OLD_CACHE_DIR = '';   % <-- SET THIS, e.g. '/home/braik/.../old_atlas_cache'

% Directory to WRITE the reconciled cache into. Defaults to <repoRoot>/atlas_cache.
% Fine to point this straight at your real atlas_cache/ IF you have already
% emptied it (recommended) — this script will not delete anything itself.
NEW_CACHE_DIR = fullfile(repoRoot, 'atlas_cache');

% Relative tolerance for comparing numeric config fields parsed out of an
% old fingerprint against TARGET_CFG below.
TOL_REL = 1e-9;

% Set true to rebuild every family from scratch regardless of what's found
% in OLD_CACHE_DIR (e.g. to force a fully clean baseline).
FORCE_REBUILD_ALL = false;

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

% ══════════════════════════════════════════════════════════════════════════════
%  TARGET CONFIG  — Table 3 ("Nominal discretization and numerical parameters
%  used in the reachable-set construction") of the paper:
%    Rdom=1.2, dx=dy=0.001 (384.4 km), dtheta=1°, ds_seed=0.01,
%    dtheta_fan=0.5°, DV_a=0.2 (204.6 m/s), Ta=pi (13.66 day),
%    RelTol=AbsTol=v2tol=1e-8.
% ══════════════════════════════════════════════════════════════════════════════
TARGET_CFG = atlas_cfg_defaults();

TARGET_CFG.families.list      = FAMILIES;
TARGET_CFG.families.test_only = false;

TARGET_CFG.grid.dx               = 0.001;
TARGET_CFG.grid.dy               = 0.001;
TARGET_CFG.grid.dtheta           = deg2rad(1);
TARGET_CFG.seed.ds_seed          = 0.01;
TARGET_CFG.propag.Tmax           = pi;      % maximum budget / base atlas
TARGET_CFG.fan.DV_cap_nd         = 0.2;     % maximum budget / base atlas
TARGET_CFG.fan.dtheta_fan        = deg2rad(0.5);
TARGET_CFG.propag.absTol         = 1e-8;
TARGET_CFG.propag.relTol         = 1e-8;
TARGET_CFG.propag.v2tol          = 1e-8;
TARGET_CFG.log.step_len_factor   = 0.75;
TARGET_CFG.log.maxstep_factor    = 2;

TARGET_CFG.cache.enable      = true;
TARGET_CFG.cache.rebuild     = false;
TARGET_CFG.cache.version_tag = 'atlas_v1_keep_masked';  % current-scheme tag

TARGET_CFG.io.save_figs   = false;
TARGET_CFG.io.save_fig    = false;
TARGET_CFG.io.fig_visible = 'off';

if exist('atlas_cfg_validate', 'file') == 2
    atlas_cfg_validate(TARGET_CFG);
end

% ══════════════════════════════════════════════════════════════════════════════
%  VALIDATE PATHS
% ══════════════════════════════════════════════════════════════════════════════
if isempty(OLD_CACHE_DIR)
    error(['run_atlas_cache_rebuild_to_target: set OLD_CACHE_DIR to the ' ...
           'existing cache directory you want to scan for reusable atlases.']);
end
if ~isfolder(OLD_CACHE_DIR)
    error('run_atlas_cache_rebuild_to_target: OLD_CACHE_DIR not found: %s', OLD_CACHE_DIR);
end
if ~exist(NEW_CACHE_DIR, 'dir'), mkdir(NEW_CACHE_DIR); end
if strcmp(OLD_CACHE_DIR, NEW_CACHE_DIR)
    warning(['run_atlas_cache_rebuild_to_target: OLD_CACHE_DIR and NEW_CACHE_DIR ' ...
             'are the same directory. The scan happens once up front so this is ' ...
             'safe, but any REBUILT family will overwrite files in place.']);
end

% ══════════════════════════════════════════════════════════════════════════════
%  SCAN OLD CACHE  (once, read-only)
% ══════════════════════════════════════════════════════════════════════════════
fprintf('[migrate] Scanning OLD_CACHE_DIR (read-only):\n  %s\n', OLD_CACHE_DIR);
entries = atlas_cache_inspect(OLD_CACHE_DIR, false);
if isempty(entries)
    fprintf('[migrate] OLD_CACHE_DIR is empty or has no valid .mat files — every family will be rebuilt.\n');
    entries = [];
else
    nValid = sum(logical([entries.valid]));
    fprintf('[migrate] Found %d .mat file(s), %d parsed successfully.\n', numel(entries), nValid);
end

grid3_target = atlas_grid_make(TARGET_CFG);

% ══════════════════════════════════════════════════════════════════════════════
%  PER-FAMILY: reuse or rebuild
% ══════════════════════════════════════════════════════════════════════════════
N = numel(FAMILIES);
action    = repmat({''}, N, 1);
detail    = repmat({''}, N, 1);

for k = 1:N
    fam = FAMILIES{k};
    fprintf('\n[migrate] === %s ===\n', fam);

    [mu_t, CJ_t, ~, ~] = cr3bp_family_ic(fam);

    match_entry = [];
    if ~FORCE_REBUILD_ALL && ~isempty(entries)
        match_entry = local_find_match(entries, fam, mu_t, CJ_t, TARGET_CFG, TOL_REL);
    end

    if ~isempty(match_entry)
        % ── REUSE: load as-is, re-stamp under the current scheme ───────────
        fprintf('[migrate]   MATCH found in old cache: %s (%.1f MB) — reusing, no recompute.\n', ...
            match_entry.fname, match_entry.MB);
        try
            tmp = load(match_entry.fpath, 'S');
            if ~isfield(tmp, 'S'), error('no S field in file'); end
            S = tmp.S;
            % NOTE: keep S.grid3 exactly as loaded (geometry + Keep mask) —
            % local_find_match already verified dx/dy/dtheta/Rdom/mu/CJ match
            % TARGET_CFG byte-for-byte, so the stored Keep mask is already
            % correct. Do NOT reassign it to atlas_grid_make(TARGET_CFG),
            % which returns geometry only and would silently drop .Keep.

            cfg_out = TARGET_CFG;
            cfg_out.cache.dir = NEW_CACHE_DIR;
            info = atlas_cache_save(S, cfg_out);

            fprintf('[migrate]   Saved (reused) -> %s\n', info.fpath);
            action{k} = 'reused';
            detail{k} = match_entry.fname;
        catch ME
            fprintf('[migrate]   ERROR reusing match (%s) — falling back to rebuild.\n', ME.message);
            match_entry = [];
        end
    end

    if isempty(match_entry)
        % ── REBUILD from scratch under TARGET_CFG ───────────────────────────
        fprintf('[migrate]   No matching cached atlas found — rebuilding from scratch (this can take a while)...\n');
        cfg_out = TARGET_CFG;
        cfg_out.cache.dir = NEW_CACHE_DIR;

        tBuild = tic;
        S = atlas_family_prepare_seeds(fam, cfg_out, grid3_target);
        S = atlas_family_build_hits(S, cfg_out, S.grid3);
        info = atlas_cache_save(S, cfg_out);
        fprintf('[migrate]   Rebuilt and saved in %.1f s -> %s\n', toc(tBuild), info.fpath);

        action{k} = 'rebuilt';
        detail{k} = '';
    end
end

% ══════════════════════════════════════════════════════════════════════════════
%  SUMMARY
% ══════════════════════════════════════════════════════════════════════════════
fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('  SUMMARY — target cache: %s\n', NEW_CACHE_DIR);
fprintf('%s\n', repmat('=', 1, 78));
for k = 1:N
    fprintf('  %-28s  %-8s  %s\n', FAMILIES{k}, action{k}, detail{k});
end
nReused  = sum(strcmp(action, 'reused'));
nRebuilt = sum(strcmp(action, 'rebuilt'));
fprintf('%s\n', repmat('-', 1, 78));
fprintf('  %d reused (no recompute), %d rebuilt from scratch.\n', nReused, nRebuilt);
fprintf('%s\n\n', repmat('=', 1, 78));

if ~strcmp(NEW_CACHE_DIR, fullfile(repoRoot, 'atlas_cache'))
    fprintf(['[migrate] NEW_CACHE_DIR is NOT <repoRoot>/atlas_cache — move/rename it there\n' ...
             '  (or update cfg.cache.dir in the runner scripts) before running the preview\n' ...
             '  or sweep scripts.\n']);
end

% ══════════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════════

function match = local_find_match(entries, fam, mu_t, CJ_t, target_cfg, tol_rel)
%LOCAL_FIND_MATCH  Return the first entry for family FAM whose parsed
% fingerprint values equal TARGET_CFG (and mu/CJ) within tol_rel, or [].
match = [];
for i = 1:numel(entries)
    e = entries(i);
    if ~e.valid || ~strcmp(e.family, fam), continue; end
    p = e.parsed;

    if ~local_close(p.mu, mu_t, tol_rel) || ~local_close(p.CJ, CJ_t, tol_rel)
        continue;
    end

    ok = true;
    ok = ok && local_close(p.grid.R,      target_cfg.grid.Rdom,        tol_rel);
    ok = ok && local_close(p.grid.dx,     target_cfg.grid.dx,          tol_rel);
    ok = ok && local_close(p.grid.dy,     target_cfg.grid.dy,          tol_rel);
    ok = ok && local_close(p.grid.dtheta, target_cfg.grid.dtheta,      tol_rel);
    ok = ok && local_close(p.seed.Tf_scale,  target_cfg.seed.Tf_scale,  tol_rel);
    ok = ok && local_close(p.seed.N_dense,   target_cfg.seed.N_dense,   tol_rel);
    ok = ok && local_close(p.seed.ds_seed,   target_cfg.seed.ds_seed,   tol_rel);
    ok = ok && local_close(p.seed.y_eps,     target_cfg.seed.y_eps,     max(tol_rel, 1e-12));
    ok = ok && local_close(p.seed.minSegPts, target_cfg.seed.minSegPts, tol_rel);
    ok = ok && local_close(p.fan.dtheta_fan, target_cfg.fan.dtheta_fan, tol_rel);
    ok = ok && local_close(p.fan.DV_cap_nd,  target_cfg.fan.DV_cap_nd,  tol_rel);
    ok = ok && local_close(p.propag.Tmax,    target_cfg.propag.Tmax,    tol_rel);
    ok = ok && local_close(p.propag.relTol,  target_cfg.propag.relTol,  tol_rel);
    ok = ok && local_close(p.propag.absTol,  target_cfg.propag.absTol,  tol_rel);
    ok = ok && local_close(p.propag.v2tol,   target_cfg.propag.v2tol,   tol_rel);
    ok = ok && local_close(p.log.segwalk_enable,  double(target_cfg.log.segwalk.enable), max(tol_rel, 1e-12));
    ok = ok && local_close(p.log.segwalk_frac,    target_cfg.log.segwalk.frac,           tol_rel);
    ok = ok && local_close(p.log.step_len_factor, target_cfg.log.step_len_factor,        tol_rel);
    ok = ok && local_close(p.log.maxstep_factor,  target_cfg.log.maxstep_factor,         tol_rel);

    if ok
        match = e;
        return;
    end
end
end

function ok = local_close(a, b, tol_rel)
if isnan(a) || isnan(b)
    ok = false;
    return;
end
ok = abs(a - b) <= tol_rel * max(1, abs(b));
end
