%% RUN_ATLAS_CACHE_REBUILD_TO_TARGET
% Reconcile an existing (possibly old / mixed-version / unknown-contents)
% atlas cache directory against the TARGET config below (the current
% codebase's paper-matching config), and produce a fresh, clean cache
% directory that the current, unmodified pipeline can use directly.
%
% The actual scan/reuse/rebuild logic lives in src/atlas_ensure_cache_ready.m
% (shared with run_overlap_mintof_maxbudget_preview.m and
% run_overlap_dv_tmax_sweep.m, which call it themselves to self-heal their
% own cache — this script is a thin wrapper for a standalone, cache-only run).
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
%  TARGET CONFIG  — extended-budget atlas (beyond Table 3's nominal 0.2/pi):
%    Rdom=1.2, dx=dy=0.001 (384.4 km), dtheta=1°, ds_seed=0.01,
%    dtheta_fan=0.5°, DV_a=0.3, Ta=3*pi/2, RelTol=AbsTol=v2tol=1e-8.
%  Raising fan.DV_cap_nd / propag.Tmax changes the atlas GENERATION bounds
%  (not just a network-analysis filter), so none of the Table-3 (0.2/pi)
%  cached atlases will match this — every family gets a genuine from-
%  scratch rebuild here, not a reuse.
% ══════════════════════════════════════════════════════════════════════════════
TARGET_CFG = atlas_cfg_defaults();

TARGET_CFG.families.list      = FAMILIES;
TARGET_CFG.families.test_only = false;

TARGET_CFG.grid.dx               = 0.001;
TARGET_CFG.grid.dy               = 0.001;
TARGET_CFG.grid.dtheta           = deg2rad(1);
TARGET_CFG.seed.ds_seed          = 0.01;
TARGET_CFG.propag.Tmax           = 3*pi/2;  % extended budget / base atlas
TARGET_CFG.fan.DV_cap_nd         = 0.3;     % extended budget / base atlas
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
%  RECONCILE  (scan + reuse/rebuild loop lives in src/atlas_ensure_cache_ready.m
%  so run_overlap_mintof_maxbudget_preview.m and run_overlap_dv_tmax_sweep.m
%  can call the exact same logic themselves and self-heal their own cache)
% ══════════════════════════════════════════════════════════════════════════════
atlas_ensure_cache_ready(FAMILIES, TARGET_CFG, OLD_CACHE_DIR, NEW_CACHE_DIR, ...
    TOL_REL, FORCE_REBUILD_ALL);

if ~strcmp(NEW_CACHE_DIR, fullfile(repoRoot, 'atlas_cache'))
    fprintf(['[migrate] NEW_CACHE_DIR is NOT <repoRoot>/atlas_cache — move/rename it there\n' ...
             '  (or update cfg.cache.dir in the runner scripts) before running the preview\n' ...
             '  or sweep scripts.\n']);
end
