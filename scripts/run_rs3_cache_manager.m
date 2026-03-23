%% RUN_RS3_CACHE_MANAGER  —  inspect, review, and derive atlas caches
%
% Cache-housekeeping utility.  Supports three actions:
%
%   'inspect'     — print the full cache table (family, config key, size, date).
%                   Always safe to run; read-only.
%   'show_config' — pretty-print the exact MATLAB config block for one or more
%                   cache entries so you can verify their grid/fan/propag settings.
%   'derive'      — build reduced-parameter subset caches from an existing source
%                   (e.g. shorter Tmax, tighter DV_cap) without re-integrating.
%
% Typical workflow:
%   1. Run with ACTION = 'inspect' → read the table, note index numbers.
%   2. Set ACTION = 'show_config', SHOW_IDX = [n] → verify parameters.
%   3. Set ACTION = 'derive', SOURCE_IDX = n, edit SUBSET_OVERRIDES with the
%      parameters you want to change.  The script finds all families sharing
%      the same config_key as entry n and derives a subset for each.
%   4. The printed cfg block at the end can be pasted into any batch runner.
%
% Prerequisites:
%   - At least one cached atlas must exist in rs3_cache/.
%     Build them with run_rs3_one_family_atlas_and_plots.m.
%
% User knobs:
%   ACTION           — 'inspect' | 'show_config' | 'derive'
%   CACHE_DIR        — path to cache folder; '' = default (repo/rs3_cache)
%   SHOW_IDX         — [integers] — entry indices to display (show_config only)
%   SOURCE_IDX       — integer — representative entry to derive from
%   FAMILIES         — {} = all families with matching config; or name list
%   SUBSET_OVERRIDES — struct with fields to change (propag.Tmax, fan.DV_cap_nd,
%                      grid.dx/dy/dtheta); absent fields are inherited unchanged
%   FORCE_REBUILD    — true to overwrite an existing derived cache
%
% Outputs (derive action only, written to rs3_cache/):
%   <MD5>_<family>_rs3_v2_keep_masked.mat  — derived subset cache file
%   A copy-pasteable cfg block is printed to the console.

clear; clc;
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot,'src'));
addpath(fullfile(repoRoot,'scripts'));
rehash;

% =====================================================================
% ========================= USER KNOBS ================================
% =====================================================================

% ── ACTION ──────────────────────────────────────────────────────────
%   'inspect'     — print the full cache table (read-only, always safe)
%   'show_config' — print the complete MATLAB config block for one or
%                   more entries so you can verify their settings
%   'derive'      — build subset caches from a selected source
ACTION = 'inspect';

% ── CACHE DIRECTORY ─────────────────────────────────────────────────
CACHE_DIR = '';   % '' = default (repo/rs3_cache)

% ── show_config ─────────────────────────────────────────────────────
% Indices of entries to show (from the inspect table, first column).
SHOW_IDX = [1];

% ── derive ──────────────────────────────────────────────────────────

% ONE representative index from the inspect table.
% The script finds this entry's config_key, then processes every family
% in the cache that shares the same key (same dx/dy/dtheta/tolerances/…).
SOURCE_IDX = 1;

% Restrict to specific families (leave empty = all with matching config).
FAMILIES = {};   % e.g. {'Lyapunov L1', 'Lyapunov L2'}

% Parameters to CHANGE from the source.
% Leave a field absent to inherit the source value unchanged.
% Supported: propag.Tmax, fan.DV_cap_nd, grid.dx, grid.dy, grid.dtheta
% (grid coarsening: sub must be >= source; cannot refine)
SUBSET_OVERRIDES = struct();
SUBSET_OVERRIDES.propag.Tmax   = pi/2;     % shorten integration time
SUBSET_OVERRIDES.fan.DV_cap_nd = 0.10;     % tighter DV budget
% SUBSET_OVERRIDES.grid.dtheta = deg2rad(4); % coarser theta bins

% Overwrite an existing derived cache if it already exists?
FORCE_REBUILD = false;

% =====================================================================
% =====================================================================

if isempty(CACHE_DIR)
    CACHE_DIR = fullfile(rs3_repo_root(), 'rs3_cache');
end
if ~exist(CACHE_DIR, 'dir')
    error('Cache directory not found: %s', CACHE_DIR);
end

% ── always inspect first ────────────────────────────────────────────
fprintf('[cache_mgr] Scanning: %s\n', CACHE_DIR);
entries = rs3_cache_inspect(CACHE_DIR, true);

if isempty(entries)
    fprintf('[cache_mgr] Cache is empty — nothing to do.\n');
    return;
end

valid_mask = logical([entries.valid]);

% ── ACTION: inspect ─────────────────────────────────────────────────
if strcmpi(ACTION, 'inspect')
    fprintf('[cache_mgr] ACTION=inspect complete.\n');
    fprintf('  → Set ACTION=''show_config'' with SHOW_IDX=[n] to see a full config block.\n');
    fprintf('  → Set ACTION=''derive'' with SOURCE_IDX=n and SUBSET_OVERRIDES to create a subset.\n\n');
    return;
end

% ── ACTION: show_config ─────────────────────────────────────────────
if strcmpi(ACTION, 'show_config')
    for k = 1:numel(SHOW_IDX)
        idx = SHOW_IDX(k);
        if idx < 1 || idx > numel(entries)
            fprintf('[cache_mgr] SHOW_IDX %d out of range (1..%d)\n', idx, numel(entries));
            continue;
        end
        e = entries(idx);
        if ~e.valid
            fprintf('[cache_mgr] Entry %d is invalid: %s\n', idx, e.error);
            continue;
        end
        local_print_config_block(e, entries, idx);
    end
    return;
end

% ── ACTION: derive ──────────────────────────────────────────────────
if ~strcmpi(ACTION, 'derive')
    error('[cache_mgr] Unknown ACTION ''%s''. Use ''inspect'', ''show_config'', or ''derive''.', ACTION);
end

% -- validate SOURCE_IDX
if SOURCE_IDX < 1 || SOURCE_IDX > numel(entries)
    error('[cache_mgr] SOURCE_IDX=%d out of range (1..%d).', SOURCE_IDX, numel(entries));
end
src_rep = entries(SOURCE_IDX);
if ~src_rep.valid
    error('[cache_mgr] Entry SOURCE_IDX=%d is invalid: %s', SOURCE_IDX, src_rep.error);
end

fprintf('[cache_mgr] Source representative: #%d  %s  (config_key=%s)\n', ...
    SOURCE_IDX, src_rep.family, src_rep.config_key);
local_print_config_block(src_rep, entries, SOURCE_IDX);

% -- find all entries sharing this config_key
same_key = valid_mask & strcmp({entries.config_key}, src_rep.config_key);
candidate_entries = entries(same_key);

if ~isempty(FAMILIES)
    filt = false(size(candidate_entries));
    for k = 1:numel(candidate_entries)
        filt(k) = any(strcmpi(candidate_entries(k).family, FAMILIES));
    end
    candidate_entries = candidate_entries(filt);
end

if isempty(candidate_entries)
    fprintf('[cache_mgr] No matching entries found for config_key=%s.\n', src_rep.config_key);
    return;
end

fprintf('[cache_mgr] Found %d entries with config_key=%s:\n', numel(candidate_entries), src_rep.config_key);
for k = 1:numel(candidate_entries)
    fprintf('    #%-3d  %s\n', candidate_entries(k).idx, candidate_entries(k).family);
end
fprintf('\n');

% -- build cfg_base from defaults
cfg_base = rs3_cfg_defaults();

nDerived = 0;
nSkipped = 0;

for ei = 1:numel(candidate_entries)
    e_src = candidate_entries(ei);
    fam   = e_src.family;

    fprintf('[cache_mgr] === %s ===\n', fam);

    % Reconstruct source cfg from parsed fingerprint
    cfg_src = local_parsed_to_cfg(cfg_base, e_src.parsed);

    % Build sub-config by applying overrides
    cfg_sub = local_deep_merge(cfg_src, SUBSET_OVERRIDES);

    % Check if this is actually a change
    sub_fp  = rs3_cache_fingerprint_family(fam, e_src.parsed.mu, e_src.parsed.CJ, cfg_sub);
    if strcmp(sub_fp, e_src.parsed.raw)
        fprintf('[cache_mgr]   SKIP — SUBSET_OVERRIDES produce identical config (nothing to derive)\n\n');
        nSkipped = nSkipped + 1;
        continue;
    end

    % Check if derived cache already exists
    info_sub = rs3_cache_get_path(fam, e_src.parsed.mu, e_src.parsed.CJ, cfg_sub);
    if exist(info_sub.fpath, 'file') == 2 && ~FORCE_REBUILD
        fprintf('[cache_mgr]   SKIP — derived cache already exists: %s\n\n', info_sub.fname);
        nSkipped = nSkipped + 1;
        continue;
    end

    % Load full atlas
    fprintf('[cache_mgr]   Loading source (%.1f MB): %s\n', e_src.MB, e_src.fname);
    try
        tmp = load(e_src.fpath, 'S');
    catch ME
        fprintf('[cache_mgr]   ERROR loading: %s\n\n', ME.message);
        nSkipped = nSkipped + 1;
        continue;
    end
    if ~isfield(tmp, 'S')
        fprintf('[cache_mgr]   ERROR: no S field in file\n\n');
        nSkipped = nSkipped + 1;
        continue;
    end
    S_src = tmp.S;
    S_src.grid3 = rs3_grid_make(cfg_src);   % rebuild function handles

    % Derive subset
    try
        S_sub = rs3_atlas_derive_subset(S_src, cfg_sub);
    catch ME
        fprintf('[cache_mgr]   ERROR during derive: %s\n\n', ME.message);
        nSkipped = nSkipped + 1;
        continue;
    end

    % Save
    info = rs3_cache_save_family(S_sub, cfg_sub);
    fprintf('[cache_mgr]   Saved: %s\n\n', info.fpath);
    nDerived = nDerived + 1;
end

fprintf('[cache_mgr] Done — %d derived, %d skipped.\n\n', nDerived, nSkipped);

if nDerived > 0
    % Print the full config block to paste into the batch runner
    fprintf('%s\n', repmat('─', 1, 72));
    fprintf('  PASTE THIS INTO run_rs4_all_pairs_summary.m / run_rs4_overlap_and_visuals.m:\n');
    fprintf('%s\n', repmat('─', 1, 72));

    % Build a temporary parsed struct for the sub-config using representative
    p = rs3_cache_fingerprint_parse( ...
            rs3_cache_fingerprint_family(src_rep.family, src_rep.parsed.mu, ...
                src_rep.parsed.CJ, local_deep_merge(local_parsed_to_cfg(cfg_base, src_rep.parsed), SUBSET_OVERRIDES)));
    local_print_matlab_block(p);
    fprintf('%s\n\n', repmat('─', 1, 72));
end

% =====================================================================
%  Local helpers
% =====================================================================

function local_print_config_block(e, entries, idx)
% Print a full MATLAB config block for entry e.
p = e.parsed;

% find all families sharing the same config_key
same = strcmp({entries.config_key}, e.config_key) & logical([entries.valid]);
fams = {entries(same).family};

fprintf('\n%s\n', repmat('─', 1, 72));
fprintf('  Entry #%d — %s\n', idx, e.family);
fprintf('  config_key: %s  (%d entries share this config: %s)\n', ...
    e.config_key, sum(same), strjoin(fams, ', '));
fprintf('  File: %s  (%.1f MB)\n', e.fname, e.MB);
fprintf('  Created: %s\n', e.created);
fprintf('%s\n', repmat('─', 1, 72));
local_print_matlab_block(p);
fprintf('%s\n\n', repmat('─', 1, 72));
end

function local_print_matlab_block(p)
% Print a copy-pasteable MATLAB config block from a parsed fingerprint.
fprintf('  %% grid settings (must match for both atlases)\n');
fprintf('  cfg.grid.dx     = %.6g;\n',         p.grid.dx);
fprintf('  cfg.grid.dy     = %.6g;\n',         p.grid.dy);
fprintf('  cfg.grid.dtheta = deg2rad(%.6g);  %% = %.4f°\n', ...
    rad2deg(p.grid.dtheta), rad2deg(p.grid.dtheta));
fprintf('  cfg.seed.ds_seed   = %.6g;\n',      p.seed.ds_seed);
fprintf('  %% propagation/fan\n');
fprintf('  cfg.propag.Tmax    = %s;  %% = %s\n', local_tmax_matlab_expr(p.propag.Tmax), local_tmax_str(p.propag.Tmax));
fprintf('  cfg.fan.DV_cap_nd  = %.6g;\n',      p.fan.DV_cap_nd);
fprintf('  cfg.fan.dtheta_fan = deg2rad(%.6g);  %% = %.4f°\n', ...
    rad2deg(p.fan.dtheta_fan), rad2deg(p.fan.dtheta_fan));
fprintf('  cfg.propag.absTol  = %.6g;\n',      p.propag.absTol);
fprintf('  cfg.propag.relTol  = %.6g;\n',      p.propag.relTol);
fprintf('  cfg.propag.v2tol   = %.6g;\n',      p.propag.v2tol);
fprintf('  cfg.log.step_len_factor = %.6g;\n', p.log.step_len_factor);
fprintf('  cfg.log.maxstep_factor  = %.6g;\n', p.log.maxstep_factor);
end

function s = local_tmax_str(T)
if isnan(T), s = '?'; return; end
frac = T/pi;
if     abs(frac-1)   < 1e-6, s = 'π';
elseif abs(frac-0.5) < 1e-6, s = 'π/2';
elseif abs(frac-2)   < 1e-6, s = '2π';
elseif abs(frac-0.25)< 1e-6, s = 'π/4';
else,                          s = sprintf('%.4gπ', frac);
end
end

function expr = local_tmax_matlab_expr(T)
% Return a pasteable MATLAB expression for Tmax — exact pi fractions where
% possible, otherwise %.16g so the fingerprint round-trips correctly.
if isnan(T), expr = 'NaN'; return; end
frac = T/pi;
if     abs(frac-1)   < 1e-6, expr = 'pi';
elseif abs(frac-0.5) < 1e-6, expr = 'pi/2';
elseif abs(frac-2)   < 1e-6, expr = '2*pi';
elseif abs(frac-0.25)< 1e-6, expr = 'pi/4';
else,                          expr = sprintf('%.16g', T);
end
end

function cfg = local_parsed_to_cfg(cfg_base, p)
% Reconstruct a full cfg from a parsed fingerprint + cfg_base defaults.
cfg = cfg_base;
if ~isstruct(p), return; end
if isfield(p,'grid')
    if isfield(p.grid,'R'),      cfg.grid.Rdom    = p.grid.R;      end
    if isfield(p.grid,'dx'),     cfg.grid.dx      = p.grid.dx;     end
    if isfield(p.grid,'dy'),     cfg.grid.dy      = p.grid.dy;     end
    if isfield(p.grid,'dtheta'), cfg.grid.dtheta  = p.grid.dtheta; end
end
if isfield(p,'seed')
    if isfield(p.seed,'Tf_scale'),  cfg.seed.Tf_scale  = p.seed.Tf_scale;        end
    if isfield(p.seed,'N_dense'),   cfg.seed.N_dense   = round(p.seed.N_dense);  end
    if isfield(p.seed,'ds_seed'),   cfg.seed.ds_seed   = p.seed.ds_seed;         end
    if isfield(p.seed,'y_eps'),     cfg.seed.y_eps     = p.seed.y_eps;           end
    if isfield(p.seed,'minSegPts'), cfg.seed.minSegPts = round(p.seed.minSegPts);end
end
if isfield(p,'fan')
    if isfield(p.fan,'dtheta_fan'), cfg.fan.dtheta_fan = p.fan.dtheta_fan; end
    if isfield(p.fan,'DV_cap_nd'),  cfg.fan.DV_cap_nd  = p.fan.DV_cap_nd;  end
end
if isfield(p,'propag')
    if isfield(p.propag,'Tmax'),   cfg.propag.Tmax   = p.propag.Tmax;   end
    if isfield(p.propag,'relTol'), cfg.propag.relTol = p.propag.relTol; end
    if isfield(p.propag,'absTol'), cfg.propag.absTol = p.propag.absTol; end
    if isfield(p.propag,'v2tol'),  cfg.propag.v2tol  = p.propag.v2tol;  end
end
if isfield(p,'log')
    if isfield(p.log,'segwalk_enable'),  cfg.log.segwalk.enable    = logical(p.log.segwalk_enable); end
    if isfield(p.log,'segwalk_frac'),    cfg.log.segwalk.frac      = p.log.segwalk_frac;  end
    if isfield(p.log,'step_len_factor'), cfg.log.step_len_factor   = p.log.step_len_factor; end
    if isfield(p.log,'maxstep_factor'),  cfg.log.maxstep_factor    = p.log.maxstep_factor; end
end
if isfield(p,'version_tag') && ~isempty(p.version_tag)
    cfg.cache.version_tag = p.version_tag;
end
end

function out = local_deep_merge(base, over)
% Recursively override fields of base struct with fields from over struct.
out = base;
if ~isstruct(over), return; end
fields = fieldnames(over);
for i = 1:numel(fields)
    f = fields{i};
    if isstruct(over.(f)) && isfield(out, f) && isstruct(out.(f))
        out.(f) = local_deep_merge(out.(f), over.(f));
    else
        out.(f) = over.(f);
    end
end
end
