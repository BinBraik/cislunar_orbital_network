function [action, detail] = atlas_ensure_cache_ready( ...
        FAMILIES, TARGET_CFG, OLD_CACHE_DIR, NEW_CACHE_DIR, TOL_REL, FORCE_REBUILD_ALL)
%ATLAS_ENSURE_CACHE_READY  Reconcile OLD_CACHE_DIR against TARGET_CFG and
% make sure NEW_CACHE_DIR ends up holding a valid, current-scheme cached
% atlas for every family in FAMILIES — reusing what already matches
% (no recompute), and rebuilding from scratch (atlas_family_prepare_seeds
% + atlas_family_build_hits, same as atlas_prepare_or_load's cache-miss
% path) whatever doesn't.
%
% For each family:
%   1. Scan OLD_CACHE_DIR for any cached atlas belonging to it (any
%      fingerprint scheme / version_tag / filename convention — matched by
%      VALUE via atlas_cache_fingerprint_parse.m, not filename).
%   2. If a cached atlas's grid/seed/fan/propag/log parameters match
%      TARGET_CFG within TOL_REL -> REUSE it: load as-is (no recompute)
%      and re-save into NEW_CACHE_DIR under the current fingerprint scheme.
%   3. Otherwise -> REBUILD from scratch under TARGET_CFG and save into
%      NEW_CACHE_DIR.
%
% OLD_CACHE_DIR and NEW_CACHE_DIR may be the SAME directory — a common
% call pattern is "make sure this cache dir has everything TARGET_CFG
% needs, reusing whatever's already there": pass the same path for both.
% OLD_CACHE_DIR is only ever read from; safe even when it equals
% NEW_CACHE_DIR because the scan happens once, up front, before any writes.
%
% Inputs
%   FAMILIES           {1×K cell of char}  family names
%   TARGET_CFG          struct              full pipeline config to match/build
%   OLD_CACHE_DIR       char                directory to scan (read-only)
%   NEW_CACHE_DIR       char                directory to ensure is populated
%   TOL_REL             scalar   (optional, default 1e-9) relative tolerance
%                        for comparing parsed fingerprint values
%   FORCE_REBUILD_ALL   logical  (optional, default false) ignore OLD_CACHE_DIR
%                        entirely and rebuild every family from scratch
%
% Outputs
%   action  {K×1 cell of char}  'reused' or 'rebuilt' per family
%   detail  {K×1 cell of char}  source filename for 'reused', '' otherwise

if nargin < 5 || isempty(TOL_REL),           TOL_REL = 1e-9;      end
if nargin < 6 || isempty(FORCE_REBUILD_ALL), FORCE_REBUILD_ALL = false; end

if ~isfolder(OLD_CACHE_DIR)
    error('atlas_ensure_cache_ready: OLD_CACHE_DIR not found: %s', OLD_CACHE_DIR);
end
if ~exist(NEW_CACHE_DIR, 'dir'), mkdir(NEW_CACHE_DIR); end

fprintf('[atlas_ensure_cache_ready] Scanning OLD_CACHE_DIR (read-only):\n  %s\n', OLD_CACHE_DIR);
entries = atlas_cache_inspect(OLD_CACHE_DIR, false);
if isempty(entries)
    fprintf('[atlas_ensure_cache_ready] OLD_CACHE_DIR is empty or has no valid .mat files — every family will be built.\n');
else
    nValid = sum(logical([entries.valid]));
    fprintf('[atlas_ensure_cache_ready] Found %d .mat file(s), %d parsed successfully.\n', numel(entries), nValid);
end

grid3_target = atlas_grid_make(TARGET_CFG);

N      = numel(FAMILIES);
action = repmat({''}, N, 1);
detail = repmat({''}, N, 1);

for k = 1:N
    fam = FAMILIES{k};
    fprintf('\n[atlas_ensure_cache_ready] === %s ===\n', fam);

    [mu_t, CJ_t, ~, ~] = cr3bp_family_ic(fam);

    match_entry = [];
    if ~FORCE_REBUILD_ALL && ~isempty(entries)
        match_entry = local_find_match(entries, fam, mu_t, CJ_t, TARGET_CFG, TOL_REL);
    end

    if ~isempty(match_entry)
        % ── REUSE: load as-is, re-stamp under the current scheme ───────────
        fprintf('[atlas_ensure_cache_ready]   MATCH found: %s (%.1f MB) — reusing, no recompute.\n', ...
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

            fprintf('[atlas_ensure_cache_ready]   Saved (reused) -> %s\n', info.fpath);
            action{k} = 'reused';
            detail{k} = match_entry.fname;
        catch ME
            fprintf('[atlas_ensure_cache_ready]   ERROR reusing match (%s) — falling back to build.\n', ME.message);
            match_entry = [];
        end
    end

    if isempty(match_entry)
        % ── BUILD from scratch under TARGET_CFG ─────────────────────────────
        fprintf('[atlas_ensure_cache_ready]   No matching cached atlas found — building from scratch (this can take a while)...\n');
        cfg_out = TARGET_CFG;
        cfg_out.cache.dir = NEW_CACHE_DIR;

        tBuild = tic;
        S = atlas_family_prepare_seeds(fam, cfg_out, grid3_target);
        S = atlas_family_build_hits(S, cfg_out, S.grid3);
        info = atlas_cache_save(S, cfg_out);
        fprintf('[atlas_ensure_cache_ready]   Built and saved in %.1f s -> %s\n', toc(tBuild), info.fpath);

        action{k} = 'rebuilt';
        detail{k} = '';
    end
end

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('  atlas_ensure_cache_ready SUMMARY — target cache: %s\n', NEW_CACHE_DIR);
fprintf('%s\n', repmat('=', 1, 78));
for k = 1:N
    fprintf('  %-28s  %-8s  %s\n', FAMILIES{k}, action{k}, detail{k});
end
nReused  = sum(strcmp(action, 'reused'));
nRebuilt = sum(strcmp(action, 'rebuilt'));
fprintf('%s\n', repmat('-', 1, 78));
fprintf('  %d reused (no recompute), %d rebuilt from scratch.\n', nReused, nRebuilt);
fprintf('%s\n\n', repmat('=', 1, 78));

end

% ─────────────────────────────────────────────────────────────────────────────
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

% ─────────────────────────────────────────────────────────────────────────────
function ok = local_close(a, b, tol_rel)
if isnan(a) || isnan(b)
    ok = false;
    return;
end
ok = abs(a - b) <= tol_rel * max(1, abs(b));
end
