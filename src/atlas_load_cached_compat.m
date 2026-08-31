function [S, scheme, info] = atlas_load_cached_compat(familyName, cfg)
%ATLAS_LOAD_CACHED_COMPAT  Load a cached family atlas, trying the current
% ('atlas|...') fingerprint scheme first, then falling back to the
% pre-rename ('rs3|...') legacy scheme. NEVER builds from scratch — errors
% if neither scheme finds a valid cache hit, so callers that require a
% cache-only run (no multi-hour rebuild) can rely on this.
%
% Outputs
%   S       atlas struct (as stored under field 'S' in the cache file)
%   scheme  'current' | 'legacy_rs3'  — which fingerprint scheme hit
%   info    cache path/fingerprint info struct from the scheme that hit

[mu, CJ, ~, ~] = cr3bp_family_ic(char(familyName));

[S, hit, info] = atlas_cache_try_load(familyName, mu, CJ, cfg);
if hit
    scheme = 'current';
    return;
end

[S, hit, info] = atlas_cache_try_load_legacy_rs3(familyName, mu, CJ, cfg);
if hit
    scheme = 'legacy_rs3';
    return;
end

info_cur = atlas_cache_get_path(familyName, mu, CJ, cfg);
info_leg = atlas_cache_get_path_legacy_rs3(familyName, mu, CJ, cfg);
error(['atlas_load_cached_compat: no cache hit for family "%s" under ' ...
       'either the current or legacy_rs3 fingerprint scheme. ' ...
       'Checked:\n  current:    %s\n  legacy_rs3: %s'], ...
       familyName, info_cur.fpath, info_leg.fpath);
end
