function [S, cacheInfo] = atlas_prepare_or_load(familyName, cfg, grid3)
%RS3_PREPARE_OR_LOAD_FAMILY  Step 5: load cached atlas or build and cache it.

% Always query family_ic first to get CJ/mu for deterministic cache key.
[mu, CJ, ~, ~] = cr3bp_family_ic(char(familyName));

[S, hit, info] = atlas_cache_try_load(familyName, mu, CJ, cfg);
if hit && ~cfg.cache.rebuild
    cacheInfo = struct('hit', true, 'path', info.fpath, 'hash', info.short_hash, 'info', info);
    return;
end

% Build from scratch (Step3 + Step4)
S = atlas_family_prepare_seeds(familyName, cfg, grid3);
S = atlas_family_build_hits(S, cfg, S.grid3);

% Save to cache
info = atlas_cache_save(S, cfg);
cacheInfo = struct('hit', false, 'path', info.fpath, 'hash', info.short_hash, 'info', info);
end
