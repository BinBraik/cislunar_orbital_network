function [S, cacheInfo] = rs3_prepare_or_load_family(familyName, cfg, grid3)
%RS3_PREPARE_OR_LOAD_FAMILY  Step 5: load cached atlas or build and cache it.

% Always query family_ic first to get CJ/mu for deterministic cache key.
[mu, CJ, ~, ~] = rs3_core_family_ic(char(familyName));

[S, hit, info] = rs3_cache_try_load_family(familyName, mu, CJ, cfg);
if hit && ~cfg.cache.rebuild
    cacheInfo = struct('hit', true, 'path', info.fpath, 'hash', info.short_hash, 'info', info);
    return;
end

% Build from scratch (Step3 + Step4)
S = rs3_family_prepare_seeds(familyName, cfg, grid3);
S = rs3_family_build_hits(S, cfg, S.grid3);

% Save to cache
info = rs3_cache_save_family(S, cfg);
cacheInfo = struct('hit', false, 'path', info.fpath, 'hash', info.short_hash, 'info', info);
end
