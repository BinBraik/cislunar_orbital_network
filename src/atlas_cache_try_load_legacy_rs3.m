function [S, hit, info] = atlas_cache_try_load_legacy_rs3(familyName, mu, CJ, cfg)
%ATLAS_CACHE_TRY_LOAD_LEGACY_RS3  Load a cached family atlas built under the
% PRE-RENAME ('rs3|...') fingerprint format, if present and valid. Mirrors
% atlas_cache_try_load.m exactly except for the fingerprint/path source.
%
% Use this ONLY as a fallback after atlas_cache_try_load.m reports a miss —
% it exists to read atlases built before the rs3_* → atlas_* rename
% (see CHANGELOG.md) without needing to rebuild them.

S = [];
hit = false;
info = atlas_cache_get_path_legacy_rs3(familyName, mu, CJ, cfg);

if ~cfg.cache.enable
    return;
end

if exist(info.fpath, 'file') ~= 2
    return;
end

tmp = load(info.fpath, 'S', 'cache_meta');
if ~isfield(tmp, 'S') || ~isfield(tmp, 'cache_meta')
    return;
end

meta = tmp.cache_meta;
if ~isfield(meta, 'fingerprint')
    return;
end
if ~strcmp(meta.fingerprint, info.fingerprint)
    return;
end

S = tmp.S;
if ~isfield(S,'PO_name')
    S.PO_name = S.name;
end
hit = true;
info.cache_meta = meta;
end
