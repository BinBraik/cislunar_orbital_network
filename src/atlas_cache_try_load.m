function [S, hit, info] = atlas_cache_try_load(familyName, mu, CJ, cfg)
%ATLAS_CACHE_TRY_LOAD  Load cached family atlas if present and valid.

S = [];
hit = false;
info = atlas_cache_get_path(familyName, mu, CJ, cfg);

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
% Validate fingerprint match
if ~isfield(meta, 'fingerprint')
    return;
end
if ~strcmp(meta.fingerprint, info.fingerprint)
    return;
end

S = tmp.S;
% Ensure compatibility fields
if ~isfield(S,'PO_name')
    S.PO_name = S.name;
end
hit = true;
info.cache_meta = meta;
end
