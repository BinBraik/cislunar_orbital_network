function [S, hit, info] = rs3_cache_try_load_family(familyName, mu, CJ, cfg)
%RS3_CACHE_TRY_LOAD_FAMILY  Load cached family atlas if present and valid.

S = [];
hit = false;
info = rs3_cache_get_path(familyName, mu, CJ, cfg);

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
