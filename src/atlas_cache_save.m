function info = atlas_cache_save(S, cfg)
%ATLAS_CACHE_SAVE  Save minimal family atlas to cache (deterministic key).

if ~cfg.cache.enable
    info = struct();
    return;
end

if ~exist(cfg.cache.dir, 'dir')
    mkdir(cfg.cache.dir);
end

% Ensure compatibility field
if ~isfield(S,'PO_name')
    S.PO_name = S.name;
end

info = atlas_cache_get_path(S.name, S.mu, S.CJ, cfg);

cache_meta = struct();
cache_meta.fingerprint = info.fingerprint;
cache_meta.hash = info.hash;
cache_meta.created = datestr(now, 31);
cache_meta.version_tag = cfg.cache.version_tag;
cache_meta.stats = atlas_cache_stats(S, cfg);

% Strip heavy fields unless explicitly allowed
% Keep a lightweight XY trace of the periodic orbit for story plots.
% (Dense PO is typically stripped for memory; this keeps a small polyline.)
if ~isfield(S,'PO_xy')
    if isfield(S,'Xpo') && ~isempty(S.Xpo)
        stride = 8;
        if isfield(cfg,'diag') && isfield(cfg.diag,'po_stride') && ~isempty(cfg.diag.po_stride)
            stride = max(1, round(cfg.diag.po_stride));
        end
        S.PO_xy = S.Xpo(1:stride:end, 1:2);
        S.PO_stride = stride;
    else
        S.PO_xy = zeros(0,2);
        S.PO_stride = NaN;
    end
end

if isfield(S,'Xpo') && ~cfg.cache.store_dense_po
    S = rmfield(S,'Xpo');
end
if isfield(S,'t_dense') && ~cfg.cache.store_dense_po
    S = rmfield(S,'t_dense');
end

% Keep filename short for Windows path limits; choose save format automatically
use73 = cache_meta.stats.bytes_est > 1.5e9; % ~1.5GB threshold
if use73
    save(info.fpath, 'S', 'cache_meta', '-v7.3');
else
    save(info.fpath, 'S', 'cache_meta', '-v7');
end
end
