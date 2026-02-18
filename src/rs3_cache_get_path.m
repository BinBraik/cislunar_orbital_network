function info = rs3_cache_get_path(familyName, mu, CJ, cfg)
%RS3_CACHE_GET_PATH  Compute short cache filename + metadata.

if isstring(familyName), familyName = char(familyName); end

fp = rs3_cache_fingerprint_family(familyName, mu, CJ, cfg);
h = rs3_md5(fp);
short = h(1:10);

tag = rs3_family_short_tag(familyName);

% Keep filename short for Windows path limits
fname = sprintf('%s_%s_%s.mat', cfg.cache.version_tag, tag, short);

info = struct();
info.familyName = familyName;
info.mu = mu;
info.CJ = CJ;
info.fingerprint = fp;
info.hash = h;
info.short_hash = short;
info.fname = fname;
info.fpath = fullfile(cfg.cache.dir, fname);
end
