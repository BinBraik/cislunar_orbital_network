function info = atlas_cache_get_path_legacy_rs3(familyName, mu, CJ, cfg)
%ATLAS_CACHE_GET_PATH_LEGACY_RS3  Compute cache filename/path using the
% PRE-RENAME ('rs3|...') fingerprint format. See
% atlas_cache_fingerprint_legacy_rs3.m. Mirrors atlas_cache_get_path.m
% exactly except for the fingerprint source.

if isstring(familyName), familyName = char(familyName); end

fp = atlas_cache_fingerprint_legacy_rs3(familyName, mu, CJ, cfg);
h = compute_md5(fp);
short = h(1:10);

tag = atlas_family_short_tag(familyName);

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
