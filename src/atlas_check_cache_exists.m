function [all_present, report] = atlas_check_cache_exists(families, cfg)
%ATLAS_CHECK_CACHE_EXISTS  Pre-flight check: verify a cached atlas file exists
% for every family in FAMILIES under the current cfg fingerprint, without
% loading any of them.
%
% This mirrors the fingerprint/path logic used by atlas_prepare_or_load /
% atlas_cache_try_load, so a "present" result here guarantees
% atlas_prepare_or_load will hit the cache rather than rebuild from scratch.
%
% Inputs
%   families  {1×K cell of char}  family names (see cr3bp_family_ic)
%   cfg       struct              full pipeline config (cfg.cache.* used)
%
% Outputs
%   all_present  logical          true iff every family has a cache hit
%   report       table            one row per family:
%                                    name, present (logical), fpath (string)

K = numel(families);
name    = cell(K, 1);
present = false(K, 1);
fpath   = cell(K, 1);

for k = 1:K
    fam = families{k};
    [mu, CJ, ~, ~] = cr3bp_family_ic(char(fam));
    info = atlas_cache_get_path(fam, mu, CJ, cfg);

    name{k}    = fam;
    fpath{k}   = info.fpath;
    present(k) = cfg.cache.enable && (exist(info.fpath, 'file') == 2);
end

report = table(name(:), present(:), fpath(:), ...
    'VariableNames', {'family', 'present', 'fpath'});

all_present = all(present);

if ~all_present
    missing = report.family(~report.present);
    fprintf('[atlas_check_cache_exists] Missing cached atlas for %d/%d families:\n', ...
        numel(missing), K);
    for k = 1:numel(missing)
        fprintf('    - %s\n', missing{k});
    end
    fprintf('  (cfg.cache.enable=%d, cfg.cache.dir=%s)\n', ...
        cfg.cache.enable, cfg.cache.dir);
else
    fprintf('[atlas_check_cache_exists] All %d family atlases found in cache:\n  %s\n', ...
        K, cfg.cache.dir);
end

end
