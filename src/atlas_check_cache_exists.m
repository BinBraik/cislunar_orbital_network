function [all_present, report] = atlas_check_cache_exists(families, cfg)
%ATLAS_CHECK_CACHE_EXISTS  Pre-flight check: verify a cached atlas file exists
% for every family in FAMILIES, without loading any of them.
%
% Checks TWO fingerprint schemes per family, in order:
%   1. current    — atlas_cache_get_path.m ('atlas|...' fingerprint prefix)
%   2. legacy_rs3 — atlas_cache_get_path_legacy_rs3.m ('rs3|...' prefix,
%                   for atlases built before the rs3_* -> atlas_* rename,
%                   see CHANGELOG.md)
% A "present" result under EITHER scheme guarantees
% atlas_load_cached_compat.m will hit the cache rather than rebuild.
%
% Inputs
%   families  {1×K cell of char}  family names (see cr3bp_family_ic)
%   cfg       struct              full pipeline config (cfg.cache.* used)
%
% Outputs
%   all_present  logical          true iff every family has a cache hit
%                                  under at least one scheme
%   report       table            one row per family:
%                                    name, present (logical), scheme
%                                    ('current'/'legacy_rs3'/'none'),
%                                    fpath_current, fpath_legacy

K = numel(families);
name          = cell(K, 1);
present       = false(K, 1);
scheme        = repmat({'none'}, K, 1);
fpath_current = cell(K, 1);
fpath_legacy  = cell(K, 1);

for k = 1:K
    fam = families{k};
    [mu, CJ, ~, ~] = cr3bp_family_ic(char(fam));

    info_cur = atlas_cache_get_path(fam, mu, CJ, cfg);
    info_leg = atlas_cache_get_path_legacy_rs3(fam, mu, CJ, cfg);

    name{k}          = fam;
    fpath_current{k} = info_cur.fpath;
    fpath_legacy{k}  = info_leg.fpath;

    if cfg.cache.enable && (exist(info_cur.fpath, 'file') == 2)
        present(k) = true;  scheme{k} = 'current';
    elseif cfg.cache.enable && (exist(info_leg.fpath, 'file') == 2)
        present(k) = true;  scheme{k} = 'legacy_rs3';
    end
end

report = table(name(:), present(:), scheme(:), fpath_current(:), fpath_legacy(:), ...
    'VariableNames', {'family', 'present', 'scheme', 'fpath_current', 'fpath_legacy'});

all_present = all(present);

if ~all_present
    missing = report.family(~report.present);
    fprintf('[atlas_check_cache_exists] Missing cached atlas for %d/%d families ' ...
        '(checked BOTH current and legacy_rs3 fingerprint schemes):\n', ...
        numel(missing), K);
    for k = 1:numel(missing)
        fprintf('    - %s\n', missing{k});
    end
    fprintf('  (cfg.cache.enable=%d, cfg.cache.dir=%s)\n', ...
        cfg.cache.enable, cfg.cache.dir);
else
    nLegacy = sum(strcmp(report.scheme, 'legacy_rs3'));
    fprintf('[atlas_check_cache_exists] All %d family atlases found in cache:\n  %s\n', ...
        K, cfg.cache.dir);
    if nLegacy > 0
        fprintf('  (%d/%d hit via the legacy_rs3 fingerprint scheme)\n', nLegacy, K);
    end
end

end
