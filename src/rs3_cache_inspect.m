function entries = rs3_cache_inspect(cacheDir, verbose)
%RS3_CACHE_INSPECT  Scan a cache directory and return parsed metadata for every entry.
%
% Usage:
%   entries = rs3_cache_inspect()                    % uses default cache dir
%   entries = rs3_cache_inspect(cacheDir)
%   entries = rs3_cache_inspect(cacheDir, false)     % suppress table printout
%
% Returns:
%   entries  — struct array, one entry per .mat file with fields:
%       .idx           row index (for selection)
%       .fname         filename
%       .fpath         full path
%       .bytes         file size in bytes
%       .MB            file size in MB
%       .family        family name
%       .created       creation timestamp string
%       .version_tag   cache version tag
%       .Tmax          propag.Tmax (nd)
%       .Tmax_str      e.g. 'π' or 'π/2'
%       .DV_cap_nd     fan.DV_cap_nd
%       .dtheta_fan_deg fan.dtheta_fan in degrees
%       .dx            grid.dx
%       .dy            grid.dy
%       .dtheta_deg    grid.dtheta in degrees
%       .ds_seed       seed.ds_seed
%       .absTol        propag.absTol
%       .rows_FRS      rows in FRS upper (from cache_meta.stats)
%       .rows_BRS      rows in BRS upper (from cache_meta.stats)
%       .parsed        full parsed-fingerprint struct (from rs3_cache_fingerprint_parse)
%       .valid         true if cache_meta loaded and fingerprint parsed OK
%       .error         error message if valid=false

if nargin < 1 || isempty(cacheDir)
    cacheDir = fullfile(rs3_repo_root(), 'rs3_cache');
end
if nargin < 2
    verbose = true;
end

% --- scan ---
mats = dir(fullfile(cacheDir, '*.mat'));
if isempty(mats)
    entries = [];
    if verbose
        fprintf('[rs3_cache_inspect] No .mat files found in: %s\n', cacheDir);
    end
    return;
end

% Pre-allocate
N = numel(mats);
entries = repmat(local_empty_entry(), N, 1);

for k = 1:N
    e = local_empty_entry();
    e.idx   = k;
    e.fname = mats(k).name;
    e.fpath = fullfile(cacheDir, mats(k).name);
    e.bytes = mats(k).bytes;
    e.MB    = mats(k).bytes / 1e6;

    try
        tmp = load(e.fpath, 'cache_meta');
        if ~isfield(tmp, 'cache_meta')
            e.valid = false;
            e.error = 'no cache_meta field';
            entries(k) = e;
            continue;
        end
        meta = tmp.cache_meta;
        e.created     = local_getfield(meta, 'created', '');
        e.version_tag = local_getfield(meta, 'version_tag', '');

        % parse fingerprint
        fp = local_getfield(meta, 'fingerprint', '');
        if isempty(fp)
            e.valid = false;
            e.error = 'empty fingerprint';
            entries(k) = e;
            continue;
        end

        p = rs3_cache_fingerprint_parse(fp);
        e.parsed        = p;
        e.family        = p.family;
        e.Tmax          = p.propag.Tmax;
        e.Tmax_str      = local_tmax_str(p.propag.Tmax);
        e.DV_cap_nd     = p.fan.DV_cap_nd;
        e.dtheta_fan_deg = rad2deg(p.fan.dtheta_fan);
        e.dx            = p.grid.dx;
        e.dy            = p.grid.dy;
        e.dtheta_deg    = rad2deg(p.grid.dtheta);
        e.ds_seed       = p.seed.ds_seed;
        e.absTol        = p.propag.absTol;

        % row counts from stats
        if isfield(meta, 'stats')
            st = meta.stats;
            e.rows_FRS = local_getfield(st, 'rows_FRS_upper', NaN);
            e.rows_BRS = local_getfield(st, 'rows_BRS_upper', NaN);
        end
        e.valid = true;

    catch ME
        e.valid = false;
        e.error = ME.message;
    end

    entries(k) = e;
end

% --- print table ---
if verbose
    fprintf('\n%s\n', repmat('=', 1, 110));
    fprintf('  Cache directory: %s\n', cacheDir);
    fprintf('  %d entries found\n', N);
    fprintf('%s\n', repmat('=', 1, 110));
    fprintf('  %-3s  %-28s  %-8s  %-6s  %-8s  %-8s  %-7s  %-7s  %-7s  %-8s  %-8s  %s\n', ...
        '#', 'Family', 'Tmax', 'DV_cap', 'dfan(°)', 'dth(°)', 'dx', 'dy', 'ds', 'FRS rows', 'BRS rows', 'Created');
    fprintf('  %s\n', repmat('-', 1, 106));
    for k = 1:N
        e = entries(k);
        if ~e.valid
            fprintf('  %-3d  [INVALID: %s]  %s\n', k, e.error, e.fname);
            continue;
        end
        fprintf('  %-3d  %-28s  %-8s  %-6.3f  %-8.2f  %-8.2f  %-7.4f  %-7.4f  %-7.4f  %-8s  %-8s  %s\n', ...
            k, ...
            e.family, ...
            e.Tmax_str, ...
            e.DV_cap_nd, ...
            e.dtheta_fan_deg, ...
            e.dtheta_deg, ...
            e.dx, ...
            e.dy, ...
            e.ds_seed, ...
            local_fmt_rows(e.rows_FRS), ...
            local_fmt_rows(e.rows_BRS), ...
            e.created);
    end
    fprintf('%s\n\n', repmat('=', 1, 110));
end

end

% =========================================================================
function e = local_empty_entry()
e = struct( ...
    'idx', NaN, 'fname', '', 'fpath', '', 'bytes', NaN, 'MB', NaN, ...
    'valid', false, 'error', '', ...
    'family', '', 'created', '', 'version_tag', '', ...
    'Tmax', NaN, 'Tmax_str', '', ...
    'DV_cap_nd', NaN, 'dtheta_fan_deg', NaN, ...
    'dx', NaN, 'dy', NaN, 'dtheta_deg', NaN, 'ds_seed', NaN, 'absTol', NaN, ...
    'rows_FRS', NaN, 'rows_BRS', NaN, ...
    'parsed', struct());
end

function v = local_getfield(s, f, def)
if isfield(s, f), v = s.(f); else, v = def; end
end

function s = local_tmax_str(T)
% Pretty-print Tmax as multiples of pi.
if isnan(T)
    s = 'NaN';
    return;
end
frac = T / pi;
if     abs(frac - 1)   < 1e-6, s = 'π';
elseif abs(frac - 0.5) < 1e-6, s = 'π/2';
elseif abs(frac - 2)   < 1e-6, s = '2π';
elseif abs(frac - 1.5) < 1e-6, s = '3π/2';
elseif abs(frac - 0.25)< 1e-6, s = 'π/4';
else,                           s = sprintf('%.3fπ', frac);
end
end

function s = local_fmt_rows(n)
if isnan(n)
    s = '?';
elseif n >= 1e6
    s = sprintf('%.1fM', n/1e6);
elseif n >= 1e3
    s = sprintf('%.1fk', n/1e3);
else
    s = sprintf('%d', round(n));
end
end
