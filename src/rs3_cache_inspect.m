function entries = rs3_cache_inspect(cacheDir, verbose)
%RS3_CACHE_INSPECT  Scan a cache directory and return parsed metadata for every entry.
%
% Usage:
%   entries = rs3_cache_inspect()                    % uses default cache dir
%   entries = rs3_cache_inspect(cacheDir)
%   entries = rs3_cache_inspect(cacheDir, false)     % suppress table printout
%
% Returns struct array with one entry per .mat file.  Key fields:
%   .idx            row index (use for SOURCE_IDX in cache manager)
%   .family         family name
%   .config_key     8-char hash of grid+fan+propag+log config (same = same build params)
%   .Tmax / .Tmax_str
%   .DV_cap_nd, .dtheta_fan_deg
%   .dx, .dy, .dtheta_deg, .ds_seed
%   .absTol, .relTol, .step_len_factor, .maxstep_factor
%   .rows_FRS, .rows_BRS
%   .MB             file size (MB)
%   .created        creation timestamp
%   .parsed         full parsed fingerprint struct
%   .valid / .error

if nargin < 1 || isempty(cacheDir)
    cacheDir = fullfile(rs3_repo_root(), 'rs3_cache');
end
if nargin < 2
    verbose = true;
end

mats = dir(fullfile(cacheDir, '*.mat'));
if isempty(mats)
    entries = [];
    if verbose
        fprintf('[rs3_cache_inspect] No .mat files in: %s\n', cacheDir);
    end
    return;
end

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
            e.error = 'no cache_meta'; entries(k) = e; continue;
        end
        meta = tmp.cache_meta;
        e.created     = local_getf(meta, 'created', '');
        e.version_tag = local_getf(meta, 'version_tag', '');

        fp = local_getf(meta, 'fingerprint', '');
        if isempty(fp)
            e.error = 'empty fingerprint'; entries(k) = e; continue;
        end

        p = rs3_cache_fingerprint_parse(fp);
        e.parsed = p;

        % --- core params ---
        e.family           = p.family;
        e.Tmax             = p.propag.Tmax;
        e.Tmax_str         = local_tmax_str(p.propag.Tmax);
        e.DV_cap_nd        = p.fan.DV_cap_nd;
        e.dtheta_fan_deg   = rad2deg(p.fan.dtheta_fan);
        e.dx               = p.grid.dx;
        e.dy               = p.grid.dy;
        e.dtheta_deg       = rad2deg(p.grid.dtheta);
        e.ds_seed          = p.seed.ds_seed;
        e.absTol           = p.propag.absTol;
        e.relTol           = p.propag.relTol;
        e.step_len_factor  = p.log.step_len_factor;
        e.maxstep_factor   = p.log.maxstep_factor;

        % --- config_key: hash of everything EXCEPT family/mu/CJ ---
        % Entries with the same config_key used the same build settings.
        parts = strsplit(fp, '|');
        if numel(parts) >= 5
            config_str  = strjoin(parts(5:end), '|');
            full_hash   = rs3_md5(config_str);
            e.config_key = full_hash(1:8);
        end

        % --- row counts ---
        if isfield(meta, 'stats')
            st = meta.stats;
            e.rows_FRS = local_getf(st, 'rows_FRS_upper', NaN);
            e.rows_BRS = local_getf(st, 'rows_BRS_upper', NaN);
        end
        e.valid = true;

    catch ME
        e.valid = false;
        e.error = ME.message;
    end
    entries(k) = e;
end

if ~verbose, return; end

% ---- print table ----
fprintf('\n%s\n', repmat('=', 1, 118));
fprintf('  Cache directory: %s\n', cacheDir);
fprintf('  %d entries found\n', N);
fprintf('%s\n', repmat('=', 1, 118));

hdr = sprintf('  %-3s  %-28s  %-5s  %-5s  %-6s  %-6s  %-7s  %-6s  %-7s  %-7s  %-8s  %-8s  %s', ...
    '#', 'Family', 'Tmax', 'DV', 'dfan°', 'dth°', 'dx', 'ds', 'absTol', 'slf/msf', 'FRS', 'BRS', 'Created (config_key)');
fprintf('%s\n', hdr);
fprintf('  %s\n', repmat('-', 1, 114));

for k = 1:N
    e = entries(k);
    if ~e.valid
        fprintf('  %-3d  [INVALID: %s]  %s\n', k, e.error, e.fname);
        continue;
    end
    fprintf('  %-3d  %-28s  %-5s  %-5.2f  %-6.2f  %-6.2f  %-7.4f  %-6.4f  %-7s  %-7s  %-8s  %-8s  %s (%s)\n', ...
        k, ...
        e.family, ...
        e.Tmax_str, ...
        e.DV_cap_nd, ...
        e.dtheta_fan_deg, ...
        e.dtheta_deg, ...
        e.dx, ...
        e.ds_seed, ...
        local_fmt_tol(e.absTol), ...
        sprintf('%.2f/%.2f', e.step_len_factor, e.maxstep_factor), ...
        local_fmt_rows(e.rows_FRS), ...
        local_fmt_rows(e.rows_BRS), ...
        e.created, ...
        e.config_key);
end
fprintf('%s\n\n', repmat('=', 1, 118));

% Print config-key grouping summary
keys = {entries(logical([entries.valid])).config_key};
[ukeys, ~, kidx] = unique(keys);
if numel(ukeys) > 1
    fprintf('  Config-key groups (%d distinct configurations):\n', numel(ukeys));
    valid_entries = entries(logical([entries.valid]));
    for g = 1:numel(ukeys)
        gmask = kidx == g;
        gfams = {valid_entries(gmask).family};
        gentry = valid_entries(find(gmask,1));
        fprintf('    [%s] dx=%.4f dtheta=%.2f° Tmax=%s absTol=%s slf=%.2f | %d entries: %s\n', ...
            ukeys{g}, gentry.dx, gentry.dtheta_deg, gentry.Tmax_str, ...
            local_fmt_tol(gentry.absTol), gentry.step_len_factor, ...
            sum(gmask), strjoin(gfams, ', '));
    end
    fprintf('\n');
end
end

% =========================================================================
function e = local_empty_entry()
e = struct('idx',NaN,'fname','','fpath','','bytes',NaN,'MB',NaN, ...
    'valid',false,'error','', ...
    'family','','created','','version_tag','','config_key','', ...
    'Tmax',NaN,'Tmax_str','','DV_cap_nd',NaN,'dtheta_fan_deg',NaN, ...
    'dx',NaN,'dy',NaN,'dtheta_deg',NaN,'ds_seed',NaN, ...
    'absTol',NaN,'relTol',NaN,'step_len_factor',NaN,'maxstep_factor',NaN, ...
    'rows_FRS',NaN,'rows_BRS',NaN,'parsed',struct());
end

function v = local_getf(s, f, def)
if isfield(s,f), v = s.(f); else, v = def; end
end

function s = local_tmax_str(T)
if isnan(T), s = '?'; return; end
frac = T/pi;
if     abs(frac-1)   < 1e-6, s = 'π';
elseif abs(frac-0.5) < 1e-6, s = 'π/2';
elseif abs(frac-2)   < 1e-6, s = '2π';
elseif abs(frac-1.5) < 1e-6, s = '3π/2';
elseif abs(frac-0.25)< 1e-6, s = 'π/4';
else,                          s = sprintf('%.2fπ', frac);
end
end

function s = local_fmt_rows(n)
if isnan(n),    s = '?';
elseif n>=1e6,  s = sprintf('%.1fM', n/1e6);
elseif n>=1e3,  s = sprintf('%.1fk', n/1e3);
else,           s = sprintf('%d', round(n));
end
end

function s = local_fmt_tol(t)
if isnan(t), s = '?'; return; end
e = floor(log10(t));
s = sprintf('1e%d', e);
end
