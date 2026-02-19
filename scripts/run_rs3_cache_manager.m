%% RUN_RS3_CACHE_MANAGER
% Interactive cache inspector + subset deriver.
%
% WORKFLOW:
%   Step 1 — Run the script once to print the cache table (ACTION = 'inspect').
%   Step 2 — Edit the USER KNOBS: set ACTION='derive', choose SOURCE_FAMILIES
%            and set subset parameters.
%   Step 3 — Run again. Derived caches are saved and immediately usable by
%            the batch runner when you set cfg to match SUBSET_CFG.
%
% The derived caches are real cache files — run_rs4_all_pairs_summary with
% cfg matching SUBSET_CFG will find them as cache hits (no rebuild needed).

clear; clc;

% --- hard-pin repo paths ---
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
restoredefaultpath; rehash toolboxcache;
addpath(repoRoot);
addpath(fullfile(repoRoot,'src'));
addpath(fullfile(repoRoot,'scripts'));
rehash;

% =====================================================================
% ========================= USER KNOBS ================================
% =====================================================================

% What to do:
%   'inspect'  — just print the cache table (read-only, safe to run anytime)
%   'derive'   — derive subset caches from selected source entries
ACTION = 'inspect';

% Cache directory to scan (leave empty for default)
CACHE_DIR = '';   % '' = repo/rs3_cache

% -- DERIVE options (only used when ACTION = 'derive') ----------------

% Which families to process. Use either:
%   [] = all families found in cache
%   {'Lyapunov L1', 'Lyapunov L2', ...} = specific subset
SOURCE_FAMILIES = {};   % {} = all

% Source selection strategy (when multiple caches exist for the same family):
%   'largest_Tmax'  — use the entry with the biggest Tmax as the source
%   'latest'        — use the most recently created entry
%   'all'           — derive from every matching entry (may produce duplicates)
SOURCE_STRATEGY = 'largest_Tmax';

% Desired subset configuration.
% Only set the fields you want to RESTRICT. Any field stricter than the
% source will be applied; if a field equals the source it's a no-op.
% Fields not listed here inherit from the source.
SUBSET_CFG = struct();
SUBSET_CFG.propag.Tmax      = pi/2;          % new max integration time (nd)
SUBSET_CFG.fan.DV_cap_nd    = 0.10;          % new DV cap (nd)
% SUBSET_CFG.grid.dtheta    = deg2rad(4);    % (optional) coarser theta grid
% SUBSET_CFG.fan.dtheta_fan = deg2rad(2);    % (optional) coarser fan — NOT applied by
%                                            %   subset (dtheta_fan just rescales DV_cap)

% =====================================================================
% =====================================================================

% --- resolve cache dir ---
if isempty(CACHE_DIR)
    CACHE_DIR = fullfile(rs3_repo_root(), 'rs3_cache');
end
if ~exist(CACHE_DIR, 'dir')
    error('Cache directory does not exist: %s', CACHE_DIR);
end

% =====================================================================
%  INSPECT — scan and print
% =====================================================================
fprintf('[cache_mgr] Scanning: %s\n', CACHE_DIR);
entries = rs3_cache_inspect(CACHE_DIR, true);

if isempty(entries)
    fprintf('[cache_mgr] Nothing to do — cache is empty.\n');
    return;
end

if strcmpi(ACTION, 'inspect')
    fprintf('[cache_mgr] ACTION=inspect done. Edit ACTION=''derive'' to create subsets.\n');
    return;
end

if ~strcmpi(ACTION, 'derive')
    error('[cache_mgr] Unknown ACTION: ''%s''. Use ''inspect'' or ''derive''.', ACTION);
end

% =====================================================================
%  DERIVE — build subset caches
% =====================================================================

% Filter to valid entries only
valid_idx = find([entries.valid]);
if isempty(valid_idx)
    error('[cache_mgr] No valid cache entries found.');
end

% Filter by family if requested
if ~isempty(SOURCE_FAMILIES)
    keep = false(size(valid_idx));
    for k = 1:numel(valid_idx)
        keep(k) = any(strcmpi(entries(valid_idx(k)).family, SOURCE_FAMILIES));
    end
    valid_idx = valid_idx(keep);
    if isempty(valid_idx)
        error('[cache_mgr] No valid entries found for the requested families.');
    end
end

% Group by family, select source per strategy
families_present = unique({entries(valid_idx).family});
fprintf('\n[cache_mgr] Families to process (%d): %s\n', numel(families_present), ...
    strjoin(families_present, ', '));

% Build full cfg from the source atlas's parsed fingerprint + SUBSET_CFG overrides
cfg_base = rs3_cfg_defaults();

nDerived = 0;
nSkipped = 0;

for fi = 1:numel(families_present)
    fam = families_present{fi};

    % Find all entries for this family
    fam_mask = strcmpi({entries(valid_idx).family}, fam);
    fam_idx  = valid_idx(fam_mask);

    % Select source
    src_entry = local_select_source(entries(fam_idx), SOURCE_STRATEGY);
    if isempty(src_entry)
        fprintf('[cache_mgr] SKIP %s — could not select source entry.\n', fam);
        nSkipped = nSkipped + 1;
        continue;
    end

    fprintf('\n[cache_mgr] === %s ===\n', fam);
    fprintf('[cache_mgr]   source: %s (Tmax=%s, DV_cap=%.3f)\n', ...
        src_entry.fname, src_entry.Tmax_str, src_entry.DV_cap_nd);

    % Build sub-config: start from base, apply parsed source params, then SUBSET_CFG
    cfg_src = local_apply_parsed(cfg_base, src_entry.parsed);
    cfg_sub = local_apply_subset_cfg(cfg_src, SUBSET_CFG);

    % Check if sub-config is actually different from source
    sub_fp = rs3_cache_fingerprint_family(fam, src_entry.parsed.mu, src_entry.parsed.CJ, cfg_sub);
    src_fp = src_entry.parsed.raw;
    if strcmp(sub_fp, src_fp)
        fprintf('[cache_mgr]   SKIP (sub-config identical to source — nothing to derive)\n');
        nSkipped = nSkipped + 1;
        continue;
    end

    % Check if derived cache already exists
    info_sub = rs3_cache_get_path(fam, src_entry.parsed.mu, src_entry.parsed.CJ, cfg_sub);
    if exist(info_sub.fpath, 'file') == 2 && ~local_force_rebuild()
        fprintf('[cache_mgr]   SKIP (derived cache already exists: %s)\n', info_sub.fname);
        nSkipped = nSkipped + 1;
        continue;
    end

    % Load full atlas from source
    fprintf('[cache_mgr]   Loading source atlas (%.1f MB) ...\n', src_entry.MB);
    tmp = load(src_entry.fpath, 'S');
    if ~isfield(tmp, 'S')
        fprintf('[cache_mgr]   ERROR: no S field in %s\n', src_entry.fname);
        nSkipped = nSkipped + 1;
        continue;
    end
    S_src = tmp.S;
    % Ensure grid3 helpers are rebuilt (function handles not saved)
    S_src.grid3 = rs3_grid_make(cfg_src);

    % Derive subset
    try
        S_sub = rs3_atlas_derive_subset(S_src, cfg_sub);
    catch ME
        fprintf('[cache_mgr]   ERROR during derive: %s\n', ME.message);
        nSkipped = nSkipped + 1;
        continue;
    end

    % Save as new cache
    info = rs3_cache_save_family(S_sub, cfg_sub);
    fprintf('[cache_mgr]   Saved: %s\n', info.fpath);
    nDerived = nDerived + 1;
end

fprintf('\n[cache_mgr] Done. %d derived, %d skipped.\n', nDerived, nSkipped);

if nDerived > 0
    fprintf('\n[cache_mgr] To use the derived caches, run the batch runner with:\n');
    fprintf('   cfg.propag.Tmax   = %.6g;\n', SUBSET_CFG.propag.Tmax);
    if isfield(SUBSET_CFG, 'fan') && isfield(SUBSET_CFG.fan, 'DV_cap_nd')
        fprintf('   cfg.fan.DV_cap_nd = %.6g;\n', SUBSET_CFG.fan.DV_cap_nd);
    end
    if isfield(SUBSET_CFG, 'grid') && isfield(SUBSET_CFG.grid, 'dtheta')
        fprintf('   cfg.grid.dtheta   = deg2rad(%.4g);\n', rad2deg(SUBSET_CFG.grid.dtheta));
    end
    fprintf('   (other params inherited from defaults / batch runner knobs)\n\n');
end

% =====================================================================
%  Local helpers
% =====================================================================

function src = local_select_source(fam_entries, strategy)
src = [];
if isempty(fam_entries), return; end
valid = fam_entries([fam_entries.valid]);
if isempty(valid), return; end

switch lower(strategy)
    case 'largest_tmax'
        tmaxVals = [valid.Tmax];
        [~, idx] = max(tmaxVals);
        src = valid(idx);
    case 'latest'
        % Sort by created string (ISO format sorts lexicographically)
        dates = {valid.created};
        [~, idx] = sort(dates);
        src = valid(idx(end));
    case 'all'
        src = valid;   % caller handles array
    otherwise
        src = valid(1);
end
end

function cfg_out = local_apply_parsed(cfg_in, p)
% Apply parsed fingerprint parameters to cfg.
cfg_out = cfg_in;
if ~isstruct(p), return; end

if isfield(p,'grid')
    if isfield(p.grid,'R'),      cfg_out.grid.Rdom    = p.grid.R;       end
    if isfield(p.grid,'dx'),     cfg_out.grid.dx      = p.grid.dx;      end
    if isfield(p.grid,'dy'),     cfg_out.grid.dy      = p.grid.dy;      end
    if isfield(p.grid,'dtheta'), cfg_out.grid.dtheta  = p.grid.dtheta;  end
end
if isfield(p,'seed')
    if isfield(p.seed,'Tf_scale'),  cfg_out.seed.Tf_scale  = p.seed.Tf_scale;  end
    if isfield(p.seed,'N_dense'),   cfg_out.seed.N_dense   = round(p.seed.N_dense); end
    if isfield(p.seed,'ds_seed'),   cfg_out.seed.ds_seed   = p.seed.ds_seed;   end
    if isfield(p.seed,'y_eps'),     cfg_out.seed.y_eps     = p.seed.y_eps;     end
    if isfield(p.seed,'minSegPts'), cfg_out.seed.minSegPts = round(p.seed.minSegPts); end
end
if isfield(p,'fan')
    if isfield(p.fan,'dtheta_fan'), cfg_out.fan.dtheta_fan = p.fan.dtheta_fan; end
    if isfield(p.fan,'DV_cap_nd'),  cfg_out.fan.DV_cap_nd  = p.fan.DV_cap_nd;  end
end
if isfield(p,'propag')
    if isfield(p.propag,'Tmax'),   cfg_out.propag.Tmax   = p.propag.Tmax;   end
    if isfield(p.propag,'relTol'), cfg_out.propag.relTol = p.propag.relTol; end
    if isfield(p.propag,'absTol'), cfg_out.propag.absTol = p.propag.absTol; end
    if isfield(p.propag,'v2tol'),  cfg_out.propag.v2tol  = p.propag.v2tol;  end
end
if isfield(p,'log')
    if isfield(p.log,'segwalk_enable'),  cfg_out.log.segwalk.enable      = logical(p.log.segwalk_enable); end
    if isfield(p.log,'segwalk_frac'),    cfg_out.log.segwalk.frac        = p.log.segwalk_frac;   end
    if isfield(p.log,'step_len_factor'), cfg_out.log.step_len_factor     = p.log.step_len_factor; end
    if isfield(p.log,'maxstep_factor'),  cfg_out.log.maxstep_factor      = p.log.maxstep_factor;  end
end
if isfield(p,'version_tag') && ~isempty(p.version_tag)
    cfg_out.cache.version_tag = p.version_tag;
end
end

function cfg_out = local_apply_subset_cfg(cfg_in, sub)
% Override cfg fields from the SUBSET_CFG struct (deep merge).
cfg_out = cfg_in;
if ~isstruct(sub), return; end
fields = fieldnames(sub);
for i = 1:numel(fields)
    f = fields{i};
    if isstruct(sub.(f)) && isfield(cfg_out, f) && isstruct(cfg_out.(f))
        cfg_out.(f) = local_apply_subset_cfg(cfg_out.(f), sub.(f));
    else
        cfg_out.(f) = sub.(f);
    end
end
end

function tf = local_force_rebuild()
% Set to true to always overwrite existing derived caches.
tf = false;
end
