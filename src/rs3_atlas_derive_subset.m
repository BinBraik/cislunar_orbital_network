function S_sub = rs3_atlas_derive_subset(S, cfg_sub)
%RS3_ATLAS_DERIVE_SUBSET  Derive a leaner atlas from a cached "fat" atlas without re-running ODE.
%
% Given a cached atlas S built with a "fat" config (large Tmax, large DV_cap),
% produces a new atlas S_sub that is valid for a STRICTER config cfg_sub:
%
%   Tmax_sub   <= Tmax_super   : drops rows where abs(t) > Tmax_sub
%   DV_cap_sub <= DV_cap_super : drops headings (per seed) beyond new DV budget
%   Grid coarsening (optional) : re-bins rows into coarser grid when
%                                cfg_sub.grid.{dx,dy,dtheta} > original
%
% S_sub has its own fingerprint matching cfg_sub and can be saved as a new
% cache entry via rs3_cache_save_family(S_sub, cfg_sub).
%
% Usage:
%   S_sub = rs3_atlas_derive_subset(S, cfg_sub)
%   info  = rs3_cache_save_family(S_sub, cfg_sub)   % cache it
%
% Limitations:
%   - Cannot expand: Tmax_sub > Tmax_super or DV_cap_sub > DV_cap_super → error
%   - Cannot refine grid: cfg_sub.grid.dx < original dx → error
%   - Tolerances (absTol, relTol) are inherited from the parent; the subset
%     is physically valid but the fingerprint will carry cfg_sub's tolerance values.

assert(isstruct(S) && isfield(S,'Step4'), 'S must be a valid atlas struct with Step4.');
assert(isstruct(cfg_sub), 'cfg_sub must be a config struct.');

grid3_orig = S.grid3;
Tmax_orig  = S.Step4.Tmax;
VU_mps     = local_cfg_get(cfg_sub, 'units.VU_mps',  1.0);

Tmax_sub   = cfg_sub.propag.Tmax;
DV_sub     = cfg_sub.fan.DV_cap_nd;
DV_orig    = S.Step4.fanStats.DV_cap_nd;   % may not exist; fall back gracefully

% --- check feasibility ---
if Tmax_sub > Tmax_orig * (1 + 1e-9)
    error('rs3:subset:expansion', ...
        'Cannot subset: Tmax_sub (%.4g) > Tmax_orig (%.4g). Use expansion instead.', ...
        Tmax_sub, Tmax_orig);
end

if isfield(S.Step4, 'fanStats') && isfield(S.Step4.fanStats, 'DV_cap_nd')
    DV_orig = S.Step4.fanStats.DV_cap_nd;
else
    DV_orig = Inf;   % unknown — allow
end

if DV_sub > DV_orig * (1 + 1e-9)
    error('rs3:subset:expansion', ...
        'Cannot subset: DV_cap_sub (%.4g) > DV_cap_orig (%.4g). Use expansion instead.', ...
        DV_sub, DV_orig);
end

% --- grid change check ---
doRegrid = false;
dx_sub  = cfg_sub.grid.dx;
dy_sub  = cfg_sub.grid.dy;
dth_sub = cfg_sub.grid.dtheta;
dx_orig  = grid3_orig.dx;
dy_orig  = grid3_orig.dy;
dth_orig = grid3_orig.dtheta;

tol = 1e-10;
if dx_sub  < dx_orig  - tol || dy_sub  < dy_orig  - tol || dth_sub < dth_orig - tol
    error('rs3:subset:refine', ...
        ['Cannot refine grid (sub must be coarser than orig).\n' ...
         '  orig: dx=%.4g dy=%.4g dtheta=%.4g°\n' ...
         '  sub:  dx=%.4g dy=%.4g dtheta=%.4g°'], ...
        dx_orig, dy_orig, rad2deg(dth_orig), dx_sub, dy_sub, rad2deg(dth_sub));
end
gridChanged = abs(dx_sub - dx_orig) > tol || abs(dy_sub - dy_orig) > tol || abs(dth_sub - dth_orig) > tol;
if gridChanged
    doRegrid = true;
    fprintf('[rs3_subset] Grid coarsening requested: dx %.4g→%.4g, dy %.4g→%.4g, dtheta %.2f°→%.2f°\n', ...
        dx_orig, dx_sub, dy_orig, dy_sub, rad2deg(dth_orig), rad2deg(dth_sub));
end

fprintf('[rs3_subset] Source: Tmax=%.4g, DV_cap=%.4g, dtheta=%.2f°  →  Sub: Tmax=%.4g, DV_cap=%.4g, dtheta=%.2f°\n', ...
    Tmax_orig, DV_orig, rad2deg(dth_orig), Tmax_sub, DV_sub, rad2deg(dth_sub));

% =====================================================================
%  1. DV_cap filter on delta_lists + build a per-seed DV mask
% =====================================================================
delta_lists_orig = S.Step4.delta_lists;
Nseeds = numel(delta_lists_orig);

% Compute v0 per seed (vectorized over SeedsUpper)
pot_U = rs3_core_cr3bp_U_and_derivs(S.SeedsUpper(:,1), S.SeedsUpper(:,2), S.mu);
v0_upper = sqrt(max(2*pot_U.U - S.CJ, 0));

% For each seed: valid heading mask
% valid_head{s}(h) = true if heading h of seed s is within DV_sub
valid_head = cell(Nseeds, 1);
delta_lists_sub = cell(Nseeds, 1);
for s = 1:Nseeds
    v0s = v0_upper(s);
    if v0s < 1e-12
        delta_max = 0;
    else
        delta_max = 2 * asin(min(1, DV_sub / (2 * v0s)));
    end
    d = double(delta_lists_orig{s});
    mask = abs(d) <= delta_max + 1e-12;
    valid_head{s} = mask;
    delta_lists_sub{s} = d(mask);   % trimmed list (re-indexed!)
end

% Build iHead remapping: old iHead → new iHead (0 = dropped)
% new_ihead(s, h_old) = new index in delta_lists_sub{s}, or 0 if dropped
ihead_remap = cell(Nseeds, 1);
for s = 1:Nseeds
    mask = valid_head{s};
    nH = numel(mask);
    remap = zeros(nH, 1, 'uint16');
    newIdx = uint16(0);
    for h = 1:nH
        if mask(h)
            newIdx = newIdx + uint16(1);
            remap(h) = newIdx;
        end
    end
    ihead_remap{s} = remap;
end

% =====================================================================
%  2. Filter rows: Tmax + DV_cap
% =====================================================================
rows_FRS_sub = local_filter_rows(S.Step4.rows_FRS_upper, Tmax_sub, valid_head, ihead_remap);
rows_BRS_sub = local_filter_rows(S.Step4.rows_BRS_upper, Tmax_sub, valid_head, ihead_remap);

fprintf('[rs3_subset] FRS: %d → %d rows (%.1f%%)\n', ...
    double(S.Step4.rows_FRS_upper.n), double(rows_FRS_sub.n), ...
    100*double(rows_FRS_sub.n)/max(1,double(S.Step4.rows_FRS_upper.n)));
fprintf('[rs3_subset] BRS: %d → %d rows (%.1f%%)\n', ...
    double(S.Step4.rows_BRS_upper.n), double(rows_BRS_sub.n), ...
    100*double(rows_BRS_sub.n)/max(1,double(S.Step4.rows_BRS_upper.n)));

% =====================================================================
%  3. (Optional) Re-grid to coarser grid
% =====================================================================
if doRegrid
    grid3_sub = rs3_grid_make(cfg_sub);
    rows_FRS_sub = local_regrid_rows(rows_FRS_sub, grid3_orig, grid3_sub);
    rows_BRS_sub = local_regrid_rows(rows_BRS_sub, grid3_orig, grid3_sub);
    fprintf('[rs3_subset] After re-grid: FRS %d rows, BRS %d rows\n', ...
        double(rows_FRS_sub.n), double(rows_BRS_sub.n));
else
    grid3_sub = grid3_orig;
end

% =====================================================================
%  4. Assemble S_sub
% =====================================================================
S_sub = S;   % inherit identity fields (name, mu, CJ, SeedsUpper, SeedsLower, PO_xy, ...)
S_sub.grid3 = grid3_sub;

step4_sub = S.Step4;
step4_sub.Tmax              = Tmax_sub;
step4_sub.rows_FRS_upper    = rows_FRS_sub;
step4_sub.rows_BRS_upper    = rows_BRS_sub;
step4_sub.delta_lists       = delta_lists_sub;
step4_sub.nJobs             = uint32(Nseeds * max(cellfun(@numel, delta_lists_sub)));

% Update fanStats if present
if isfield(step4_sub, 'fanStats')
    step4_sub.fanStats.DV_cap_nd = DV_sub;
end

S_sub.Step4 = step4_sub;

fprintf('[rs3_subset] Done. Derived atlas ready (family: %s).\n', S_sub.name);
end

% =========================================================================
function rows_out = local_filter_rows(rows_in, Tmax_sub, valid_head, ihead_remap)
% Apply Tmax and DV_cap masks to a packed rows struct.
n = double(rows_in.n);
if n == 0
    rows_out = rows_in;
    return;
end

t    = double(rows_in.t(1:n));
iS   = double(rows_in.iSeed(1:n));
iH   = double(rows_in.iHead(1:n));

% Tmax mask
mask_t = abs(t) <= Tmax_sub + 1e-12;

% DV_cap mask: per row, check if its heading is still valid
mask_dv = false(n, 1);
for r = 1:n
    s = iS(r);  h = iH(r);
    if s >= 1 && s <= numel(valid_head)
        vh = valid_head{s};
        if h >= 1 && h <= numel(vh)
            mask_dv(r) = vh(h);
        end
    end
end

keep = mask_t & mask_dv;
nKeep = sum(keep);

rows_out = rs3_rows_empty(nKeep);
rows_out.n = uint32(nKeep);
if nKeep == 0, return; end

idx = find(keep);
rows_out.iSeed(1:nKeep)    = rows_in.iSeed(idx);
rows_out.iHead(1:nKeep)    = local_remap_ihead(rows_in.iSeed(idx), rows_in.iHead(idx), ihead_remap);
rows_out.leg(1:nKeep)      = rows_in.leg(idx);
rows_out.halfFlag(1:nKeep) = rows_in.halfFlag(idx);
rows_out.t(1:nKeep)        = rows_in.t(idx);
rows_out.ix(1:nKeep)       = rows_in.ix(idx);
rows_out.iy(1:nKeep)       = rows_in.iy(idx);
rows_out.it(1:nKeep)       = rows_in.it(idx);
end

function new_ihead = local_remap_ihead(iSeed_kept, iHead_kept, ihead_remap)
% Remap iHead values after DV_cap filtering (delta_lists are now re-indexed).
n = numel(iSeed_kept);
new_ihead = zeros(n, 1, 'uint16');
for r = 1:n
    s = double(iSeed_kept(r));
    h = double(iHead_kept(r));
    if s >= 1 && s <= numel(ihead_remap)
        rm = ihead_remap{s};
        if h >= 1 && h <= numel(rm)
            new_ihead(r) = rm(h);
        end
    end
end
end

% =========================================================================
function rows_out = local_regrid_rows(rows_in, grid3_orig, grid3_new)
% Re-bin rows from original grid to a new (coarser) grid.
% Uses voxel-center approximation; keeps only first hit per (seed,head,new_voxel).
n = double(rows_in.n);
if n == 0
    rows_out = rows_in;
    return;
end

ix = double(rows_in.ix(1:n));
iy = double(rows_in.iy(1:n));
it = double(rows_in.it(1:n));

% look up voxel centers in original grid (clamp to valid range)
ix = max(1, min(numel(grid3_orig.x_centers), ix));
iy = max(1, min(numel(grid3_orig.y_centers), iy));
it = max(1, min(numel(grid3_orig.th_centers), it));

x_c  = grid3_orig.x_centers(ix);
y_c  = grid3_orig.y_centers(iy);
th_c = grid3_orig.th_centers(it);

% bin into new grid
Ny_new = numel(grid3_new.y_centers);
Nx_new = numel(grid3_new.x_centers);
Nt_new = numel(grid3_new.th_centers);

ix_new = max(1, min(Nx_new, floor((x_c  - grid3_new.x_edges(1))  / grid3_new.dx)  + 1));
iy_new = max(1, min(Ny_new, floor((y_c  - grid3_new.y_edges(1))  / grid3_new.dy)  + 1));
it_new = max(1, min(Nt_new, floor((th_c - grid3_new.th_edges(1)) / grid3_new.dtheta) + 1));

% Check Keep mask for new grid
if isfield(grid3_new, 'Keep') && ~isempty(grid3_new.Keep)
    inKeep = grid3_new.Keep(sub2ind(size(grid3_new.Keep), iy_new, ix_new));
else
    inKeep = true(n, 1);
end

% New voxel linear index
vid_new = sub2ind([Ny_new, Nx_new, Nt_new], iy_new, ix_new, it_new);

% Dedup: keep first hit (minimum |t|) per (iSeed, iHead, new_voxel)
iS  = double(rows_in.iSeed(1:n));
iH  = double(rows_in.iHead(1:n));
t_s = double(rows_in.t(1:n));

% sort by abs(t) ascending so first occurrence = smallest t
[~, sortOrd] = sort(abs(t_s));

% group key = (iSeed, iHead, vid_new) — use a combined hash
% iSeed/iHead are uint16 (0..65535), vid_new can be large
key = int64(iS) * int64(2^32) + int64(iH) * int64(2^20) + int64(vid_new);

seen = containers.Map('KeyType','int64','ValueType','logical');
keepMask = false(n, 1);
for ri = 1:n
    r = sortOrd(ri);
    if ~inKeep(r), continue; end
    k = key(r);
    if ~isKey(seen, k)
        seen(k) = true;
        keepMask(r) = true;
    end
end

idx = find(keepMask);
nKeep = numel(idx);

rows_out = rs3_rows_empty(nKeep);
rows_out.n = uint32(nKeep);
if nKeep == 0, return; end

rows_out.iSeed(1:nKeep)    = rows_in.iSeed(idx);
rows_out.iHead(1:nKeep)    = rows_in.iHead(idx);
rows_out.leg(1:nKeep)      = rows_in.leg(idx);
rows_out.halfFlag(1:nKeep) = rows_in.halfFlag(idx);
rows_out.t(1:nKeep)        = rows_in.t(idx);
rows_out.ix(1:nKeep)       = uint16(ix_new(idx));
rows_out.iy(1:nKeep)       = uint16(iy_new(idx));
rows_out.it(1:nKeep)       = uint16(it_new(idx));
end

% =========================================================================
function v = local_cfg_get(cfg, path, def)
v = def;
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        if ~isstruct(cur) || ~isfield(cur, parts{i}), return; end
        cur = cur.(parts{i});
    end
    if ~isempty(cur), v = cur; end
catch, end
end
