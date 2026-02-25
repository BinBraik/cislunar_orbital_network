function V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg)
%RS4_OVERLAP_EXTRACT_VOXEL_INFO  Build voxel-wise overlap candidate metadata.
%
% Vectorized implementation: avoids per-voxel O(N_rows) search and struct
% array construction. Uses ismember + accumarray to compute per-voxel
% statistics in O(N_rows log K) total instead of O(N_rows * K).
%
% Output V is a flat-array struct (not a struct array):
%   V.ids                    [K,1]  overlap voxel linear indices
%   V.ix, V.iy, V.it         [K,1]  grid index triplets
%   V.x,  V.y,  V.th         [K,1]  voxel center coordinates
%   V.nRowsA, V.nRowsB       [K,1]  row counts per voxel
%   V.uniqueSeedsA/B         [K,1]  unique seed count per voxel
%   V.uniqueHeadsA/B         [K,1]  unique heading count per voxel
%   V.dv_turn_mps_min/max/mean_A/B  [K,1]  DV-turn statistics (m/s)
%   V.t_days_min/max/mean_A/B       [K,1]  time-of-flight statistics (days)
%   V.summary                       scalar metadata struct

if nargin < 4
    cfg = struct();
end

assert(isstruct(O) && isfield(O,'ids'), 'O must be overlap struct from rs4_overlap_pair');

grid3 = SA.grid3;
Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

VU_mps  = local_cfg_get(cfg, 'units.VU_mps',  1.0);
TU_days = local_cfg_get(cfg, 'units.TU_days', 1.0);
usePar = local_cfg_get(cfg, 'rs4.extract.parallel', false);

if usePar && isempty(gcp('nocreate'))
    warning('[rs4] cfg.rs4.extract.parallel=true but no parpool exists. Running serial extraction.');
    usePar = false;
end

% ---------- reconstruct full sets (BRS = R(FRS)) ----------------------
% A.FRS_full: upper + lower stored directly from integration
rowsA_F = local_rows_cat(SA.Step4.rows_FRS_upper, SA.Step4.rows_FRS_lower);

% B.BRS_full = R(B.FRS_full): BRS_upper=R(FRS_lower), BRS_lower=R(FRS_upper)
rowsB_B_u = rs3_rows_mirror_lower(SB.Step4.rows_FRS_lower, grid3, 2);
rowsB_B_l = rs3_rows_mirror_lower(SB.Step4.rows_FRS_upper, grid3, 2);
rowsB_B = local_rows_cat(rowsB_B_u, rowsB_B_l);

idsA = double(rs3_rows_to_voxelid(rowsA_F, grid3));
idsB = double(rs3_rows_to_voxelid(rowsB_B, grid3));

idsO = double(O.ids(:));
K = numel(idsO);

% ---------- voxel center coordinates (vectorized ind2sub) ----------
[iy_vec, ix_vec, it_vec] = ind2sub([Ny, Nx, Nt], idsO);
x_vec  = grid3.x_centers(ix_vec);
y_vec  = grid3.y_centers(iy_vec);
th_vec = grid3.th_centers(it_vec);

% ---------- map each row to its overlap voxel (O(N log K)) ----------
[inO_A, rankA] = ismember(idsA, idsO);
[inO_B, rankB] = ismember(idsB, idsO);

rows_inO_A = find(inO_A);
rows_inO_B = find(inO_B);
group_A = rankA(rows_inO_A);   % 1..K index into overlap voxels
group_B = rankB(rows_inO_B);

% ---------- compute per-row dv_turn for side A ----------
[dv_mps_A, t_days_A, iSeed_A, iHead_A] = ...
    local_side_stats(rowsA_F, rows_inO_A, SA, +1, VU_mps, TU_days);

% ---------- compute per-row dv_turn for side B ----------
[dv_mps_B, t_days_B, iSeed_B, iHead_B] = ...
    local_side_stats(rowsB_B, rows_inO_B, SB, -1, VU_mps, TU_days);

% ---------- accumarray: per-voxel statistics ----------
V.ids = idsO;
V.ix  = ix_vec(:);
V.iy  = iy_vec(:);
V.it  = it_vec(:);
V.x   = x_vec(:);
V.y   = y_vec(:);
V.th  = th_vec(:);

V.nRowsA = accumarray(group_A, 1, [K,1], @sum, 0);
V.nRowsB = accumarray(group_B, 1, [K,1], @sum, 0);

V.dv_turn_mps_min_A  = accumarray(group_A, dv_mps_A, [K,1], @min,  NaN);
V.dv_turn_mps_max_A  = accumarray(group_A, dv_mps_A, [K,1], @max,  NaN);
V.dv_turn_mps_mean_A = accumarray(group_A, dv_mps_A, [K,1], @mean, NaN);

V.dv_turn_mps_min_B  = accumarray(group_B, dv_mps_B, [K,1], @min,  NaN);
V.dv_turn_mps_max_B  = accumarray(group_B, dv_mps_B, [K,1], @max,  NaN);
V.dv_turn_mps_mean_B = accumarray(group_B, dv_mps_B, [K,1], @mean, NaN);

V.t_days_min_A  = accumarray(group_A, t_days_A, [K,1], @min,  NaN);
V.t_days_max_A  = accumarray(group_A, t_days_A, [K,1], @max,  NaN);
V.t_days_mean_A = accumarray(group_A, t_days_A, [K,1], @mean, NaN);

V.t_days_min_B  = accumarray(group_B, t_days_B, [K,1], @min,  NaN);
V.t_days_max_B  = accumarray(group_B, t_days_B, [K,1], @max,  NaN);
V.t_days_mean_B = accumarray(group_B, t_days_B, [K,1], @mean, NaN);

% unique counts (custom accumulator — slower but only called once)
V.uniqueSeedsA = accumarray(group_A, iSeed_A, [K,1], @(x) numel(unique(x)), 0);
V.uniqueHeadsA = accumarray(group_A, iHead_A, [K,1], @(x) numel(unique(x)), 0);
V.uniqueSeedsB = accumarray(group_B, iSeed_B, [K,1], @(x) numel(unique(x)), 0);
V.uniqueHeadsB = accumarray(group_B, iHead_B, [K,1], @(x) numel(unique(x)), 0);

% ---------- summary ----------
summary = struct();
summary.nVoxels  = K;
summary.familyA  = SA.name;
summary.familyB  = SB.name;
summary.grid     = struct('Nx',Nx,'Ny',Ny,'Nt',Nt,'dx',grid3.dx,'dy',grid3.dy,'dtheta',grid3.dtheta);
summary.generated = datestr(now, 31);
summary.note     = 'Flat-array voxel metadata (vectorized; use V.field(k) not V.voxels(k)).';
summary.units    = struct('dv_turn_mps','m/s','tof_days','days');
summary.VU_mps   = VU_mps;
summary.TU_days  = TU_days;
if isfield(cfg,'fan') && isfield(cfg.fan,'DV_cap_nd')
    summary.fan_DV_cap_nd = cfg.fan.DV_cap_nd;
end
V.summary = summary;

fprintf('[rs4] extracted voxel metadata: %d overlap voxels\n', K);
end

% =========================================================================
function [dv_mps, t_days, iSeed_out, iHead_out] = ...
        local_side_stats(rows, rowIdx, S, sideSign, VU_mps, TU_days) %#ok<INUSD>
% Compute per-row dv_turn and tof for all in-overlap rows of one side.

if isempty(rowIdx)
    dv_mps    = zeros(0,1);
    t_days    = zeros(0,1);
    iSeed_out = zeros(0,1);
    iHead_out = zeros(0,1);
    return;
end

iSeed = double(rows.iSeed(rowIdx));
iHead = double(rows.iHead(rowIdx));
t_nd  = double(rows.t(rowIdx));
hf    = double(rows.halfFlag(rowIdx));

n = numel(rowIdx);

% --- delta lookup (tight loop; just array indexing into cell/vector) ---
delta_lists = S.Step4.delta_lists;
delta = zeros(n, 1);
for i = 1:n
    s = iSeed(i);  h = iHead(i);
    if s >= 1 && s <= numel(delta_lists)
        v = delta_lists{s};
        if h >= 1 && h <= numel(v)
            delta(i) = double(v(h));
        end
    end
end

% --- seed position lookup (vectorized matrix row-indexing) ---
seed_x = zeros(n, 1);
seed_y = zeros(n, 1);
is_upper = hf >= 0;

idxU = find(is_upper);
idxL = find(~is_upper);

if ~isempty(idxU)
    su = iSeed(idxU);
    seed_x(idxU) = S.SeedsUpper(su, 1);
    seed_y(idxU) = S.SeedsUpper(su, 2);
end
if ~isempty(idxL) && isfield(S,'SeedsLower') && ~isempty(S.SeedsLower)
    sl = iSeed(idxL);
    seed_x(idxL) = S.SeedsLower(sl, 1);
    seed_y(idxL) = S.SeedsLower(sl, 2);
end

% --- vectorized potential and v0 (one call for all rows) ---
pot  = rs3_core_cr3bp_U_and_derivs(seed_x, seed_y, S.mu);
v0   = sqrt(max(2*pot.U - S.CJ, 0));

dv_nd  = 2 * v0 .* sin(abs(delta) / 2);
dv_mps = dv_nd * VU_mps;
t_days = abs(t_nd) * TU_days;   % abs: BRS stores t < 0 (backward propagation)

iSeed_out = iSeed;
iHead_out = iHead;
end

% =========================================================================
function r = local_rows_cat(a, b)
if isstruct(a) && isstruct(b)
    nA = double(a.n); nB = double(b.n);
    if nA==0, r=b; return; end
    if nB==0, r=a; return; end

    r = rs3_rows_empty(nA+nB);
    r.n = uint32(nA+nB);

    r.iSeed(1:nA)    = a.iSeed(1:nA);    r.iSeed(nA+1:end)    = b.iSeed(1:nB);
    r.iHead(1:nA)    = a.iHead(1:nA);    r.iHead(nA+1:end)    = b.iHead(1:nB);
    r.leg(1:nA)      = a.leg(1:nA);      r.leg(nA+1:end)      = b.leg(1:nB);
    r.halfFlag(1:nA) = a.halfFlag(1:nA); r.halfFlag(nA+1:end) = b.halfFlag(1:nB);
    r.t(1:nA)        = a.t(1:nA);        r.t(nA+1:end)        = b.t(1:nB);
    r.ix(1:nA)       = a.ix(1:nA);       r.ix(nA+1:end)       = b.ix(1:nB);
    r.iy(1:nA)       = a.iy(1:nA);       r.iy(nA+1:end)       = b.iy(1:nB);
    r.it(1:nA)       = a.it(1:nA);       r.it(nA+1:end)       = b.it(1:nB);
    return;
end

if isempty(a), r=b; return; end
if isempty(b), r=a; return; end
r = [a; b];
end

% =========================================================================
function v = local_cfg_get(cfg, path, defaultVal)
v = defaultVal;
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        k = parts{i};
        if ~isstruct(cur) || ~isfield(cur, k)
            return;
        end
        cur = cur.(k);
    end
    if ~isempty(cur)
        v = cur;
    end
catch
    v = defaultVal;
end
end
