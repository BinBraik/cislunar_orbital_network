function V = rs4_overlap_extract_voxel_info(SA, SB, O, cfg)
%RS4_OVERLAP_EXTRACT_VOXEL_INFO  Build voxel-wise overlap candidate metadata.
%
% For each overlap voxel, extract all contributing rows from:
%   - A.FRS_full (rows from SA)
%   - B.BRS_full (rows from SB)
% and derive per-row seed/heading/delta/DV-turn information so voxel-level
% ranking can be done later without rerunning propagation.
%
% Output
%   V : struct with fields
%       .summary   : scalar metadata
%       .voxels    : [K x 1] struct array, one per overlap voxel
%
% Per-voxel fields include:
%   id, ix,iy,it, x,y,th, countA,countB,
%   nRowsA, nRowsB,
%   uniqueSeedsA, uniqueSeedsB,
%   uniqueHeadsA, uniqueHeadsB,
%   uniqueSeedHeadPairsA, uniqueSeedHeadPairsB,
%   A, B (detailed row-level metadata structs)

if nargin < 4
    cfg = struct();
end

assert(isstruct(O) && isfield(O,'ids'), 'O must be overlap struct from rs4_overlap_pair');

grid3 = SA.grid3;
Ny = numel(grid3.y_centers);
Nx = numel(grid3.x_centers);
Nt = numel(grid3.th_centers);

VU_mps = local_cfg_get(cfg, 'units.VU_mps', 1.0);
TU_days = local_cfg_get(cfg, 'units.TU_days', 1.0);
usePar = local_cfg_get(cfg, 'rs4.extract.parallel', false);

if usePar && isempty(gcp('nocreate'))
    warning('[rs4] cfg.rs4.extract.parallel=true but no parpool exists. Running serial extraction.');
    usePar = false;
end

% ---------- reconstruct full sets exactly as in overlap stage ----------
rowsA_F_u = SA.Step4.rows_FRS_upper;
rowsA_B_u = SA.Step4.rows_BRS_upper;
rowsA_F_l = rs3_rows_mirror_lower(rowsA_B_u, grid3, 1);
rowsA_F = local_rows_cat(rowsA_F_u, rowsA_F_l);

rowsB_B_u = SB.Step4.rows_BRS_upper;
rowsB_F_u = SB.Step4.rows_FRS_upper;
rowsB_B_l = rs3_rows_mirror_lower(rowsB_F_u, grid3, 2);
rowsB_B = local_rows_cat(rowsB_B_u, rowsB_B_l);

idsA = rs3_rows_to_voxelid(rowsA_F, grid3);
idsB = rs3_rows_to_voxelid(rowsB_B, grid3);

[uA, idxCellsA] = local_group_indices(idsA);
[uB, idxCellsB] = local_group_indices(idsB);

idsO = O.ids(:);
K = numel(idsO);

[tfA, locA] = ismember(idsO, uA);
[tfB, locB] = ismember(idsO, uB);

cacheA = local_build_side_cache(SA);
cacheB = local_build_side_cache(SB);

voxels = repmat(local_empty_voxel(), K, 1);

if usePar
    parfor k = 1:K
        idxA_k = local_group_fetch(tfA(k), locA(k), idxCellsA);
        idxB_k = local_group_fetch(tfB(k), locB(k), idxCellsB);
        voxels(k) = local_build_voxel_entry(k, idsO, O, idxA_k, idxB_k, rowsA_F, rowsB_B, ...
            cacheA, cacheB, VU_mps, TU_days, grid3, Ny, Nx, Nt);
    end
else
    for k = 1:K
        idxA_k = local_group_fetch(tfA(k), locA(k), idxCellsA);
        idxB_k = local_group_fetch(tfB(k), locB(k), idxCellsB);
        voxels(k) = local_build_voxel_entry(k, idsO, O, idxA_k, idxB_k, rowsA_F, rowsB_B, ...
            cacheA, cacheB, VU_mps, TU_days, grid3, Ny, Nx, Nt);
    end
end

summary = struct();
summary.nVoxels = K;
summary.familyA = SA.name;
summary.familyB = SB.name;
summary.grid = struct('Nx',Nx,'Ny',Ny,'Nt',Nt,'dx',grid3.dx,'dy',grid3.dy,'dtheta',grid3.dtheta);
summary.generated = datestr(now, 31);
summary.note = ['Voxel-wise overlap candidate extraction (seeds/headings/delta/DV-turn/time) ' ...
                'for post-overlap ranking workflows.'];
summary.units = struct('dv_turn','m/s','tof','days','dv_turn_nd','nondimensional','tof_nd','nondimensional');
summary.VU_mps = VU_mps;
summary.TU_days = TU_days;
if isfield(cfg,'fan') && isfield(cfg.fan,'DV_cap_nd')
    summary.fan_DV_cap_nd = cfg.fan.DV_cap_nd;
end

V = struct();
V.summary = summary;
V.voxels = voxels;

fprintf('[rs4] extracted voxel metadata: %d overlap voxels\n', K);

end

% -------------------------------------------------------------------------
function info = local_side_info(M, cache, sideSign, VU_mps, TU_days)
% M columns: [iSeed iHead leg t ix iy it halfFlag]

info = struct();
if isempty(M)
    info.nRows = 0;
    info.uniqueSeeds = 0;
    info.uniqueHeads = 0;
    info.uniqueSeedHeadPairs = 0;
    info.iSeed = zeros(0,1);
    info.iHead = zeros(0,1);
    info.leg = zeros(0,1);
    info.t = zeros(0,1);
    info.ix = zeros(0,1);
    info.iy = zeros(0,1);
    info.it = zeros(0,1);
    info.halfFlag = zeros(0,1);
    info.delta = zeros(0,1);
    info.seed_x = zeros(0,1);
    info.seed_y = zeros(0,1);
    info.seed_th = zeros(0,1);
    info.heading_th = zeros(0,1);
    info.v0 = zeros(0,1);
    info.dv_turn = zeros(0,1);
    info.t_days = zeros(0,1);
    info.dv_turn_mps = zeros(0,1);
    info.t_min = NaN; info.t_max = NaN; info.t_mean = NaN;
    info.t_days_min = NaN; info.t_days_max = NaN; info.t_days_mean = NaN;
    info.dv_turn_min = NaN; info.dv_turn_max = NaN; info.dv_turn_mean = NaN;
    info.dv_turn_mps_min = NaN; info.dv_turn_mps_max = NaN; info.dv_turn_mps_mean = NaN;
    info.sideSign = sideSign;
    return;
end

iSeed = M(:,1);
iHead = M(:,2);
leg = M(:,3);
t = M(:,4);
ix = M(:,5);
iy = M(:,6);
it = M(:,7);
halfFlag = M(:,8);

n = size(M,1);

delta = local_lookup_delta_vec(cache.delta_lists, iSeed, iHead);

isUpper = halfFlag >= 0;
seed_x = zeros(n,1);
seed_y = zeros(n,1);
seed_th = zeros(n,1);
v0 = zeros(n,1);

if any(isUpper)
    iu = iSeed(isUpper);
    seed_x(isUpper) = cache.SeedsUpper(iu,1);
    seed_y(isUpper) = cache.SeedsUpper(iu,2);
    seed_th(isUpper) = cache.SeedsUpper(iu,3);
    v0(isUpper) = cache.v0_upper(iu);
end

if any(~isUpper)
    il = iSeed(~isUpper);
    seed_x(~isUpper) = cache.SeedsLower(il,1);
    seed_y(~isUpper) = cache.SeedsLower(il,2);
    seed_th(~isUpper) = cache.SeedsLower(il,3);
    v0(~isUpper) = cache.v0_lower(il);
end

heading_th = rs3_wrapToPi(seed_th + delta);
dv_turn = 2*v0.*sin(abs(delta)/2);

pairs = [iSeed, iHead];

info.nRows = n;
info.uniqueSeeds = numel(unique(iSeed));
info.uniqueHeads = numel(unique(iHead));
info.uniqueSeedHeadPairs = size(unique(pairs,'rows'),1);

info.iSeed = iSeed;
info.iHead = iHead;
info.leg = leg;
info.t = t;
info.ix = ix;
info.iy = iy;
info.it = it;
info.halfFlag = halfFlag;

info.delta = delta;
info.seed_x = seed_x;
info.seed_y = seed_y;
info.seed_th = seed_th;
info.heading_th = heading_th;
info.v0 = v0;
info.dv_turn = dv_turn;
info.t_days = t * TU_days;
info.dv_turn_mps = dv_turn * VU_mps;

info.t_min = min(t); info.t_max = max(t); info.t_mean = mean(t);
info.t_days_min = min(info.t_days); info.t_days_max = max(info.t_days); info.t_days_mean = mean(info.t_days);
info.dv_turn_min = min(dv_turn); info.dv_turn_max = max(dv_turn); info.dv_turn_mean = mean(dv_turn);
info.dv_turn_mps_min = min(info.dv_turn_mps); info.dv_turn_mps_max = max(info.dv_turn_mps); info.dv_turn_mps_mean = mean(info.dv_turn_mps);
info.sideSign = sideSign;
end

% -------------------------------------------------------------------------
function d = local_lookup_delta(delta_lists, iSeed, iHead)
if iSeed < 1 || iSeed > numel(delta_lists)
    d = NaN;
    return;
end
v = delta_lists{iSeed};
if isempty(v) || iHead < 1 || iHead > numel(v)
    d = NaN;
    return;
end
d = double(v(iHead));
end

% -------------------------------------------------------------------------
function d = local_lookup_delta_vec(delta_lists, iSeed, iHead)
n = numel(iSeed);
d = NaN(n,1);
for i = 1:n
    d(i) = local_lookup_delta(delta_lists, iSeed(i), iHead(i));
end
end

% -------------------------------------------------------------------------
function cache = local_build_side_cache(S)
cache = struct();
cache.delta_lists = S.Step4.delta_lists;
cache.SeedsUpper = S.SeedsUpper;
cache.SeedsLower = S.SeedsLower;

potU = rs3_core_cr3bp_U_and_derivs(cache.SeedsUpper(:,1), cache.SeedsUpper(:,2), S.mu);
potL = rs3_core_cr3bp_U_and_derivs(cache.SeedsLower(:,1), cache.SeedsLower(:,2), S.mu);
cache.v0_upper = sqrt(max(2*potU.U - S.CJ, 0));
cache.v0_lower = sqrt(max(2*potL.U - S.CJ, 0));
end

% -------------------------------------------------------------------------
function [u, idxCells] = local_group_indices(ids)
[idsSorted, ord] = sort(ids(:));
[u, ~, g] = unique(idsSorted);
idxCells = accumarray(g, ord, [numel(u), 1], @(v){v(:)});
end

% -------------------------------------------------------------------------
function idx = local_group_fetch(tf, loc, idxCells)
if ~tf || loc < 1 || loc > numel(idxCells)
    idx = zeros(0,1);
    return;
end
idx = idxCells{loc};
end

% -------------------------------------------------------------------------
function v = local_build_voxel_entry(k, idsO, O, idxA, idxB, rowsA_F, rowsB_B, ...
    cacheA, cacheB, VU_mps, TU_days, grid3, Ny, Nx, Nt)

id = idsO(k);
[iy, ix, it] = ind2sub([Ny, Nx, Nt], id);

rowsA_k = rs3_rows_subset(rowsA_F, idxA);
rowsB_k = rs3_rows_subset(rowsB_B, idxB);

MA = rs3_rows_to_matrix(rowsA_k);
MB = rs3_rows_to_matrix(rowsB_k);

Ainfo = local_side_info(MA, cacheA, +1, VU_mps, TU_days);
Binfo = local_side_info(MB, cacheB, -1, VU_mps, TU_days);

v = local_empty_voxel();
v.id = id;
v.ix = ix; v.iy = iy; v.it = it;
v.x = grid3.x_centers(ix);
v.y = grid3.y_centers(iy);
v.th = grid3.th_centers(it);

if isfield(O,'countA') && numel(O.countA) >= k
    v.countA = O.countA(k);
else
    v.countA = numel(idxA);
end
if isfield(O,'countB') && numel(O.countB) >= k
    v.countB = O.countB(k);
else
    v.countB = numel(idxB);
end

v.nRowsA = size(MA,1);
v.nRowsB = size(MB,1);
v.uniqueSeedsA = Ainfo.uniqueSeeds;
v.uniqueSeedsB = Binfo.uniqueSeeds;
v.uniqueHeadsA = Ainfo.uniqueHeads;
v.uniqueHeadsB = Binfo.uniqueHeads;
v.uniqueSeedHeadPairsA = Ainfo.uniqueSeedHeadPairs;
v.uniqueSeedHeadPairsB = Binfo.uniqueSeedHeadPairs;
v.A = Ainfo;
v.B = Binfo;
end

% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
function v = local_empty_voxel()
v = struct('id',NaN, 'ix',NaN, 'iy',NaN, 'it',NaN, 'x',NaN, 'y',NaN, 'th',NaN, ...
    'countA',0, 'countB',0, 'nRowsA',0, 'nRowsB',0, ...
    'uniqueSeedsA',0, 'uniqueSeedsB',0, 'uniqueHeadsA',0, 'uniqueHeadsB',0, ...
    'uniqueSeedHeadPairsA',0, 'uniqueSeedHeadPairsB',0, ...
    'A',struct(), 'B',struct());
end


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
