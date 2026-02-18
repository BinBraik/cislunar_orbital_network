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

idsO = O.ids(:);
K = numel(idsO);

voxels = repmat(local_empty_voxel(), K, 1);

for k = 1:K
    id = idsO(k);
    [iy, ix, it] = ind2sub([Ny, Nx, Nt], id);

    idxA = find(idsA == id);
    idxB = find(idsB == id);

    rowsA_k = rs3_rows_subset(rowsA_F, idxA);
    rowsB_k = rs3_rows_subset(rowsB_B, idxB);

    MA = rs3_rows_to_matrix(rowsA_k);
    MB = rs3_rows_to_matrix(rowsB_k);

    Ainfo = local_side_info(MA, SA, +1);
    Binfo = local_side_info(MB, SB, -1);

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

    voxels(k) = v;
end

summary = struct();
summary.nVoxels = K;
summary.familyA = SA.name;
summary.familyB = SB.name;
summary.grid = struct('Nx',Nx,'Ny',Ny,'Nt',Nt,'dx',grid3.dx,'dy',grid3.dy,'dtheta',grid3.dtheta);
summary.generated = datestr(now, 31);
summary.note = ['Voxel-wise overlap candidate extraction (seeds/headings/delta/DV-turn/time) ' ...
                'for post-overlap ranking workflows.'];
if isfield(cfg,'fan') && isfield(cfg.fan,'DV_cap_nd')
    summary.fan_DV_cap_nd = cfg.fan.DV_cap_nd;
end

V = struct();
V.summary = summary;
V.voxels = voxels;

fprintf('[rs4] extracted voxel metadata: %d overlap voxels\n', K);

end

% -------------------------------------------------------------------------
function info = local_side_info(M, S, sideSign)
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
    info.t_min = NaN; info.t_max = NaN; info.t_mean = NaN;
    info.dv_turn_min = NaN; info.dv_turn_max = NaN; info.dv_turn_mean = NaN;
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
delta = zeros(n,1);
seed_x = zeros(n,1);
seed_y = zeros(n,1);
seed_th = zeros(n,1);
heading_th = zeros(n,1);
v0 = zeros(n,1);
dv_turn = zeros(n,1);

for i = 1:n
    sidx = iSeed(i);
    hidx = iHead(i);

    delta(i) = local_lookup_delta(S.Step4.delta_lists, sidx, hidx);

    if halfFlag(i) >= 0
        seed = S.SeedsUpper(sidx,:);
    else
        seed = S.SeedsLower(sidx,:);
    end

    seed_x(i) = seed(1);
    seed_y(i) = seed(2);
    seed_th(i) = seed(3);
    heading_th(i) = rs3_wrapToPi(seed_th(i) + delta(i));

    pot = rs3_core_cr3bp_U_and_derivs(seed_x(i), seed_y(i), S.mu);
    v0(i) = sqrt(max(2*pot.U - S.CJ, 0));
    dv_turn(i) = 2*v0(i)*sin(abs(delta(i))/2);
end

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

info.t_min = min(t); info.t_max = max(t); info.t_mean = mean(t);
info.dv_turn_min = min(dv_turn); info.dv_turn_max = max(dv_turn); info.dv_turn_mean = mean(dv_turn);
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
