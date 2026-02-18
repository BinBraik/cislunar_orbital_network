function O = rs4_overlap_pair(SA, SB, cfg)
%RS4_OVERLAP_PAIR  Overlap voxels between A.FRS_full and B.BRS_full,
% with conservative filtering:
%   - forbidden region removed via KeepA & KeepB
%   - buffers around primaries: r > (1+bufFrac)*R
%
% Output O:
%   ids, ix,iy,it, x,y,th, countA,countB, Nx,Ny,Nt

if nargin < 3, cfg = struct(); end

% ---------------- grid dims (robust) ----------------
grid3 = SA.grid3;

Nx = numel(grid3.x_centers);
Ny = numel(grid3.y_centers);
Nt = numel(grid3.th_centers);

gridB = SB.grid3;
assert(Nx==numel(gridB.x_centers) && Ny==numel(gridB.y_centers) && Nt==numel(gridB.th_centers), ...
    'A and B must share the same grid dims (x_centers/y_centers/th_centers).');

% ---------------- FULL sets ----------------
% A.FRS_full = FRS_upper + mirror(BRS_upper -> FRS_lower)
rowsA_F_u = SA.Step4.rows_FRS_upper;
rowsA_B_u = SA.Step4.rows_BRS_upper;
rowsA_F_l = rs3_rows_mirror_lower(rowsA_B_u, grid3, 1);
rowsA_F   = local_rows_cat(rowsA_F_u, rowsA_F_l);

% B.BRS_full = BRS_upper + mirror(FRS_upper -> BRS_lower)
rowsB_B_u = SB.Step4.rows_BRS_upper;
rowsB_F_u = SB.Step4.rows_FRS_upper;
rowsB_B_l = rs3_rows_mirror_lower(rowsB_F_u, grid3, 2);
rowsB_B   = local_rows_cat(rowsB_B_u, rowsB_B_l);

% ---------------- voxel IDs ----------------
idsA = local_rows_to_vid(rowsA_F, Ny, Nx, Nt);
idsB = local_rows_to_vid(rowsB_B, Ny, Nx, Nt);

idsO = intersect(unique(idsA), unique(idsB));
idsO = idsO(:);  % column

O = struct('ids',idsO,'ix',[],'iy',[],'it',[],'x',[],'y',[],'th',[], ...
           'countA',[],'countB',[],'Nx',Nx,'Ny',Ny,'Nt',Nt);

if isempty(idsO)
    fprintf('[rs4] overlap voxels: 0\n');
    return;
end

% ---------------- unpack overlap IDs ----------------
[iy, ix, it] = ind2sub([Ny, Nx, Nt], idsO);
ix = ix(:); iy = iy(:); it = it(:);  % force column

% ---------------- FILTER: Keep + buffered primaries ----------------
bufFrac = 0.05;
if isfield(cfg,'overlap') && isfield(cfg.overlap,'primary_buffer_frac') && ~isempty(cfg.overlap.primary_buffer_frac)
    bufFrac = cfg.overlap.primary_buffer_frac;
end

% Keep masks must be [Ny,Nx]
keepA = local_get_keep_xy(SA, Ny, Nx);
keepB = local_get_keep_xy(SB, Ny, Nx);
keepXY = keepA & keepB;  % conservative for CJ mismatch

% radii + mu (nd)
if ~(isfield(cfg,'sys') && isfield(cfg.sys,'RE_nd') && isfield(cfg.sys,'RM_nd'))
    error('cfg.sys.RE_nd and cfg.sys.RM_nd are required for primary buffer filtering.');
end
RE = cfg.sys.RE_nd;
RM = cfg.sys.RM_nd;

mu = SA.mu;  

x = grid3.x_centers(ix);
y = grid3.y_centers(iy);

rE = hypot(x + mu, y);           % Earth at (-mu,0)
rM = hypot(x - (1-mu), y);       % Moon  at (1-mu,0)

okKeep  = keepXY(sub2ind([Ny, Nx], iy, ix));
okEarth = rE > (1 + bufFrac)*RE;
okMoon  = rM > (1 + bufFrac)*RM;

% CRITICAL: force ALL to column vectors to avoid implicit expansion -> KxK
okKeep  = okKeep(:);
okEarth = okEarth(:);
okMoon  = okMoon(:);

ok = okKeep & okEarth & okMoon;  % now Kx1
ok = ok(:);

% Apply filter safely
idsO = idsO(ok);
ix   = ix(ok);
iy   = iy(ok);
it   = it(ok);

O.ids = idsO;

if isempty(idsO)
    O.ix=[]; O.iy=[]; O.it=[]; O.x=[]; O.y=[]; O.th=[];
    O.countA=[]; O.countB=[];
    fprintf('[rs4] overlap voxels after Keep+buffer: 0\n');
    return;
end

% ---------------- centers ----------------
O.ix = ix; O.iy = iy; O.it = it;
O.x  = grid3.x_centers(ix);
O.y  = grid3.y_centers(iy);
O.th = grid3.th_centers(it);

% ---------------- counts per kept overlap voxel ----------------
O.countA = local_counts_in_bins(idsA, idsO);
O.countB = local_counts_in_bins(idsB, idsO);

fprintf('[rs4] overlap voxels after Keep+buffer: %d\n', numel(O.ids));

end

% =======================================================================
% Helpers
% =======================================================================

function keep = local_get_keep_xy(S, Ny, Nx)
assert(isfield(S,'grid3') && isfield(S.grid3,'Keep') && ~isempty(S.grid3.Keep), ...
    'Keep mask missing: expected S.grid3.Keep (built in Step 3).');

keep = logical(S.grid3.Keep);
sz = size(keep);

if isequal(sz, [Ny, Nx])
    return;
elseif isequal(sz, [Nx, Ny])
    keep = keep.';   % common swapped dims
    return;
else
    error('Keep mask size [%d,%d], expected [%d,%d] (or transposed).', sz(1), sz(2), Ny, Nx);
end
end

function ids = local_rows_to_vid(rows, Ny, Nx, Nt)
if isempty(rows), ids = zeros(0,1); return; end

if isstruct(rows)
    n = double(rows.n);
    if n==0, ids = zeros(0,1); return; end
    iy = double(rows.iy(1:n));
    ix = double(rows.ix(1:n));
    it = double(rows.it(1:n));
else
    ix = double(rows(:,5));
    iy = double(rows(:,6));
    it = double(rows(:,7));
end

% defensive clamp
ix = max(1, min(Nx, ix));
iy = max(1, min(Ny, iy));
it = max(1, min(Nt, it));

ids = sub2ind([Ny, Nx, Nt], iy, ix, it);
ids = ids(:);
end

function c = local_counts_in_bins(idsAll, idsBins)
if isempty(idsBins)
    c = zeros(0,1);
    return;
end
idsBins = idsBins(:);
[tf, loc] = ismember(idsAll, idsBins);
loc = loc(tf);
c = accumarray(loc, 1, [numel(idsBins), 1], @sum, 0);
end

function r = local_rows_cat(a, b)
% concat packed rows without rs3_rows_vcat dependency
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
