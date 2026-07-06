function [Ax,Ay,Ath, Bx,By,Bth, Ox,Oy,Oth] = overlap_extract_xy_sets(SA, SB, O)
%OVERLAP_EXTRACT_XY_SETS  FRS-only, BRS-only, and overlap voxel centres.
%
% Shared extraction logic used by overlap_visualize_dark and
% overlap_visualize_gif (factored out of overlap_visualize.m so both new
% presentation outputs stay in sync with the standard overlap figure).
%
% Outputs:
%   (Ax,Ay,Ath) — A.FRS voxels not in the overlap
%   (Bx,By,Bth) — B.BRS voxels not in the overlap
%   (Ox,Oy,Oth) — overlap voxels

grid3 = SA.grid3;
Nx = numel(grid3.x_centers);
Ny = numel(grid3.y_centers);
Nt = numel(grid3.th_centers);

rowsA_F   = local_rows_cat(SA.Step4.rows_FRS_upper, SA.Step4.rows_FRS_lower);
rowsB_B_u = atlas_rows_mirror_lower(SB.Step4.rows_FRS_lower, grid3, 2);
rowsB_B_l = atlas_rows_mirror_lower(SB.Step4.rows_FRS_upper, grid3, 2);
rowsB_B   = local_rows_cat(rowsB_B_u, rowsB_B_l);

idsA = unique(local_rows_to_vid(rowsA_F, Ny, Nx, Nt));
idsB = unique(local_rows_to_vid(rowsB_B, Ny, Nx, Nt));
idsO = O.ids(:);

idsA_only = setdiff(idsA, idsO);
idsB_only = setdiff(idsB, idsO);

[Ax, Ay, Ath] = local_ids_to_centers(idsA_only, grid3, Ny, Nx, Nt);
[Bx, By, Bth] = local_ids_to_centers(idsB_only, grid3, Ny, Nx, Nt);
[Ox, Oy, Oth] = local_ids_to_centers(idsO,      grid3, Ny, Nx, Nt);

end

% ================= helpers (copied from overlap_visualize.m) =================

function ids = local_rows_to_vid(rows, Ny, Nx, Nt)
if isempty(rows), ids = zeros(0,1); return; end
if isstruct(rows)
    n = double(rows.n);
    if n==0, ids = zeros(0,1); return; end
    iy = double(rows.iy(1:n));
    ix = double(rows.ix(1:n));
    it = double(rows.it(1:n));
else
    ix = double(rows(:,5)); iy = double(rows(:,6)); it = double(rows(:,7));
end
ix = max(1, min(Nx, ix));
iy = max(1, min(Ny, iy));
it = max(1, min(Nt, it));
ids = sub2ind([Ny, Nx, Nt], iy, ix, it);
end

function [x,y,th] = local_ids_to_centers(ids, grid3, Ny, Nx, Nt)
if isempty(ids)
    x = []; y = []; th = [];
    return;
end
[iy, ix, it] = ind2sub([Ny, Nx, Nt], ids);
x  = grid3.x_centers(ix);
y  = grid3.y_centers(iy);
th = grid3.th_centers(it);
end

function r = local_rows_cat(a, b)
if isstruct(a) && isstruct(b)
    nA = double(a.n); nB = double(b.n);
    if nA==0, r=b; return; end
    if nB==0, r=a; return; end
    r = atlas_rows_empty(nA+nB);
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
