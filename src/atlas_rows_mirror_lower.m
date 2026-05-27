function rL = atlas_rows_mirror_lower(rU, grid3, target_leg)
%ATLAS_ROWS_MIRROR_LOWER  On-the-fly symmetry completion of lower-half rows.
%
% CR3BP symmetry: (t,x,y,th) -> (-t, x, -y, wrapToPi(pi-th))
%
% Swap logic:
%   FRS_lower = mirror(BRS_upper) with leg=1, time flipped
%   BRS_lower = mirror(FRS_upper) with leg=2, time flipped
%
% Inputs
%   rU         : packed row struct (upper-half rows)
%   grid3      : grid struct
%   target_leg : 1 for FRS_lower (input must be BRS_upper)
%                2 for BRS_lower (input must be FRS_upper)
%
% Output
%   rL : packed row struct (lower-half rows)

n = double(rU.n);
if n == 0
    rL = atlas_rows_empty();
    return;
end

Ny = grid3.Ny;

rL = atlas_rows_empty(n);
rL.n        = rU.n;
rL.iSeed    = rU.iSeed(1:n);
rL.iHead    = rU.iHead(1:n);
rL.leg      = repmat(uint8(target_leg), n, 1);
rL.halfFlag = repmat(int8(-1), n, 1);
rL.t        = -rU.t(1:n);                              % time flip
rL.ix       = rU.ix(1:n);                              % x unchanged
rL.iy       = uint16(Ny - double(rU.iy(1:n)) + 1);    % iy mirror

% Theta mirror: wrapToPi(pi - th)
th = grid3.th_centers(double(rU.it(1:n)));
thm = wrap_to_pi(pi - th);
itm = discretize(thm, grid3.th_edges);

% Drop any NaN theta bins (shouldn't happen with well-formed grids)
bad = isnan(itm);
if any(bad)
    keep = ~bad;
    rL = atlas_rows_subset(rL, keep);
    itm = itm(keep);
end
rL.it = uint16(itm);
end
