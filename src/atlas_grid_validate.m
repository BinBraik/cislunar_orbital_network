function atlas_grid_validate(grid3, cfg)
%RS3_GRID_VALIDATE  Unit/sanity tests for symmetric grid and theta binning.
%
% Checks
%   1) y=0 is an edge and no voxel straddles it
%   2) y edges are symmetric
%   3) theta seam periodic handling: bins across seam are adjacent in circular sense
%   4) representative in-range points return finite indices

assert(isstruct(grid3), 'grid3 must be a struct');

% Required fields
req = {'x_edges','y_edges','th_edges','x_centers','y_centers','th_centers','Nx','Ny','Nth'};
for k = 1:numel(req)
    assert(isfield(grid3, req{k}), 'grid3 missing required field: %s', req{k});
end

% ---- y=0 edge enforcement ----
if cfg.grid.enforce_y0_edge
    y_edges = grid3.y_edges;
    assert(any(abs(y_edges) < 10*eps), 'Validation failed: y=0 is not an edge.');

    % Ensure no voxel straddles 0: i.e., edges include 0 and there is a split
    % such that one bin ends at 0 and the next starts at 0.
    i0 = find(abs(y_edges) < 10*eps, 1, 'first');
    assert(~isempty(i0), 'Validation failed: could not locate y=0 edge index.');
    assert(i0 > 1 && i0 < numel(y_edges), 'Validation failed: y=0 at boundary, not interior edge.');
end

% ---- y symmetry check ----
if cfg.grid.enforce_xy_symmetry
    y_edges = grid3.y_edges;
    y_flip = -fliplr(y_edges);
    assert(numel(y_edges) == numel(y_flip), 'Validation failed: y edge length mismatch');
    assert(max(abs(y_edges - y_flip)) < 1e-12, 'Validation failed: y edges are not symmetric about 0');
end

% ---- theta periodic seam test ----
% Goal: points near -pi and +pi should map to bins whose centers are adjacent
% in circular distance (not separated by ~2*pi).
eps_th = 1e-12;
thA = -pi + eps_th;
thB =  pi - eps_th;

[~,~,itA] = grid3.bin_xyth(0, 0, thA);
[~,~,itB] = grid3.bin_xyth(0, 0, thB);

assert(isfinite(itA) && isfinite(itB), 'Theta seam test failed: bin indices are NaN');

% Compute circular distance between the corresponding bin centers
cA = grid3.th_centers(itA);
cB = grid3.th_centers(itB);
dc = abs(circ_diff(cA, cB));

% Adjacent bins should have circular center distance about one bin width
% Allow a bit of slack for discretization and the center definition.
assert(dc <= 1.5*grid3.dtheta + 1e-10, ...
    'Theta seam test failed: bins across seam are not adjacent (dc=%.3g rad)', dc);

% ---- representative in-range points test ----
xs = [0, 0.1*grid3.Rdom, -0.1*grid3.Rdom];
ys = [0.25*grid3.Rdom, -0.25*grid3.Rdom, 0.01*grid3.Rdom];
ths = [0, 0.5, -2.0];

[ix,iy,it] = grid3.bin_xyth(xs, ys, ths);
assert(all(isfinite(ix)) && all(isfinite(iy)) && all(isfinite(it)), 'Binning test failed: got NaN for in-range points');
assert(all(ix >= 1 & ix <= grid3.Nx), 'Binning test failed: ix out of range');
assert(all(iy >= 1 & iy <= grid3.Ny), 'Binning test failed: iy out of range');
assert(all(it >= 1 & it <= grid3.Nth), 'Binning test failed: it out of range');

end
