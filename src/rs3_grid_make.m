function grid3 = rs3_grid_make(cfg)
%RS3_GRID_MAKE  Create symmetric (x,y,theta) voxel grid.
%
% Requirements enforced:
% - y=0 is an EDGE (no voxel straddles y=0)
% - y edges symmetric about 0
% - theta periodic seam-safe binning using [-pi,pi) convention
%
% Notes
% - x symmetry is not required for the x-axis symmetry constraint, but we
%   still build x edges symmetrically for cleanliness.

assert(isfield(cfg,'grid') && isstruct(cfg.grid), 'cfg.grid struct is required');
Rdom = cfg.grid.Rdom;
dx   = cfg.grid.dx;
dy   = cfg.grid.dy;
dth  = cfg.grid.dtheta;

assert(isfinite(Rdom) && Rdom > 0, 'cfg.grid.Rdom must be > 0');
assert(isfinite(dx) && dx > 0, 'cfg.grid.dx must be > 0');
assert(isfinite(dy) && dy > 0, 'cfg.grid.dy must be > 0');
assert(isfinite(dth) && dth > 0 && dth <= 2*pi, 'cfg.grid.dtheta must be in (0,2pi]');

% ---- x edges (symmetric; include 0 as edge) ----
x_edges = rs3_symmetric_edges(Rdom, dx);

% ---- y edges (MUST have y=0 as edge) ----
y_edges = rs3_symmetric_edges(Rdom, dy);

% Ensure y=0 is exactly an edge (robust check)
if cfg.grid.enforce_y0_edge
    assert(any(abs(y_edges) < 10*eps), 'y=0 is not an edge (numerical failure)');
end

% ---- theta edges (periodic seam-safe) ----
% Build edges over [-pi, pi] and wrap angles to [-pi,pi) before discretize.
Nth = max(1, round(2*pi / dth));
% Recompute dtheta to fit exactly into 2*pi (avoids accumulation)
dth_eff = 2*pi / Nth;
th_edges = linspace(-pi, pi, Nth + 1);

% Centers
x_centers = 0.5 * (x_edges(1:end-1) + x_edges(2:end));
y_centers = 0.5 * (y_edges(1:end-1) + y_edges(2:end));
th_centers = 0.5 * (th_edges(1:end-1) + th_edges(2:end));

% Sizes
Nx = numel(x_centers);
Ny = numel(y_centers);

grid3 = struct();
grid3.Rdom = Rdom;
grid3.dx = dx;
grid3.dy = dy;
grid3.dtheta = dth_eff;

grid3.x_edges = x_edges;
grid3.y_edges = y_edges;
grid3.th_edges = th_edges;

grid3.x_centers = x_centers;
grid3.y_centers = y_centers;
grid3.th_centers = th_centers;
grid3.xc = x_centers; % alias for compatibility
grid3.yc = y_centers; % alias for compatibility

grid3.Nx = Nx;
grid3.Ny = Ny;
grid3.Nth = Nth;

% Placeholder Keep mask (2D). The CJ-allowed region will be implemented later.
grid3.Keep = true(Ny, Nx);

% Helper handles
% Bin (x,y,theta) to indices (ix,iy,it). Returns NaN where out of bounds.
grid3.bin_xyth = @(x,y,th) rs3_bin_xyth(x,y,th,grid3);

% Convert indices to linear voxel id for 3D storage
% Note: MATLAB sub2ind uses (row, col, page) => (iy, ix, it)
grid3.voxel_id = @(ix,iy,it) sub2ind([Ny, Nx, Nth], iy, ix, it);

end

% -------------------------------------------------------------------------
function edges = rs3_symmetric_edges(R, d)
% Build edges symmetric about 0 and including 0 exactly.
%
% We avoid requiring that R/d is an integer. If it isn't, the last bin on
% each side will be smaller than d.

pos = 0:d:R;
if isempty(pos)
    pos = 0;
end
if abs(pos(end) - R) > 100*eps
    pos = [pos, R];
end
% Mirror, removing duplicate 0
neg = -fliplr(pos);
edges = [neg(1:end-1), pos];

% Sanity: edges monotonic increasing
assert(all(diff(edges) > 0), 'Edge construction failed (non-monotonic).');

end

% -------------------------------------------------------------------------
function [ix, iy, it] = rs3_bin_xyth(x, y, th, grid3)
% Bin inputs to voxel indices. Out-of-bounds => NaN.

% Wrap theta to [-pi,pi)
thw = rs3_wrapToPi(th);

ix = discretize(x, grid3.x_edges);
iy = discretize(y, grid3.y_edges);
it = discretize(thw, grid3.th_edges);

end