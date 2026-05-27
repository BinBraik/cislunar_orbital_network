function edges = atlas_symmetric_edges(R, d)
%RS3_SYMMETRIC_EDGES  Build edges symmetric about 0 and including 0 exactly.
%
% edges = atlas_symmetric_edges(R, d) returns a monotone increasing edge vector
% spanning [-R, +R] with nominal spacing d, guaranteed to include 0 exactly.
%
% Note: If R/d is not an integer, the last bin on each side may be smaller
% than d (same behavior as the Step-2 grid builder).

assert(isfinite(R) && R>0, 'R must be positive.');
assert(isfinite(d) && d>0, 'd must be positive.');

pos = 0:d:R;
if isempty(pos), pos = 0; end
if abs(pos(end) - R) > 100*eps
    pos = [pos, R];
end

% Mirror, removing duplicate 0
neg = -fliplr(pos);
edges = [neg(1:end-1), pos];

assert(all(diff(edges) > 0), 'Edge construction failed (non-monotonic).');
end
