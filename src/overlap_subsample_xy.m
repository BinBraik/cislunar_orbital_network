function [xs, ys, keepIdx] = overlap_subsample_xy(x, y, capN, seed)
%OVERLAP_SUBSAMPLE_XY  Randomly cap a point cloud for fast/lightweight rendering.
%
% Presentation scatters (dark overlap figure, process GIF) don't need every
% single voxel to read correctly — at typical marker sizes a random subsample
% is visually indistinguishable from the full cloud, but draws and exports
% far faster (this is the main lever against slow GIF rendering and
% oversized vector PDF/EPS files from huge point counts).
%
% Usage:
%   [xs, ys] = overlap_subsample_xy(x, y, capN, seed)
%
% If numel(x) <= capN, x and y are returned unchanged.

if nargin < 4 || isempty(seed), seed = 1; end

x = x(:); y = y(:);
n = numel(x);
if n <= capN
    xs = x; ys = y; keepIdx = (1:n)';
    return;
end

rs = RandStream('twister', 'Seed', seed);
keepIdx = sort(randperm(rs, n, capN))';
xs = x(keepIdx);
ys = y(keepIdx);
end
