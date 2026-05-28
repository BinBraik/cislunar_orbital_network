function M = atlas_rows_to_matrix(r)
%ATLAS_ROWS_TO_MATRIX  Convert packed row struct to [N x 8] double matrix.
%
% Used when small subsets of rows are needed in double format (e.g.,
% per-voxel candidate extracts for Step 8 scoring).

if ~isstruct(r) 
    M = r;   % already a matrix
    return;
end

n = double(r.n);
if n == 0
    M = zeros(0, 8);
    return;
end

M = zeros(n, 8);
M(:,1) = double(r.iSeed(1:n));
M(:,2) = double(r.iHead(1:n));
M(:,3) = double(r.leg(1:n));
M(:,4) = double(r.t(1:n));
M(:,5) = double(r.ix(1:n));
M(:,6) = double(r.iy(1:n));
M(:,7) = double(r.it(1:n));
M(:,8) = double(r.halfFlag(1:n));
end
