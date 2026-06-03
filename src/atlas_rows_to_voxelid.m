function vid = atlas_rows_to_voxelid(rows, grid3)
%ATLAS_ROWS_TO_VOXELID  Convert hit rows to linear voxel IDs.
% Handles both packed struct and legacy [N x 8] double matrix.

if isstruct(rows)
    n = double(rows.n);
    if n == 0
        vid = zeros(0,1);
        return;
    end
    ix = double(rows.ix(1:n));
    iy = double(rows.iy(1:n));
    it = double(rows.it(1:n));
else
    if isempty(rows)
        vid = zeros(0,1);
        return;
    end
    ix = rows(:,5);
    iy = rows(:,6);
    it = rows(:,7);
end

Ny = grid3.Ny;
Nx = grid3.Nx;
Nth = grid3.Nth;
vid = sub2ind([Ny, Nx, Nth], iy, ix, it);
end
