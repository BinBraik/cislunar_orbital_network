function atlas_hits_validate(rowsF_u, rowsF_l, grid3)
%RS3_HITS_VALIDATE  Sanity checks for voxel hit rows.
% Handles both packed struct and legacy [N x 8] double matrix formats.
%
% rowsF_u : rows_FRS_upper (leg=1, halfFlag=+1, t>=0)
% rowsF_l : rows_FRS_lower (leg=1, halfFlag=-1, t>=0)

if isstruct(rowsF_u)
    nF = double(rowsF_u.n);
    if nF > 0
        assert(all(rowsF_u.ix(1:nF) >= 1 & rowsF_u.ix(1:nF) <= grid3.Nx), 'FRS_upper ix out of bounds');
        assert(all(rowsF_u.iy(1:nF) >= 1 & rowsF_u.iy(1:nF) <= grid3.Ny), 'FRS_upper iy out of bounds');
        assert(all(rowsF_u.it(1:nF) >= 1 & rowsF_u.it(1:nF) <= grid3.Nth), 'FRS_upper it out of bounds');
        assert(all(rowsF_u.leg(1:nF) == 1), 'FRS_upper leg must be 1');
        assert(all(rowsF_u.halfFlag(1:nF) == 1), 'FRS_upper halfFlag must be +1');
        assert(all(rowsF_u.t(1:nF) >= -1e-6), 'FRS_upper times should be >=0 (forward)');
    end
else
    assert(size(rowsF_u,2)==8, 'rows_FRS_upper must be Nx8');
    if ~isempty(rowsF_u)
        assert(all(rowsF_u(:,5) >= 1 & rowsF_u(:,5) <= grid3.Nx), 'FRS_upper ix out of bounds');
        assert(all(rowsF_u(:,6) >= 1 & rowsF_u(:,6) <= grid3.Ny), 'FRS_upper iy out of bounds');
        assert(all(rowsF_u(:,7) >= 1 & rowsF_u(:,7) <= grid3.Nth), 'FRS_upper it out of bounds');
        assert(all(rowsF_u(:,3) == 1), 'FRS_upper leg must be 1');
        assert(all(rowsF_u(:,8) == 1), 'FRS_upper halfFlag must be +1');
        assert(all(rowsF_u(:,4) >= -1e-12), 'FRS_upper times should be >=0 (forward)');
    end
end

if isstruct(rowsF_l)
    nL = double(rowsF_l.n);
    if nL > 0
        assert(all(rowsF_l.ix(1:nL) >= 1 & rowsF_l.ix(1:nL) <= grid3.Nx), 'FRS_lower ix out of bounds');
        assert(all(rowsF_l.iy(1:nL) >= 1 & rowsF_l.iy(1:nL) <= grid3.Ny), 'FRS_lower iy out of bounds');
        assert(all(rowsF_l.it(1:nL) >= 1 & rowsF_l.it(1:nL) <= grid3.Nth), 'FRS_lower it out of bounds');
        assert(all(rowsF_l.leg(1:nL) == 1), 'FRS_lower leg must be 1');
        assert(all(rowsF_l.halfFlag(1:nL) == -1), 'FRS_lower halfFlag must be -1');
        assert(all(rowsF_l.t(1:nL) >= -1e-6), 'FRS_lower times should be >=0 (forward)');
    end
else
    assert(size(rowsF_l,2)==8, 'rows_FRS_lower must be Nx8');
    if ~isempty(rowsF_l)
        assert(all(rowsF_l(:,5) >= 1 & rowsF_l(:,5) <= grid3.Nx), 'FRS_lower ix out of bounds');
        assert(all(rowsF_l(:,6) >= 1 & rowsF_l(:,6) <= grid3.Ny), 'FRS_lower iy out of bounds');
        assert(all(rowsF_l(:,7) >= 1 & rowsF_l(:,7) <= grid3.Nth), 'FRS_lower it out of bounds');
        assert(all(rowsF_l(:,3) == 1), 'FRS_lower leg must be 1');
        assert(all(rowsF_l(:,8) == -1), 'FRS_lower halfFlag must be -1');
        assert(all(rowsF_l(:,4) >= -1e-12), 'FRS_lower times should be >=0 (forward)');
    end
end

end
