function rs3_hits_validate(rowsF, rowsB, grid3)
%RS3_HITS_VALIDATE  Sanity checks for voxel hit rows.
% Handles both packed struct and legacy [N x 8] double matrix formats.

if isstruct(rowsF)
    nF = double(rowsF.n);
    if nF > 0
        assert(all(rowsF.ix(1:nF) >= 1 & rowsF.ix(1:nF) <= grid3.Nx), 'FRS ix out of bounds');
        assert(all(rowsF.iy(1:nF) >= 1 & rowsF.iy(1:nF) <= grid3.Ny), 'FRS iy out of bounds');
        assert(all(rowsF.it(1:nF) >= 1 & rowsF.it(1:nF) <= grid3.Nth), 'FRS it out of bounds');
        assert(all(rowsF.leg(1:nF) == 1), 'FRS leg must be 1');
        assert(all(rowsF.halfFlag(1:nF) == 1), 'FRS upper halfFlag must be +1');
        assert(all(rowsF.t(1:nF) >= -1e-6), 'FRS times should be >=0 (forward)');
    end
else
    assert(size(rowsF,2)==8, 'rows_FRS_upper must be Nx8');
    if ~isempty(rowsF)
        assert(all(rowsF(:,5) >= 1 & rowsF(:,5) <= grid3.Nx), 'FRS ix out of bounds');
        assert(all(rowsF(:,6) >= 1 & rowsF(:,6) <= grid3.Ny), 'FRS iy out of bounds');
        assert(all(rowsF(:,7) >= 1 & rowsF(:,7) <= grid3.Nth), 'FRS it out of bounds');
        assert(all(rowsF(:,3) == 1), 'FRS leg must be 1');
        assert(all(rowsF(:,8) == 1), 'FRS upper halfFlag must be +1');
        assert(all(rowsF(:,4) >= -1e-12), 'FRS times should be >=0 (forward)');
    end
end

if isstruct(rowsB)
    nB = double(rowsB.n);
    if nB > 0
        assert(all(rowsB.ix(1:nB) >= 1 & rowsB.ix(1:nB) <= grid3.Nx), 'BRS ix out of bounds');
        assert(all(rowsB.iy(1:nB) >= 1 & rowsB.iy(1:nB) <= grid3.Ny), 'BRS iy out of bounds');
        assert(all(rowsB.it(1:nB) >= 1 & rowsB.it(1:nB) <= grid3.Nth), 'BRS it out of bounds');
        assert(all(rowsB.leg(1:nB) == 2), 'BRS leg must be 2');
        assert(all(rowsB.halfFlag(1:nB) == 1), 'BRS upper halfFlag must be +1');
        assert(all(rowsB.t(1:nB) <= 1e-6), 'BRS times should be <=0 (backward)');
    end
else
    assert(size(rowsB,2)==8, 'rows_BRS_upper must be Nx8');
    if ~isempty(rowsB)
        assert(all(rowsB(:,5) >= 1 & rowsB(:,5) <= grid3.Nx), 'BRS ix out of bounds');
        assert(all(rowsB(:,6) >= 1 & rowsB(:,6) <= grid3.Ny), 'BRS iy out of bounds');
        assert(all(rowsB(:,7) >= 1 & rowsB(:,7) <= grid3.Nth), 'BRS it out of bounds');
        assert(all(rowsB(:,3) == 2), 'BRS leg must be 2');
        assert(all(rowsB(:,8) == 1), 'BRS upper halfFlag must be +1');
        assert(all(rowsB(:,4) <= 1e-12), 'BRS times should be <=0 (backward)');
    end
end

end
