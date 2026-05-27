function r = atlas_rows_subset(src, idx)
%RS3_ROWS_SUBSET  Extract a subset of packed rows by logical mask or indices.
%
% idx can be a logical vector or an integer index vector.

if ~isstruct(src)
    % Fallback for double matrix
    r = src(idx, :);
    return;
end

n = double(src.n);

% Apply indexing to 1:n portion
r = atlas_rows_empty();
r.iSeed    = src.iSeed(idx);
r.iHead    = src.iHead(idx);
r.leg      = src.leg(idx);
r.halfFlag = src.halfFlag(idx);
r.t        = src.t(idx);
r.ix       = src.ix(idx);
r.iy       = src.iy(idx);
r.it       = src.it(idx);
r.n        = uint32(numel(r.iSeed));
end
