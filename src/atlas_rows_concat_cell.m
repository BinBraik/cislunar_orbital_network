function r = atlas_rows_concat_cell(cellArr)
%RS3_ROWS_CONCAT_CELL  Concatenate cell array of packed row structs.
%
% Equivalent to vertcat(cellArr{:}) for the old double-matrix format.
% Pre-allocates output for efficiency.

if isempty(cellArr)
    r = atlas_rows_empty();
    return;
end

% Compute total row count
ns = zeros(numel(cellArr), 1, 'uint32');
for i = 1:numel(cellArr)
    if ~isempty(cellArr{i}) && isstruct(cellArr{i})
        ns(i) = cellArr{i}.n;
    end
end
total = sum(ns);

if total == 0
    r = atlas_rows_empty();
    return;
end

% Pre-allocate
r = atlas_rows_empty(double(total));
r.n = total;

pos = uint32(0);
for i = 1:numel(cellArr)
    ni = ns(i);
    if ni == 0, continue; end
    c = cellArr{i};
    idx = pos+1 : pos+ni;
    r.iSeed(idx)    = c.iSeed(1:ni);
    r.iHead(idx)    = c.iHead(1:ni);
    r.leg(idx)      = c.leg(1:ni);
    r.halfFlag(idx) = c.halfFlag(1:ni);
    r.t(idx)        = c.t(1:ni);
    r.ix(idx)       = c.ix(1:ni);
    r.iy(idx)       = c.iy(1:ni);
    r.it(idx)       = c.it(1:ni);
    pos = pos + ni;
end
end
