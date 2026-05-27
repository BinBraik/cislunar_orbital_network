function r = atlas_rows_empty(nAlloc)
%RS3_ROWS_EMPTY  Create an empty (or pre-allocated) packed row struct.
%
% Packed row format (struct-of-arrays, ~14 bytes/row vs 64 for double):
%   iSeed    uint16
%   iHead    uint16
%   leg      uint8
%   halfFlag int8
%   t        single
%   ix       uint16
%   iy       uint16
%   it       uint16
%   n        uint32   (number of valid rows)

if nargin < 1 || isempty(nAlloc), nAlloc = 0; end
nAlloc = max(0, nAlloc);

r = struct();
r.n        = uint32(0);
r.iSeed    = zeros(nAlloc, 1, 'uint16');
r.iHead    = zeros(nAlloc, 1, 'uint16');
r.leg      = zeros(nAlloc, 1, 'uint8');
r.halfFlag = zeros(nAlloc, 1, 'int8');
r.t        = zeros(nAlloc, 1, 'single');
r.ix       = zeros(nAlloc, 1, 'uint16');
r.iy       = zeros(nAlloc, 1, 'uint16');
r.it       = zeros(nAlloc, 1, 'uint16');
end
