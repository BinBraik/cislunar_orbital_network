function r = atlas_rows_vcat(a, b)
%RS3_ROWS_VCAT  Concatenate two packed row structs.

na = double(a.n);
nb = double(b.n);
total = na + nb;

if total == 0
    r = atlas_rows_empty();
    return;
end

r = atlas_rows_empty(total);
r.n = uint32(total);

if na > 0
    r.iSeed(1:na)    = a.iSeed(1:na);
    r.iHead(1:na)    = a.iHead(1:na);
    r.leg(1:na)      = a.leg(1:na);
    r.halfFlag(1:na) = a.halfFlag(1:na);
    r.t(1:na)        = a.t(1:na);
    r.ix(1:na)       = a.ix(1:na);
    r.iy(1:na)       = a.iy(1:na);
    r.it(1:na)       = a.it(1:na);
end

if nb > 0
    idx = na+1:total;
    r.iSeed(idx)    = b.iSeed(1:nb);
    r.iHead(idx)    = b.iHead(1:nb);
    r.leg(idx)      = b.leg(1:nb);
    r.halfFlag(idx) = b.halfFlag(1:nb);
    r.t(idx)        = b.t(1:nb);
    r.ix(idx)       = b.ix(1:nb);
    r.iy(idx)       = b.iy(1:nb);
    r.it(idx)       = b.it(1:nb);
end
end
