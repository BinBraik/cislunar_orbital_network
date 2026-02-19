function d = rs3_circ_diff(a, b)
%RS3_CIRC_DIFF  Smallest signed circular difference a-b in radians.
% Result is in [-pi, pi), consistent with rs3_wrapToPi convention.

d = rs3_wrapToPi(a - b);

end
