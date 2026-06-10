function d = circ_diff(a, b)
%CIRC_DIFF  Smallest signed circular difference a-b in radians.
% Result is in [-pi, pi), consistent with wrap_to_pi convention.

d = wrap_to_pi(a - b);

end
