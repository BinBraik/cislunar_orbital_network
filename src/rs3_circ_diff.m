function d = rs3_circ_diff(a, b)
%RS3_CIRC_DIFF  Smallest signed circular difference a-b in radians.
% Result is in (-pi, pi].

d = rs3_wrapToPi(a - b);

% Map -pi to +pi for symmetry (optional, but helps some asserts)
mask = (d <= -pi);
if any(mask(:))
    d(mask) = d(mask) + 2*pi;
end

end
