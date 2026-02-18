function a = rs3_wrapToPi(a)
%RS3_WRAPTOPI  Wrap angle(s) to [-pi, pi).
a = mod(a + pi, 2*pi) - pi;
% Optional: keep +pi out of range to avoid discretize edge ambiguity
mask = (a == pi);
a(mask) = pi - 1e-15;
end
