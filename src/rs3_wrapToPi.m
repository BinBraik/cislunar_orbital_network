function a = rs3_wrapToPi(a)
%RS3_WRAPTOPI  Wrap angle(s) to [-pi, pi).
a = mod(a + pi, 2*pi) - pi;
end
