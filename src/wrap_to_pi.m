function a = wrap_to_pi(a)
%RS3_WRAPTOPI  Wrap angle(s) to [-pi, pi).
a = mod(a + pi, 2*pi) - pi;
end
