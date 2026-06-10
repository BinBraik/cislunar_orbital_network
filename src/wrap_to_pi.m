function a = wrap_to_pi(a)
%WRAP_TO_PI  Wrap angle(s) to [-pi, pi).
a = mod(a + pi, 2*pi) - pi;
end
