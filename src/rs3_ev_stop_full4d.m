function [value,isterminal,direction] = rs3_ev_stop_full4d(~, X, mu, RE, RM, Rdom)
%RS3_EV_STOP_FULL4D  Event function for full 4D model integration.
% Stops when:
%  - enters Earth/Moon (rE<RE or rM<RM),
%  - exits domain (r>Rdom).
x = X(1); y = X(2);
xE = -mu; yE = 0;
xM = 1-mu; yM = 0;

rE = hypot(x-xE, y-yE);
rM = hypot(x-xM, y-yM);
rD = hypot(x, y);

value = [rE-RE; rM-RM; Rdom-rD];
isterminal = [1; 1; 1];
direction  = [-1; -1; -1];
end
