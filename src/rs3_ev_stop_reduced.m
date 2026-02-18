function [value,isterminal,direction] = rs3_ev_stop_reduced(~, X, mu, RE, RM, Rdom, CJ, v2tol)
%RS3_EV_STOP_REDUCED  Event function for reduced model integration.
% Stops when:
%  - enters Earth/Moon (rE<RE or rM<RM),
%  - exits domain (r>Rdom),
%  - violates CJ-allowed region by dropping below v2tol (v^2 < v2tol).
x = X(1); y = X(2);
xE = -mu; yE = 0;
xM = 1-mu; yM = 0;

rE = hypot(x-xE, y-yE);
rM = hypot(x-xM, y-yM);
rD = hypot(x, y);

pot = rs3_core_cr3bp_U_and_derivs(x, y, mu);
v2  = 2*pot.U - CJ;

value = [rE-RE; rM-RM; Rdom-rD; v2 - v2tol];
isterminal = [1; 1; 1; 1];
direction  = [-1; -1; -1; -1];
end
