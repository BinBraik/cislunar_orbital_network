function Xdot = cr3bp_full_ode(~, X, mu)
%CR3BP_FULL_ODE  Full planar CR3BP model in (x,y,xd,yd).
% Reused from baseline nested full4d_cr3bp_model.
    x = X(1); y = X(2); xd = X(3); yd = X(4);
    pot = cr3bp_potential(x,y,mu);
    xdd =  2*yd + pot.Ux;
    ydd = -2*xd + pot.Uy;
    Xdot = [xd; yd; xdd; ydd];
end

function [value,isterminal,direction] = ev_stop_reduced(~, X, mu, RE, RM, Rdom, CJ, v2tol)
    x = X(1); y = X(2);
    xE = -mu; yE = 0;
    xM = 1-mu; yM = 0;

    rE = hypot(x-xE, y-yE);
    rM = hypot(x-xM, y-yM);
    rD = hypot(x, y);

    pot = cr3bp_potential(x,y,mu);
    v2  = 2*pot.U - CJ;

    value = [rE-RE; rM-RM; Rdom-rD; v2 - v2tol];
    isterminal = [1; 1; 1; 1];
    direction  = [-1; -1; -1; -1];
end

function [value,isterminal,direction] = ev_stop_full4d(~, X, mu, RE, RM, Rdom)
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
