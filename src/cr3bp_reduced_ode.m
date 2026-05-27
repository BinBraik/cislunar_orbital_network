function Xdot = cr3bp_reduced_ode(~,X,CJ,mu,hard_fail)
%RS3_CORE_REDUCED_CR3BP_MODEL  Reduced planar CR3BP model in (x,y,theta).
% Reused from baseline nested reduced_cr3bp_model.
    x  = X(1); y  = X(2); th = X(3);
    pot = cr3bp_potential(x,y,mu);
    v2  = 2*pot.U - CJ;
    if v2 <= 0
        if hard_fail
            error('v^2=2U-CJ <= 0 encountered (x=%.6g,y=%.6g,CJ=%.6g).', x, y, CJ);
        else
            Xdot = [0;0;0]; return;
        end
    end
    v = sqrt(v2);
    xdot = v*cos(th);
    ydot = v*sin(th);
    thetadot = -2 + (pot.Uy*xdot - pot.Ux*ydot) / v2;
    Xdot = [xdot; ydot; thetadot];
end

