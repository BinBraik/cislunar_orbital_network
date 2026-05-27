function [t, X] = cr3bp_integrate(X0red, tspan, CJ, mu, v2tol, optsR, opts4)
%RS3_CORE_INTEGRATE_REDUCED_WITH_FALLBACK
% Reduced integration with event stop; if stop due to v2tol, redo from same
% IC in full planar 4D.
%
% Inputs
%   X0red : [3x1] reduced IC [x;y;theta]
%   tspan : [t0 tf] (can be decreasing for backward integration)
%   CJ, mu, v2tol : system params
%   optsR : odeset options for reduced model (with ev_stop_reduced)
%   opts4 : odeset options for full 4D model (with ev_stop_full4d)
%
% Outputs
%   t : time vector (monotone with tspan direction)
%   X : trajectory states. If reduced integration succeeded, size is [n x 3].
%       If fallback triggered, size is [n x 4] full state [x y xd yd].

% Guard: MATLAB ODE solvers require tspan(1) ~= tspan(2)
if numel(tspan)==2 && tspan(1)==tspan(2)
    t = tspan(1);
    X = X0red(:).';
    return;
end

[tR, XR, ~, ~, ie] = ode113(@(tt,XX) cr3bp_reduced_ode(tt, XX, CJ, mu, false), ...
                            tspan, X0red, optsR);

ie_last = 0;
if ~isempty(ie)
    ie_last = ie(end);
end

% Baseline convention: event #4 corresponds to v2 - v2tol crossing
if ie_last == 4
    x  = X0red(1); y  = X0red(2); th = X0red(3);
    pot = cr3bp_potential(x, y, mu);
    v2  = max(2*pot.U - CJ, 0);
    v   = sqrt(v2);
    xd0 = v*cos(th);
    yd0 = v*sin(th);
    X0_4 = [x; y; xd0; yd0];

    [t, X] = ode113(@(tt,XX) cr3bp_full_ode(tt, XX, mu), ...
                    tspan, X0_4, opts4);
else
    t = tR;
    X = XR;
end
end
