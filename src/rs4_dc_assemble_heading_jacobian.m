function J = rs4_dc_assemble_heading_jacobian(A, B, cfg)
%RS4_DC_ASSEMBLE_HEADING_JACOBIAN  Assemble heading-residual Jacobian for two arcs.
% Preferred path uses STM/variational equations per arc; fallback uses finite differences.
%
% Inputs
%   A, B : structs with fields
%       .IC      [3x1] initial [x;y;th] in FRS frame
%       .tfinal  scalar final time
%       .CJ      Jacobi constant
%       .mu      CR3BP mass ratio
%       .mirror_heading logical; if true, output heading is pi-theta (BRS transform)
%   cfg  : config struct. Uses cfg.rs4.dc.jacobian_mode in {'stm','fd','auto'}.
%
% Output J (struct)
%   .mode_requested
%   .mode_used
%   .residual_heading_rad     wrapped heading residual in [-pi,pi)
%   .theta_A_rad              terminal heading for arc A (in requested frame)
%   .theta_B_rad              terminal heading for arc B (in requested frame)
%   .jacobian                 [1x2] = [d(res)/d(deltaA), d(res)/d(deltaB)]
%   .details                  diagnostics including per-arc sensitivities

modeReq = local_get_cfg(cfg, 'rs4.dc.jacobian_mode', 'auto');
if isstring(modeReq), modeReq = char(modeReq); end
modeReq = lower(strtrim(modeReq));
if ~ismember(modeReq, {'stm','fd','auto'})
    error('rs4:dc:badMode', 'cfg.rs4.dc.jacobian_mode must be ''stm'', ''fd'', or ''auto''.');
end

fdStep = local_get_cfg(cfg, 'rs4.dc.fd_step', 1e-6);
absTol = local_get_cfg(cfg, 'propag.absTol', 1e-9);
relTol = local_get_cfg(cfg, 'propag.relTol', 1e-9);

J = struct();
J.mode_requested = modeReq;
J.mode_used = '';
J.residual_heading_rad = NaN;
J.theta_A_rad = NaN;
J.theta_B_rad = NaN;
J.jacobian = [NaN NaN];
J.details = struct();

trySTM = ismember(modeReq, {'stm','auto'});
if trySTM
    [okSTM, outSTM] = local_try_stm(A, B, absTol, relTol);
    if okSTM
        J.mode_used = 'stm';
        J.theta_A_rad = outSTM.thetaA;
        J.theta_B_rad = outSTM.thetaB;
        J.residual_heading_rad = local_heading_residual(outSTM.thetaA, outSTM.thetaB);
        J.jacobian = [outSTM.dthA_dth0, -outSTM.dthB_dth0];
        J.details.stm = outSTM;
        fprintf('[rs4_dc] Jacobian mode used: stm\n');
        return;
    elseif strcmp(modeReq, 'stm')
        error('rs4:dc:stmFailed', 'STM Jacobian requested but unavailable/unstable: %s', outSTM.reason);
    end
end

[outFD] = local_fd(A, B, absTol, relTol, fdStep);
J.mode_used = 'fd';
J.theta_A_rad = outFD.thetaA;
J.theta_B_rad = outFD.thetaB;
J.residual_heading_rad = outFD.residual;
J.jacobian = outFD.jac;
J.details.fd = outFD;
if trySTM && isfield(J.details, 'stm') && isfield(J.details.stm, 'reason')
    fprintf('[rs4_dc] Jacobian mode used: fd (STM fallback: %s)\n', J.details.stm.reason);
else
    fprintf('[rs4_dc] Jacobian mode used: fd\n');
end

end

% -------------------------------------------------------------------------
function [ok, out] = local_try_stm(A, B, absTol, relTol)
out = struct('reason', 'unknown');
ok = false;

[okA, arcA] = local_integrate_arc_stm(A, absTol, relTol);
if ~okA
    out.reason = sprintf('Arc A STM failed: %s', arcA.reason);
    return;
end

[okB, arcB] = local_integrate_arc_stm(B, absTol, relTol);
if ~okB
    out.reason = sprintf('Arc B STM failed: %s', arcB.reason);
    return;
end

if ~(isfinite(arcA.dth_dth0) && isfinite(arcB.dth_dth0))
    out.reason = 'Non-finite STM sensitivities';
    return;
end

out.thetaA = arcA.theta_out;
out.thetaB = arcB.theta_out;
out.dthA_dth0 = arcA.dth_dth0;
out.dthB_dth0 = arcB.dth_dth0;
out.arcA = arcA;
out.arcB = arcB;
ok = true;
end

% -------------------------------------------------------------------------
function out = local_fd(A, B, absTol, relTol, h)
thA0 = local_arc_terminal_heading(A, absTol, relTol);
thB0 = local_arc_terminal_heading(B, absTol, relTol);
r0   = local_heading_residual(thA0, thB0);

A_p = A; A_m = A;
A_p.IC(3) = rs3_wrapToPi(A.IC(3) + h);
A_m.IC(3) = rs3_wrapToPi(A.IC(3) - h);
rAp = local_heading_residual(local_arc_terminal_heading(A_p, absTol, relTol), thB0);
rAm = local_heading_residual(local_arc_terminal_heading(A_m, absTol, relTol), thB0);
dr_dA = local_wrap_diff(rAp, rAm) / (2*h);

B_p = B; B_m = B;
B_p.IC(3) = rs3_wrapToPi(B.IC(3) + h);
B_m.IC(3) = rs3_wrapToPi(B.IC(3) - h);
rBp = local_heading_residual(thA0, local_arc_terminal_heading(B_p, absTol, relTol));
rBm = local_heading_residual(thA0, local_arc_terminal_heading(B_m, absTol, relTol));
dr_dB = local_wrap_diff(rBp, rBm) / (2*h);

out = struct();
out.thetaA = thA0;
out.thetaB = thB0;
out.residual = r0;
out.jac = [dr_dA, dr_dB];
out.fd_step = h;
end

% -------------------------------------------------------------------------
function [ok, arc] = local_integrate_arc_stm(arcIn, absTol, relTol)
arc = struct('reason','');
ok = false;
if arcIn.tfinal < 0
    arc.reason = 'tfinal must be >= 0';
    return;
end
IC = arcIn.IC(:);
if numel(IC) ~= 3
    arc.reason = 'IC must have 3 states';
    return;
end

Phi0 = eye(3);
Y0 = [IC; Phi0(:)];
opts = odeset('RelTol', relTol, 'AbsTol', absTol);

try
    [~, Y] = ode113(@(t,y) local_aug_rhs(t, y, arcIn.CJ, arcIn.mu), [0, arcIn.tfinal], Y0, opts);
catch ME
    arc.reason = ME.message;
    return;
end

Yf = Y(end,:).';
Xf = Yf(1:3);
Phi = reshape(Yf(4:end), 3, 3);
if any(~isfinite(Xf)) || any(~isfinite(Phi(:)))
    arc.reason = 'non-finite state/STM';
    return;
end

theta_out = Xf(3);
dth = Phi(3,3);

if isfield(arcIn, 'mirror_heading') && arcIn.mirror_heading
    theta_out = rs3_wrapToPi(pi - theta_out);
    dth = -dth;
else
    theta_out = rs3_wrapToPi(theta_out);
end

arc.theta_out = theta_out;
arc.dth_dth0 = dth;
ok = true;
end

% -------------------------------------------------------------------------
function dy = local_aug_rhs(~, y, CJ, mu)
X = y(1:3);
Phi = reshape(y(4:end), 3, 3);

pot = rs3_core_cr3bp_U_and_derivs(X(1), X(2), mu);
v2 = 2*pot.U - CJ;
if v2 <= 0
    dy = zeros(12,1);
    return;
end
v = sqrt(v2);
c = cos(X(3));
s = sin(X(3));

f = [v*c;
     v*s;
     -2 + (pot.Uy*c - pot.Ux*s)/v];

A = zeros(3,3);
A(1,1) = (pot.Ux/v) * c;
A(1,2) = (pot.Uy/v) * c;
A(1,3) = -v*s;
A(2,1) = (pot.Ux/v) * s;
A(2,2) = (pot.Uy/v) * s;
A(2,3) =  v*c;

g = pot.Uy*c - pot.Ux*s;
dgx = pot.Uxy*c - pot.Uxx*s;
dgy = pot.Uyy*c - pot.Uxy*s;
dgth = -pot.Uy*s - pot.Ux*c;

A(3,1) = dgx / v - g * pot.Ux / (v^3);
A(3,2) = dgy / v - g * pot.Uy / (v^3);
A(3,3) = dgth / v;

PhiDot = A * Phi;
dy = [f; PhiDot(:)];
end

% -------------------------------------------------------------------------
function thOut = local_arc_terminal_heading(arcIn, absTol, relTol)
opts = odeset('RelTol', relTol, 'AbsTol', absTol);
[~, X] = ode113(@(tt,XX) rs3_core_reduced_cr3bp_model(tt, XX, arcIn.CJ, arcIn.mu, false), ...
    [0, arcIn.tfinal], arcIn.IC(:), opts);
thOut = X(end,3);
if isfield(arcIn, 'mirror_heading') && arcIn.mirror_heading
    thOut = rs3_wrapToPi(pi - thOut);
else
    thOut = rs3_wrapToPi(thOut);
end
end

% -------------------------------------------------------------------------
function r = local_heading_residual(thA, thB)
% Wrap-safe heading residual using atan2(sin,cos).
d = thA - thB;
r = atan2(sin(d), cos(d));
end

function d = local_wrap_diff(a, b)
% Wrap-safe difference a-b in [-pi, pi).
d = atan2(sin(a-b), cos(a-b));
end

function v = local_get_cfg(cfg, path, defaultVal)
v = defaultVal;
try
    parts = strsplit(path, '.');
    cur = cfg;
    for i = 1:numel(parts)
        k = parts{i};
        if ~isstruct(cur) || ~isfield(cur, k)
            return;
        end
        cur = cur.(k);
    end
    if ~isempty(cur)
        v = cur;
    end
catch
    v = defaultVal;
end
end
