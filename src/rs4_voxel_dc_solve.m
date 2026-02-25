function R = rs4_voxel_dc_solve(SA, SB, T, cfg)
%RS4_VOXEL_DC_SOLVE  Differential-corrector solve on a voxel overlap pair.
%
% Solve for decision vector
%   u = [phi_A, delta_A, t_A, phi_B, delta_B, t_B]
% using damped least-squares / Levenberg-Marquardt on the patch residual
% between the A arc and mirrored-B arc states:
%   rx  = xA - xB
%   ry  = yA - yB
%   rth = wrapToPi(thetaA - thetaB)
%
% Optional soft voxel-neighborhood penalties can be appended to the residual
% via cfg.rs4.dc.soft_voxel_penalty.enable.

if nargin < 4
    cfg = struct();
end

% -------------------- configuration --------------------
maxIter     = local_cfg_get(cfg, 'rs4.dc.max_iter', 30);
lambda      = local_cfg_get(cfg, 'rs4.dc.lambda0', 1e-2);
lambdaMin   = local_cfg_get(cfg, 'rs4.dc.lambda_min', 1e-8);
lambdaMax   = local_cfg_get(cfg, 'rs4.dc.lambda_max', 1e8);
accThresh   = local_cfg_get(cfg, 'rs4.dc.accept_rho_min', 1e-4);
lambdaUp    = local_cfg_get(cfg, 'rs4.dc.lambda_up', 5.0);
lambdaDown  = local_cfg_get(cfg, 'rs4.dc.lambda_down', 0.3);
resTol      = local_cfg_get(cfg, 'rs4.dc.res_tol', 1e-8);
stepTol     = local_cfg_get(cfg, 'rs4.dc.step_tol', 1e-10);
fdEps       = local_cfg_get(cfg, 'rs4.dc.fd_eps', [1e-6,1e-6,1e-6,1e-6,1e-6,1e-6]);
Wdiag       = local_cfg_get(cfg, 'rs4.dc.Wdiag', [1, 1, 1]);

relTol      = local_cfg_get(cfg, 'propag.relTol', 1e-9);
absTol      = local_cfg_get(cfg, 'propag.absTol', 1e-9);
VU_mps      = local_cfg_get(cfg, 'units.VU_mps', 1.0);

usePenalty  = local_cfg_get(cfg, 'rs4.dc.soft_voxel_penalty.enable', false);
penaltyW    = local_cfg_get(cfg, 'rs4.dc.soft_voxel_penalty.weight', 1.0);

if isscalar(fdEps)
    fdEps = repmat(fdEps, 1, 6);
end
if numel(fdEps) ~= 6
    error('cfg.rs4.dc.fd_eps must be scalar or 6-vector.');
end

% -------------------- base data --------------------
seedA = T.seed_A(:);       % [x, y, th_nominal]
seedB = T.seed_B_frs(:);   % [x, y, th_nominal] on B FRS side

u = [0, T.delta_A, T.t_A, 0, T.delta_B, T.t_B]';

if numel(Wdiag) < 3
    Wdiag = [1,1,1];
end

[ra, sA, sB] = local_residual(u, SA, SB, seedA, seedB, T, relTol, absTol, usePenalty, penaltyW);
W = eye(numel(ra));
W(1,1) = Wdiag(1);
W(2,2) = Wdiag(2);
W(3,3) = Wdiag(3);
if numel(ra) > 3
    W(4:end,4:end) = penaltyW * eye(numel(ra)-3);
end

hist = repmat(struct('iter',0,'lambda',0,'res_norm',0,'cost',0,'step_norm',0,'accepted',false,'rho',NaN), maxIter, 1);

converged = false;
msg = '';

for k = 1:maxIter
    r = ra;
    cost = 0.5 * (r' * W * r);
    rNorm = norm(r);

    if rNorm < resTol
        converged = true;
        msg = 'Residual tolerance reached.';
        hist(k) = local_hist(k, lambda, rNorm, cost, 0, true, NaN);
        break;
    end

    J = local_fd_jacobian(u, r, fdEps, @(uu) local_residual(uu, SA, SB, seedA, seedB, T, relTol, absTol, usePenalty, penaltyW));

    A = J' * W * J + lambda * eye(6);
    g = J' * W * r;

    du = -A \ g;
    stepNorm = norm(du);

    if stepNorm < stepTol
        converged = true;
        msg = 'Step tolerance reached.';
        hist(k) = local_hist(k, lambda, rNorm, cost, stepNorm, true, NaN);
        break;
    end

    uTry = u + du;
    [rTry, sATry, sBTry] = local_residual(uTry, SA, SB, seedA, seedB, T, relTol, absTol, usePenalty, penaltyW);
    costTry = 0.5 * (rTry' * W * rTry);

    pred = 0.5 * du' * (lambda * du - g);
    if pred <= 0
        rho = -Inf;
    else
        rho = (cost - costTry) / pred;
    end

    accepted = isfinite(rho) && (rho > accThresh) && (costTry < cost);
    if accepted
        u = uTry;
        ra = rTry;
        sA = sATry;
        sB = sBTry;
        lambda = max(lambdaMin, lambda * lambdaDown);
    else
        lambda = min(lambdaMax, lambda * lambdaUp);
    end

    hist(k) = local_hist(k, lambda, norm(ra), 0.5*(ra' * W * ra), stepNorm, accepted, rho);
end

if ~converged && isempty(msg)
    if norm(ra) < resTol
        converged = true;
        msg = 'Residual tolerance reached at final iteration.';
    else
        msg = 'Maximum iterations reached.';
    end
end

% trim history
used = find([hist.iter] > 0);
if isempty(used)
    hist = struct([]);
else
    hist = hist(used);
end

% corrected DV values
phiA = u(1); deltaA = u(2); tA = u(3);
phiB = u(4); deltaB = u(5); tB = u(6);

[dvTurnA_nd, dvTurnB_nd] = local_dv_turns(SA, SB, seedA, seedB, phiA + deltaA, phiB + deltaB);

thetaDiff = abs(rs3_wrapToPi(sA(3) - sB(3)));
xPatch = 0.5 * (sA(1) + sB(1));
yPatch = 0.5 * (sA(2) + sB(2));
potPatch = rs3_core_cr3bp_U_and_derivs(xPatch, yPatch, SA.mu);
vPatch = sqrt(max(2 * potPatch.U - min(SA.CJ, SB.CJ), 0));
dvPatch_nd = 2 * vPatch * sin(thetaDiff / 2);

dvPatch_mps = dvPatch_nd * VU_mps;
dvTotal_nd = dvTurnA_nd + dvPatch_nd + dvTurnB_nd;
dvTotal_mps = dvTotal_nd * VU_mps;

R = struct();
R.converged = converged;
R.message = msg;
R.iterations = numel(hist);
R.final_residual_norm = norm(ra);

R.u0 = [0, T.delta_A, T.t_A, 0, T.delta_B, T.t_B];
R.u = u(:)';
R.phi_A = phiA; R.delta_A = deltaA; R.t_A = tA;
R.phi_B = phiB; R.delta_B = deltaB; R.t_B = tB;

R.stateA_patch = struct('x', sA(1), 'y', sA(2), 'theta', sA(3));
R.stateB_patch = struct('x', sB(1), 'y', sB(2), 'theta', sB(3));

R.DV_turn_A_dc_nd = dvTurnA_nd;
R.DV_turn_B_dc_nd = dvTurnB_nd;
R.DV_patch_dc_nd = dvPatch_nd;
R.DV_total_dc_nd = dvTotal_nd;

R.DV_turn_A_dc_mps = dvTurnA_nd * VU_mps;
R.DV_turn_B_dc_mps = dvTurnB_nd * VU_mps;
R.DV_patch_dc_mps = dvPatch_mps;
R.DV_total_dc_mps = dvTotal_mps;

R.history = hist;
R.residual = ra(:)';
end

% -------------------------------------------------------------------------
function [r, stateA, stateB] = local_residual(u, SA, SB, seedA, seedB, T, relTol, absTol, usePenalty, penaltyW)
phiA = u(1); deltaA = u(2); tA = u(3);
phiB = u(4); deltaB = u(5); tB = u(6);

stateA = local_arc_end_state(seedA, SA, phiA + deltaA, tA, relTol, absTol, false);
stateB = local_arc_end_state(seedB, SB, phiB + deltaB, tB, relTol, absTol, true);

rx  = stateA(1) - stateB(1);
ry  = stateA(2) - stateB(2);
rth = rs3_wrapToPi(stateA(3) - stateB(3));

r = [rx; ry; rth];

if usePenalty
    g = SA.grid3;
    hx = 0.5 * g.dx;
    hy = 0.5 * g.dy;
    hth = 0.5 * g.dtheta;

    pA = local_voxel_penalty_terms(stateA, T.xc, T.yc, T.thc, hx, hy, hth);
    pB = local_voxel_penalty_terms(stateB, T.xc, T.yc, T.thc, hx, hy, hth);

    r = [r; sqrt(max(penaltyW,0)) * [pA; pB]];
end
end

% -------------------------------------------------------------------------
function s = local_arc_end_state(seed, S, dth, tEnd, relTol, absTol, mirrorToBRS)
tEnd = max(tEnd, 0);
IC = [seed(1); seed(2); rs3_wrapToPi(seed(3) + dth)];

odeOpts = odeset('RelTol', relTol, 'AbsTol', absTol);
[~, X] = ode113(@(tt,XX) rs3_core_reduced_cr3bp_model(tt, XX, S.CJ, S.mu, false), [0, tEnd], IC, odeOpts);

x = X(end,1);
y = X(end,2);
th = X(end,3);

if mirrorToBRS
    y = -y;
    th = rs3_wrapToPi(pi - th);
else
    th = rs3_wrapToPi(th);
end

s = [x; y; th];
end

% -------------------------------------------------------------------------
function J = local_fd_jacobian(u, r0, epsVec, f)
m = numel(r0);
n = numel(u);
J = zeros(m, n);
for i = 1:n
    du = zeros(n,1);
    e = epsVec(i);
    du(i) = e;
    rp = f(u + du);
    rm = f(u - du);
    J(:,i) = (rp - rm) / (2*e);
end
end

% -------------------------------------------------------------------------
function h = local_hist(iter, lambda, rNorm, cost, stepNorm, accepted, rho)
h = struct('iter',iter, 'lambda',lambda, 'res_norm',rNorm, 'cost',cost, ...
    'step_norm',stepNorm, 'accepted',accepted, 'rho',rho);
end

% -------------------------------------------------------------------------
function p = local_voxel_penalty_terms(s, xc, yc, thc, hx, hy, hth)
dx = abs(s(1) - xc) - hx;
dy = abs(s(2) - yc) - hy;
dth = abs(rs3_wrapToPi(s(3) - thc)) - hth;

p = [max(dx,0); max(dy,0); max(dth,0)];
end

% -------------------------------------------------------------------------
function [dvA, dvB] = local_dv_turns(SA, SB, seedA, seedB, dthA, dthB)
potA = rs3_core_cr3bp_U_and_derivs(seedA(1), seedA(2), SA.mu);
potB = rs3_core_cr3bp_U_and_derivs(seedB(1), seedB(2), SB.mu);
vA = sqrt(max(2 * potA.U - SA.CJ, 0));
vB = sqrt(max(2 * potB.U - SB.CJ, 0));

dvA = 2 * vA * sin(abs(dthA) / 2);
dvB = 2 * vB * sin(abs(dthB) / 2);
end

% -------------------------------------------------------------------------
function v = local_cfg_get(cfg, path, defaultVal)
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
