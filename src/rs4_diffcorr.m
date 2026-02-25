function Tc = rs4_diffcorr(T, SA, SB, cfg)
%RS4_DIFFCORR  Differential correction of the voxel-overlap transfer pair.
%
% Starting from the argmin-DV arcs in T (from rs4_voxel_traj_extract),
% refines the 6-vector
%
%   z = [alpha_A, delta_A, t_A, alpha_B, delta_B, t_B]
%
% via fmincon to minimise DV_turn_A + DV_turn_B subject to the position +
% heading continuity constraint at the arc endpoints:
%
%   r(z) = [ x_A(t_A)  -  x_B(t_B)           ]
%           [ y_A(t_A)  -  y_B(t_B)           ]  =  0
%           [ wrap( th_A(t_A) - th_B(t_B) )   ]
%
% where
%   alpha_A/B   continuous departure-time parameters along PO_A / PO_B
%               (interpolated from SA.Xpo / SB.Xpo via SA.t_dense)
%   delta_A/B   heading kicks applied at departure
%   t_A / t_B   forward arc-integration times
%
% and the B arc is produced by the R-transform of the FRS arc from PO_B:
%   (x_B, y_B, th_B)  =  (x_Bfrs, -y_Bfrs, pi - th_Bfrs)
%
% The residual is normalised by the voxel grid spacing before being passed
% to fmincon so that all three components are O(1) at the warm start.
%
% Inputs
%   T   : struct from rs4_voxel_traj_extract          (warm start)
%   SA  : family A struct  (.Xpo, .t_dense, .Tf_PO, .CJ, .mu, .grid3, ...)
%   SB  : family B struct
%   cfg : config struct
%
% Output
%   Tc  : corrected trajectory struct  (same schema as T, plus .dc fields)
%
% Optional cfg fields (with defaults):
%   cfg.diffcorr.tol_patch   constraint tolerance in normalised units (1e-6)

% -------------------------------------------------------------------------
% Config
% -------------------------------------------------------------------------
relTol    = local_cfg_get(cfg, 'propag.relTol',      1e-9);
absTol    = local_cfg_get(cfg, 'propag.absTol',      1e-9);
VU_mps    = local_cfg_get(cfg, 'units.VU_mps',       1.0);
TU_days   = local_cfg_get(cfg, 'units.TU_days',      1.0);
Tmax      = local_cfg_get(cfg, 'propag.Tmax',        pi/2);
tol_patch = local_cfg_get(cfg, 'diffcorr.tol_patch', 1e-6);

odeOpts = odeset('RelTol', relTol, 'AbsTol', absTol);

% -------------------------------------------------------------------------
% Step 1 — Build initial guess from T
% -------------------------------------------------------------------------
alpha_A_0 = local_find_alpha(SA, T.seed_A(1:2));
alpha_B_0 = local_find_alpha(SB, T.seed_B_frs(1:2));
delta_A_0 = rs3_wrapToPi(T.IC_A(3)     - T.seed_A(3));
delta_B_0 = rs3_wrapToPi(T.IC_B_frs(3) - T.seed_B_frs(3));
t_A_0     = T.t_A;
t_B_0     = T.t_B;

z0 = [alpha_A_0; delta_A_0; t_A_0; alpha_B_0; delta_B_0; t_B_0];

% -------------------------------------------------------------------------
% Step 2 — Bounds
%   alpha  in [0, Tf_PO]   (wraps on PO)
%   delta  in [-pi/2, pi/2]  (generous; objective naturally limits kicks)
%   t      in [1e-4, Tmax]
% -------------------------------------------------------------------------
delta_max = pi / 2;
lb = [0;          -delta_max;  1e-4;   0;          -delta_max;  1e-4 ];
ub = [SA.Tf_PO;   +delta_max;  Tmax;   SB.Tf_PO;   +delta_max;  Tmax ];

% -------------------------------------------------------------------------
% Step 3 — Residual scaling (normalise to voxel-grid units)
%   r_scaled = r ./ [dx; dy; dtheta]  =>  O(1) at initial miss
% -------------------------------------------------------------------------
grid3   = SA.grid3;
r_scale = [grid3.dx; grid3.dy; grid3.dtheta];

% -------------------------------------------------------------------------
% Step 4 — Diagnostic: evaluate warm-start residual
% -------------------------------------------------------------------------
r0 = local_residual(z0, SA, SB, odeOpts);
fprintf('\n[diffcorr] ---- Differential Correction ----\n');
fprintf('[diffcorr] Warm start:\n');
fprintf('[diffcorr]   DV_total = %.3f m/s  (turn_A=%.3f  patch=%.3f  turn_B=%.3f)\n', ...
    T.DV_total_true_mps, T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps);
fprintf('[diffcorr]   residual (raw)     = [%.3e  %.3e  %.3e]\n', r0(1), r0(2), r0(3));
fprintf('[diffcorr]   residual (scaled)  = %.3e\n', norm(r0 ./ r_scale));

% -------------------------------------------------------------------------
% Step 5 — fmincon
% -------------------------------------------------------------------------
obj_fun = @(z) local_objective(z, SA, SB);
con_fun = @(z) local_constraints(z, SA, SB, odeOpts, r_scale);

opts_fmin = optimoptions('fmincon', ...
    'Algorithm',               'sqp',   ...
    'Display',                 'iter',  ...
    'ConstraintTolerance',     tol_patch, ...
    'OptimalityTolerance',     1e-7,    ...
    'StepTolerance',           1e-10,   ...
    'MaxIterations',           300,     ...
    'MaxFunctionEvaluations',  5000,    ...
    'FiniteDifferenceType',    'central', ...
    'FiniteDifferenceStepSize', 1e-6);

[z_star, ~, exitflag, output] = fmincon( ...
    obj_fun, z0, [], [], [], [], lb, ub, con_fun, opts_fmin);

fprintf('[diffcorr] fmincon done: exitflag=%d  iter=%d  fevals=%d\n', ...
    exitflag, output.iterations, output.funcCount);

if exitflag <= 0
    warning('[diffcorr] fmincon did not converge (exitflag=%d). Tc may be inaccurate.', exitflag);
end

% -------------------------------------------------------------------------
% Step 6 — Unpack corrected solution and re-integrate dense arcs
% -------------------------------------------------------------------------
alpha_A = z_star(1);  delta_A = z_star(2);  t_A = z_star(3);
alpha_B = z_star(4);  delta_B = z_star(5);  t_B = z_star(6);

seed_A   = local_interp_po(SA, alpha_A);
IC_A     = [seed_A(1); seed_A(2); rs3_wrapToPi(seed_A(3) + delta_A)];

seed_B   = local_interp_po(SB, alpha_B);
IC_B_frs = [seed_B(1); seed_B(2); rs3_wrapToPi(seed_B(3) + delta_B)];

[tA_vec, XA]     = ode113(@(tt,XX) rs3_core_reduced_cr3bp_model(tt,XX,SA.CJ,SA.mu,false), ...
                           [0, t_A], IC_A, odeOpts);
[tB_vec, XB_frs] = ode113(@(tt,XX) rs3_core_reduced_cr3bp_model(tt,XX,SB.CJ,SB.mu,false), ...
                           [0, t_B], IC_B_frs, odeOpts);

% R-transform: (x, y, th)_frs -> (x, -y, pi-th)_brs
x_B  =  XB_frs(:,1);
y_B  = -XB_frs(:,2);
th_B =  rs3_wrapToPi(pi - XB_frs(:,3));

% Patch point — endpoint of A (= endpoint of B within tol)
xp    = XA(end,1);
yp    = XA(end,2);
thp_A = XA(end,3);
thp_B = th_B(end);

% -------------------------------------------------------------------------
% Step 7 — DV computations at corrected solution
% -------------------------------------------------------------------------
pot_A        = rs3_core_cr3bp_U_and_derivs(seed_A(1), seed_A(2), SA.mu);
v0_A         = sqrt(max(2*pot_A.U - SA.CJ, 0));
DV_turn_A_nd = 2 * v0_A * sin(abs(delta_A) / 2);

pot_B        = rs3_core_cr3bp_U_and_derivs(seed_B(1), seed_B(2), SB.mu);
v0_B         = sqrt(max(2*pot_B.U - SB.CJ, 0));
DV_turn_B_nd = 2 * v0_B * sin(abs(delta_B) / 2);

pot_p        = rs3_core_cr3bp_U_and_derivs(xp, yp, SA.mu);
v_patch      = sqrt(max(2*pot_p.U - min(SA.CJ, SB.CJ), 0));
delta_th     = abs(rs3_wrapToPi(thp_A - thp_B));
DV_patch_nd  = 2 * v_patch * sin(delta_th / 2);

r_final = local_residual(z_star, SA, SB, odeOpts);

fprintf('[diffcorr] Corrected solution:\n');
fprintf('[diffcorr]   DV_total = %.3f m/s  (turn_A=%.3f  patch=%.3f  turn_B=%.3f)\n', ...
    (DV_turn_A_nd + DV_patch_nd + DV_turn_B_nd)*VU_mps, ...
    DV_turn_A_nd*VU_mps, DV_patch_nd*VU_mps, DV_turn_B_nd*VU_mps);
fprintf('[diffcorr]   residual (raw)     = [%.3e  %.3e  %.3e]\n', r_final(1), r_final(2), r_final(3));
fprintf('[diffcorr]   residual (scaled)  = %.3e  (tol=%.2e)\n', ...
    norm(r_final ./ r_scale), tol_patch);
fprintf('[diffcorr]   delta_th at patch  = %.4f deg\n', rad2deg(delta_th));

% -------------------------------------------------------------------------
% Pack output struct
% -------------------------------------------------------------------------
Tc = struct();

% Corrected arc A
Tc.alpha_A      = alpha_A;
Tc.delta_A      = delta_A;
Tc.t_A          = t_A;
Tc.IC_A         = IC_A;
Tc.seed_A       = seed_A;
Tc.tA_vec       = tA_vec;
Tc.XA           = XA;

% Corrected arc B
Tc.alpha_B      = alpha_B;
Tc.delta_B      = delta_B;
Tc.t_B          = t_B;
Tc.IC_B_frs     = IC_B_frs;
Tc.seed_B_frs   = seed_B;
Tc.tB_vec       = tB_vec;
Tc.x_B          = x_B;
Tc.y_B          = y_B;
Tc.th_B         = th_B;

% Patch point
Tc.xp           = xp;
Tc.yp           = yp;
Tc.thp_A        = thp_A;
Tc.thp_B        = thp_B;
Tc.delta_th_rad = delta_th;
Tc.r_final      = r_final;

% DV
Tc.DV_turn_A_nd  = DV_turn_A_nd;
Tc.DV_turn_B_nd  = DV_turn_B_nd;
Tc.DV_patch_nd   = DV_patch_nd;
Tc.DV_turn_A_mps = DV_turn_A_nd * VU_mps;
Tc.DV_turn_B_mps = DV_turn_B_nd * VU_mps;
Tc.DV_patch_mps  = DV_patch_nd  * VU_mps;
Tc.DV_total_mps  = (DV_turn_A_nd + DV_patch_nd + DV_turn_B_nd) * VU_mps;

% TOF
Tc.tof_A_days   = t_A * TU_days;
Tc.tof_B_days   = t_B * TU_days;

% Optimisation metadata
Tc.exitflag     = exitflag;
Tc.iterations   = output.iterations;
Tc.tol_patch    = tol_patch;
end

% =========================================================================
% Local helpers
% =========================================================================

function dv = local_objective(z, SA, SB)
% Minimise total turn cost (patch cost goes to zero via constraint).
    alpha_A = z(1);  delta_A = z(2);
    alpha_B = z(4);  delta_B = z(5);
    seed_A  = local_interp_po(SA, alpha_A);
    seed_B  = local_interp_po(SB, alpha_B);
    pot_A   = rs3_core_cr3bp_U_and_derivs(seed_A(1), seed_A(2), SA.mu);
    v0_A    = sqrt(max(2*pot_A.U - SA.CJ, 0));
    pot_B   = rs3_core_cr3bp_U_and_derivs(seed_B(1), seed_B(2), SB.mu);
    v0_B    = sqrt(max(2*pot_B.U - SB.CJ, 0));
    dv      = 2*v0_A*sin(abs(delta_A)/2) + 2*v0_B*sin(abs(delta_B)/2);
end

% -------------------------------------------------------------------------

function [c, ceq] = local_constraints(z, SA, SB, odeOpts, r_scale)
    c   = [];
    r   = local_residual(z, SA, SB, odeOpts);
    ceq = r ./ r_scale;      % normalised to O(1) at warm start
end

% -------------------------------------------------------------------------

function r = local_residual(z, SA, SB, odeOpts)
    alpha_A = z(1);  delta_A = z(2);  t_A = z(3);
    alpha_B = z(4);  delta_B = z(5);  t_B = z(6);

    seed_A   = local_interp_po(SA, alpha_A);
    IC_A     = [seed_A(1); seed_A(2); rs3_wrapToPi(seed_A(3) + delta_A)];

    seed_B   = local_interp_po(SB, alpha_B);
    IC_B_frs = [seed_B(1); seed_B(2); rs3_wrapToPi(seed_B(3) + delta_B)];

    [~, XA]     = ode113(@(tt,XX) rs3_core_reduced_cr3bp_model(tt,XX,SA.CJ,SA.mu,false), ...
                          [0, max(t_A, 1e-6)], IC_A, odeOpts);
    [~, XB_frs] = ode113(@(tt,XX) rs3_core_reduced_cr3bp_model(tt,XX,SB.CJ,SB.mu,false), ...
                          [0, max(t_B, 1e-6)], IC_B_frs, odeOpts);

    xA  = XA(end,1);
    yA  = XA(end,2);
    thA = XA(end,3);

    xB  =  XB_frs(end,1);
    yB  = -XB_frs(end,2);                       % R-transform y
    thB =  rs3_wrapToPi(pi - XB_frs(end,3));    % R-transform th

    r = [xA - xB; yA - yB; rs3_wrapToPi(thA - thB)];
end

% -------------------------------------------------------------------------

function seed = local_interp_po(S, alpha)
% Interpolate PO state [x, y, th] at continuous time alpha in [0, Tf_PO).
% Uses pchip for (x,y) and linear on unwrapped theta to avoid jump artefacts.
    t   = mod(alpha, S.Tf_PO);
    xy  = interp1(S.t_dense, S.Xpo(:,1:2), t, 'pchip');
    th_uw = unwrap(S.Xpo(:,3));
    th  = rs3_wrapToPi(interp1(S.t_dense, th_uw, t, 'linear'));
    seed = [xy(1), xy(2), th];
end

% -------------------------------------------------------------------------

function alpha = local_find_alpha(S, xy)
% Return the t_dense value at which PO is closest to 2-D position xy.
    d = hypot(S.Xpo(:,1) - xy(1), S.Xpo(:,2) - xy(2));
    [~, idx] = min(d);
    alpha = S.t_dense(idx);
end

% -------------------------------------------------------------------------

function v = local_cfg_get(cfg, path, defaultVal)
    v = defaultVal;
    try
        parts = strsplit(path, '.');
        cur = cfg;
        for i = 1:numel(parts)
            k = parts{i};
            if ~isstruct(cur) || ~isfield(cur, k), return; end
            cur = cur.(k);
        end
        if ~isempty(cur), v = cur; end
    catch
        v = defaultVal;
    end
end
