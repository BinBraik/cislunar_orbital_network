function Tc = traj_diffcorr(T, SA, SB, cfg)
%RS4_DIFFCORR  Differential correction of the voxel-overlap transfer pair.
%
% Starting from the argmin-DV arcs in T (from overlap_voxel_traj_extract),
% refines the 6-vector
%
%   z = [alpha_A, delta_A, t_A, alpha_B, delta_B, t_B]
%
% via fmincon.  Two modes are supported (cfg.diffcorr.relax_heading):
%
%   relax_heading = false (default)
%     Minimise DV_turn_A + DV_turn_B subject to full position + heading
%     continuity (forces DV_patch = 0 via constraint):
%       r(z) = [ x_A(t_A) - x_B(t_B)          ]
%              [ y_A(t_A) - y_B(t_B)          ]  = 0
%              [ wrap( th_A(t_A) - th_B(t_B)) ]
%
%   relax_heading = true
%     Minimise DV_turn_A + DV_patch + DV_turn_B subject to spatial
%     continuity only (heading mismatch enters the cost, not the constraint):
%       r(z) = [ x_A(t_A) - x_B(t_B) ]  = 0
%              [ y_A(t_A) - y_B(t_B) ]
%     This is a relaxation of the default mode — the feasible set is
%     strictly larger so total DV can only be equal or lower.
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
%   T   : struct from overlap_voxel_traj_extract          (warm start)
%   SA  : family A struct  (.Xpo, .t_dense, .Tf_PO, .CJ, .mu, .grid3, ...)
%   SB  : family B struct
%   cfg : config struct
%
% Output
%   Tc  : corrected trajectory struct  (same schema as T, plus .dc fields)
%
% Optional cfg fields (with defaults):
%   cfg.diffcorr.tol_patch      constraint tolerance in normalised units (1e-4)
%   cfg.diffcorr.tol_converged  "CONVERGED" label threshold — can be >= tol_patch (1e-4)
%   cfg.diffcorr.display        fmincon display mode: 'off'|'iter'|'final' ('iter')
%   cfg.diffcorr.MaxIterations  fmincon max iterations (300)
%   cfg.diffcorr.MaxFunEvals    fmincon max function evaluations (8000)
%   cfg.diffcorr.N_po_dt        target time between PO knots [ND] — scales with period (0.003)
%   cfg.diffcorr.N_po_min       minimum knot count floor (1001)
%   cfg.diffcorr.relax_heading  if true: 2-component spatial constraint +
%                               DV_patch in objective (false)

% -------------------------------------------------------------------------
% Config
% -------------------------------------------------------------------------
relTol_final = local_cfg_get(cfg, 'propag.relTol',      1e-9);
absTol_final = local_cfg_get(cfg, 'propag.absTol',      1e-9);
VU_mps       = local_cfg_get(cfg, 'units.VU_mps',       1.0);
TU_days      = local_cfg_get(cfg, 'units.TU_days',      1.0);
Tmax         = local_cfg_get(cfg, 'propag.Tmax',        pi/2);
tol_patch     = local_cfg_get(cfg, 'diffcorr.tol_patch',     1e-4);
tol_converged = local_cfg_get(cfg, 'diffcorr.tol_converged', tol_patch);
dc_display    = local_cfg_get(cfg, 'diffcorr.display',        'iter');  % 'off' for batch

% Two-tolerance strategy: use 10x looser tolerances during the many ODE
% evaluations inside the optimizer (FD gradients), then tight tolerances
% only for the final arc re-integration and residual check.
% The FD step is 1e-6, so 1e-7 arc accuracy is more than sufficient for
% gradient fidelity while reducing per-arc cost noticeably.
relTol_dc    = local_cfg_get(cfg, 'diffcorr.relTol_dc',      relTol_final * 10);
absTol_dc    = local_cfg_get(cfg, 'diffcorr.absTol_dc',      absTol_final * 10);
dc_max_iter    = local_cfg_get(cfg, 'diffcorr.MaxIterations',   300);
dc_max_feval   = local_cfg_get(cfg, 'diffcorr.MaxFunEvals',     8000);
relax_heading  = local_cfg_get(cfg, 'diffcorr.relax_heading',   false);
relTol     = relTol_final;  % alias kept for local_ensure_xpo calls
absTol     = absTol_final;

odeOpts       = odeset('RelTol', relTol_dc,    'AbsTol', absTol_dc);    % optimizer evals
odeOpts_final = odeset('RelTol', relTol_final, 'AbsTol', absTol_final); % final arcs

% -------------------------------------------------------------------------
% Step 0 — Ensure SA/SB have dense PO trace (Xpo, t_dense)
%   The cache strips these by default (cfg.cache.store_dense_po=false).
%   Re-integrate on demand using X0, CJ, mu, Tf_PO — all of which are cached.
% -------------------------------------------------------------------------
% Scale knot count to orbit period so all families get ~N_po_dt ND time between
% knots — avoids Cycler 21 (T≈21) being 7x coarser than Lyapunov L1 (T≈3).
N_po_dt  = local_cfg_get(cfg, 'diffcorr.N_po_dt',  0.003);   % target spacing [ND]
N_po_min = local_cfg_get(cfg, 'diffcorr.N_po_min', 1001);    % floor
N_po_A   = max(N_po_min, ceil(SA.Tf_PO / N_po_dt) + 1);
N_po_B   = max(N_po_min, ceil(SB.Tf_PO / N_po_dt) + 1);
SA       = local_ensure_xpo(SA, relTol, absTol, N_po_A);
SB       = local_ensure_xpo(SB, relTol, absTol, N_po_B);

% -------------------------------------------------------------------------
% Step 1 — Build initial guess from T
% -------------------------------------------------------------------------
alpha_A_0 = local_find_alpha_cont(SA, T.seed_A(1:2));
alpha_B_0 = local_find_alpha_cont(SB, T.seed_B_frs(1:2));
delta_A_0 = wrap_to_pi(T.IC_A(3)     - T.seed_A(3));
delta_B_0 = wrap_to_pi(T.IC_B_frs(3) - T.seed_B_frs(3));
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
%   When relax_heading=true the heading row is dropped from the constraint,
%   so r_scale is 2-element.
% -------------------------------------------------------------------------
grid3   = SA.grid3;
if relax_heading
    r_scale = [grid3.dx; grid3.dy];
else
    r_scale = [grid3.dx; grid3.dy; grid3.dtheta];
end

% -------------------------------------------------------------------------
% Step 4 — Diagnostic: evaluate warm-start residual
% -------------------------------------------------------------------------
r0 = local_residual(z0, SA, SB, odeOpts, relax_heading);
if ~strcmp(dc_display, 'off')
    fprintf('\n[diffcorr] ---- Differential Correction ----\n');
    if relax_heading
        fprintf('[diffcorr] Mode: relax_heading (spatial constraint only; DV_patch in cost)\n');
    end
    fprintf('[diffcorr] Warm start:\n');
    fprintf('[diffcorr]   DV_total = %.3f m/s  (turn_A=%.3f  patch=%.3f  turn_B=%.3f)\n', ...
        T.DV_total_true_mps, T.DV_turn_A_mps, T.DV_patch_mps, T.DV_turn_B_mps);
    fprintf('[diffcorr]   residual (raw)     = %s\n', mat2str(r0', 3));
    fprintf('[diffcorr]   residual (scaled)  = %.3e\n', norm(r0 ./ r_scale));
end

% -------------------------------------------------------------------------
% Step 5 — fmincon (pass 1)
% -------------------------------------------------------------------------
obj_fun = @(z) local_objective(z, SA, SB, odeOpts, relax_heading);
con_fun = @(z) local_constraints(z, SA, SB, odeOpts, r_scale, relax_heading);

% TypicalX guides the per-variable FD step size.
% alpha [0,Tf_PO] ~ Tf_PO/2,  delta [-pi/2,pi/2] ~ pi/6,  t [1e-4,Tmax] ~ Tmax/2
typX = [SA.Tf_PO/2; pi/6; Tmax/2; SB.Tf_PO/2; pi/6; Tmax/2];

opts_fmin = optimoptions('fmincon', ...
    'Algorithm',               'sqp',        ...
    'Display',                 dc_display,   ...
    'ConstraintTolerance',     tol_patch, ...
    'OptimalityTolerance',     1e-7,    ...
    'StepTolerance',           1e-10,   ...
    'MaxIterations',           dc_max_iter,   ...
    'MaxFunctionEvaluations',  dc_max_feval,  ...
    'FiniteDifferenceType',    'central', ...
    'FiniteDifferenceStepSize', 1e-6,   ...
    'TypicalX',                typX);

[z_star, ~, exitflag, output] = fmincon( ...
    obj_fun, z0, [], [], [], [], lb, ub, con_fun, opts_fmin);

if ~strcmp(dc_display, 'off')
    fprintf('[diffcorr] fmincon pass 1: exitflag=%d  iter=%d  fevals=%d\n', ...
        exitflag, output.iterations, output.funcCount);
end

% -------------------------------------------------------------------------
% Step 5b — Pre-feasibility fallback
%
%   exitflag=-2 : SQP constraint-recovery phase failed from the warm start.
%                 The miss-distance is too large for SQP to bridge with its
%                 quadratic sub-problems.
%   exitflag= 0 : budget exhausted before residual converged.
%
%   Fix: run lsqnonlin on ||r||^2 alone (no DV objective) to drive the
%        patch residual to near-zero, then restart fmincon from that
%        feasible warm point to recover DV optimality.
% -------------------------------------------------------------------------
r_pass1      = local_residual(z_star, SA, SB, odeOpts, relax_heading);
r_pass1_norm = norm(r_pass1 ./ r_scale);

z_best      = z_star;
r_best_norm = r_pass1_norm;

if r_pass1_norm > tol_patch
    if ~strcmp(dc_display, 'off')
        fprintf('[diffcorr] Pass 1 infeasible (||r_sc||=%.2e, exitflag=%d). Running lsqnonlin ...\n', ...
            r_pass1_norm, exitflag);
    end

    % Warm-start: use z_star whenever fmincon made any progress (any exitflag
    % except -2 which means it never moved), otherwise fall back to z0.
    z_lsq0 = z0;
    if exitflag ~= -2, z_lsq0 = z_star; end

    opts_lsq = optimoptions('lsqnonlin', ...
        'Display',                'off', ...
        'FunctionTolerance',      tol_patch^2, ...
        'StepTolerance',          1e-10, ...
        'MaxFunctionEvaluations', 5000, ...
        'FiniteDifferenceType',   'central');

    z_feas      = lsqnonlin(@(z) local_residual(z, SA, SB, odeOpts, relax_heading) ./ r_scale, ...
                             z_lsq0, lb, ub, opts_lsq);
    r_feas_norm = norm(local_residual(z_feas, SA, SB, odeOpts, relax_heading) ./ r_scale);
    if ~strcmp(dc_display, 'off')
        fprintf('[diffcorr] lsqnonlin: ||r_sc||=%.2e\n', r_feas_norm);
    end

    if r_feas_norm < r_best_norm
        z_best      = z_feas;
        r_best_norm = r_feas_norm;

        % Restart fmincon from the feasible point to optimise DV.
        if ~strcmp(dc_display, 'off')
            fprintf('[diffcorr] Restarting fmincon from feasible point ...\n');
        end
        opts_fmin2 = optimoptions(opts_fmin, 'Display', 'off');
        [z_star2, ~, exitflag2, output2] = fmincon( ...
            obj_fun, z_feas, [], [], [], [], lb, ub, con_fun, opts_fmin2);

        r_pass2_norm = norm(local_residual(z_star2, SA, SB, odeOpts, relax_heading) ./ r_scale);
        if ~strcmp(dc_display, 'off')
            fprintf('[diffcorr] fmincon pass 2: exitflag=%d  iter=%d  ||r_sc||=%.2e\n', ...
                exitflag2, output2.iterations, r_pass2_norm);
        end

        if r_pass2_norm <= r_best_norm + tol_patch
            z_best      = z_star2;
            exitflag    = exitflag2;
            r_best_norm = r_pass2_norm;
        end
    end
end

z_star = z_best;

% -------------------------------------------------------------------------
% Step 6 — Unpack corrected solution and re-integrate dense arcs
% -------------------------------------------------------------------------
alpha_A = z_star(1);  delta_A = z_star(2);  t_A = z_star(3);
alpha_B = z_star(4);  delta_B = z_star(5);  t_B = z_star(6);

seed_A   = local_interp_po(SA, alpha_A);
IC_A     = [seed_A(1); seed_A(2); wrap_to_pi(seed_A(3) + delta_A)];

seed_B   = local_interp_po(SB, alpha_B);
IC_B_frs = [seed_B(1); seed_B(2); wrap_to_pi(seed_B(3) + delta_B)];

[tA_vec, XA]     = ode113(@(tt,XX) cr3bp_reduced_ode(tt,XX,SA.CJ,SA.mu,false), ...
                           [0, t_A], IC_A, odeOpts_final);
[tB_vec, XB_frs] = ode113(@(tt,XX) cr3bp_reduced_ode(tt,XX,SB.CJ,SB.mu,false), ...
                           [0, t_B], IC_B_frs, odeOpts_final);

% R-transform: (x, y, th)_frs -> (x, -y, pi-th)_brs
x_B  =  XB_frs(:,1);
y_B  = -XB_frs(:,2);
th_B =  wrap_to_pi(pi - XB_frs(:,3));

% Patch point — endpoint of A (= endpoint of B within tol)
xp    = XA(end,1);
yp    = XA(end,2);
thp_A = XA(end,3);
thp_B = th_B(end);

% -------------------------------------------------------------------------
% Step 7 — DV computations at corrected solution
% -------------------------------------------------------------------------
pot_A        = cr3bp_potential(seed_A(1), seed_A(2), SA.mu);
v0_A         = sqrt(max(2*pot_A.U - SA.CJ, 0));
DV_turn_A_nd = 2 * v0_A * sin(abs(delta_A) / 2);

pot_B        = cr3bp_potential(seed_B(1), seed_B(2), SB.mu);
v0_B         = sqrt(max(2*pot_B.U - SB.CJ, 0));
DV_turn_B_nd = 2 * v0_B * sin(abs(delta_B) / 2);

pot_p        = cr3bp_potential(xp, yp, SA.mu);
v_patch      = sqrt(max(2*pot_p.U - min(SA.CJ, SB.CJ), 0));
delta_th     = abs(wrap_to_pi(thp_A - thp_B));
DV_patch_nd  = 2 * v_patch * sin(delta_th / 2);

r_final      = local_residual(z_star, SA, SB, odeOpts_final, relax_heading);
r_final_norm = norm(r_final ./ r_scale);
converged    = r_final_norm <= tol_converged;

if ~strcmp(dc_display, 'off')
    fprintf('[diffcorr] Corrected solution:\n');
    fprintf('[diffcorr]   DV_total = %.3f m/s  (turn_A=%.3f  patch=%.3f  turn_B=%.3f)\n', ...
        (DV_turn_A_nd + DV_patch_nd + DV_turn_B_nd)*VU_mps, ...
        DV_turn_A_nd*VU_mps, DV_patch_nd*VU_mps, DV_turn_B_nd*VU_mps);
    fprintf('[diffcorr]   residual (raw)     = %s\n', mat2str(r_final', 3));
    fprintf('[diffcorr]   residual (scaled)  = %.3e  (tol_conv=%.2e)  %s\n', ...
        r_final_norm, tol_converged, local_conv_str(converged));
    fprintf('[diffcorr]   delta_th at patch  = %.4f deg\n', rad2deg(delta_th));
end

if ~converged
    warning('[diffcorr] Patch constraint NOT satisfied (||r_sc||=%.2e > tol=%.2e, exitflag=%d).', ...
        r_final_norm, tol_converged, exitflag);
end

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
Tc.converged      = converged;   % true iff ||r_final ./ r_scale|| <= tol_converged
Tc.r_final_norm   = r_final_norm;
Tc.iterations     = output.iterations;
Tc.tol_patch      = tol_patch;
Tc.tol_converged  = tol_converged;
end

% =========================================================================
% Local helpers
% =========================================================================

function dv = local_objective(z, SA, SB, odeOpts, relax_heading)
% Minimise total DV cost.
% relax_heading=false : DV_turn_A + DV_turn_B  (patch zeroed via constraint)
% relax_heading=true  : DV_turn_A + DV_patch + DV_turn_B  (heading is free)
    e       = local_eval_arcs(z, SA, SB, odeOpts);
    delta_A = z(2);
    delta_B = z(5);
    pot_A   = cr3bp_potential(e.seed_A(1), e.seed_A(2), SA.mu);
    v0_A    = sqrt(max(2*pot_A.U - SA.CJ, 0));
    pot_B   = cr3bp_potential(e.seed_B(1), e.seed_B(2), SB.mu);
    v0_B    = sqrt(max(2*pot_B.U - SB.CJ, 0));
    DV_turns = 2*v0_A*sin(abs(delta_A)/2) + 2*v0_B*sin(abs(delta_B)/2);
    if relax_heading
        pot_p    = cr3bp_potential(e.xA, e.yA, SA.mu);
        v_p      = sqrt(max(2*pot_p.U - min(SA.CJ, SB.CJ), 0));
        dth      = abs(wrap_to_pi(e.thA - e.thB));
        DV_patch = 2 * v_p * sin(dth / 2);
        dv = DV_turns + DV_patch;
    else
        dv = DV_turns;
    end
end

% -------------------------------------------------------------------------

function [c, ceq] = local_constraints(z, SA, SB, odeOpts, r_scale, relax_heading)
    c   = [];
    r   = local_residual(z, SA, SB, odeOpts, relax_heading);
    ceq = r ./ r_scale;      % normalised to O(1) at warm start
end

% -------------------------------------------------------------------------

function r = local_residual(z, SA, SB, odeOpts, relax_heading)
% Returns spatial (and optionally heading) endpoint mismatch.
% relax_heading=false : r = [dx; dy; d_theta]  (3-component)
% relax_heading=true  : r = [dx; dy]           (2-component)
    e = local_eval_arcs(z, SA, SB, odeOpts);
    if relax_heading
        r = [e.xA - e.xB; e.yA - e.yB];
    else
        r = [e.xA - e.xB; e.yA - e.yB; wrap_to_pi(e.thA - e.thB)];
    end
end

% -------------------------------------------------------------------------

function e = local_eval_arcs(z, SA, SB, odeOpts)
% Integrate both arcs from current decision vector z and return their
% endpoints.  Uses a persistent cache keyed on z to avoid redundant ODE
% evaluations when fmincon queries the objective and constraint at the same
% point (which happens every iteration with central finite differences).
    persistent z_cache e_cache
    if ~isempty(z_cache) && isequal(z, z_cache)
        e = e_cache;
        return;
    end
    alpha_A = z(1);  delta_A = z(2);  t_A = z(3);
    alpha_B = z(4);  delta_B = z(5);  t_B = z(6);

    seed_A   = local_interp_po(SA, alpha_A);
    IC_A     = [seed_A(1); seed_A(2); wrap_to_pi(seed_A(3) + delta_A)];

    seed_B   = local_interp_po(SB, alpha_B);
    IC_B_frs = [seed_B(1); seed_B(2); wrap_to_pi(seed_B(3) + delta_B)];

    [~, XA]     = ode113(@(tt,XX) cr3bp_reduced_ode(tt,XX,SA.CJ,SA.mu,false), ...
                          [0, max(t_A, 1e-6)], IC_A, odeOpts);
    [~, XB_frs] = ode113(@(tt,XX) cr3bp_reduced_ode(tt,XX,SB.CJ,SB.mu,false), ...
                          [0, max(t_B, 1e-6)], IC_B_frs, odeOpts);

    e.xA    = XA(end,1);
    e.yA    = XA(end,2);
    e.thA   = XA(end,3);
    e.xB    =  XB_frs(end,1);
    e.yB    = -XB_frs(end,2);                       % R-transform y
    e.thB   =  wrap_to_pi(pi - XB_frs(end,3));    % R-transform th
    e.seed_A = seed_A;
    e.seed_B = seed_B;

    z_cache = z;
    e_cache = e;
end

% -------------------------------------------------------------------------

function seed = local_interp_po(S, alpha)
% Interpolate PO state [x, y, th] at continuous time alpha in [0, Tf_PO).
% Uses pre-built pchip pp-forms (pp_xy, pp_th) for fast scalar queries —
% avoids the per-call interp1 setup and unwrap over 1001 points.
    t  = mod(alpha, S.Tf_PO);
    if isfield(S, 'pp_xy') && ~isempty(S.pp_xy)
        xy = ppval(S.pp_xy, t)';          % [1x2]
        th = wrap_to_pi(ppval(S.pp_th, t));
    else
        % Fallback (pp-forms not yet built — should not happen after local_ensure_xpo).
        xy    = interp1(S.t_dense, S.Xpo(:,1:2), t, 'pchip');
        th_uw = unwrap(S.Xpo(:,3));
        th    = wrap_to_pi(interp1(S.t_dense, th_uw, t, 'linear'));
    end
    seed = [xy(1), xy(2), th];
end

% -------------------------------------------------------------------------

function S = local_ensure_xpo(S, relTol, absTol, N_po)
% If SA.Xpo / SA.t_dense were stripped by the cache, rebuild them by
% re-integrating the periodic orbit from S.X0 over [0, S.Tf_PO].
% Also builds pp-form splines (pp_xy, pp_th) for fast repeated scalar
% interpolation inside the optimizer — replaces the per-call interp1+unwrap.
    need_pp = ~isfield(S, 'pp_xy') || isempty(S.pp_xy);

    if isfield(S, 'Xpo') && ~isempty(S.Xpo) && ...
       isfield(S, 't_dense') && ~isempty(S.t_dense)
        % Xpo already present — only (re)build pp-forms if missing.
        if need_pp
            S = local_build_pp(S);
        end
        return;
    end
    fprintf('[diffcorr] Xpo not in cache for "%s" — re-integrating PO ...\n', S.name);
    opts   = odeset('RelTol', relTol, 'AbsTol', absTol);
    solPO  = ode113(@(t,X) cr3bp_reduced_ode(t, X, S.CJ, S.mu, true), ...
                    [0, S.Tf_PO], S.X0, opts);
    t_dense = linspace(0, S.Tf_PO, N_po)';
    Xpo     = deval(solPO, t_dense)';
    S.t_dense = t_dense;
    S.Xpo     = Xpo;
    S = local_build_pp(S);
end

% -------------------------------------------------------------------------

function S = local_build_pp(S)
% Build pchip pp-form structs for fast repeated scalar interpolation.
% pp_xy : 2-D pchip over (x,y);  pp_th : 1-D pchip over unwrapped theta.
% ppval(S.pp_xy, t) is ~5-10x faster than interp1(...,'pchip') for scalar t,
% and avoids the per-call unwrap(1001 elements) that was in local_interp_po.
    t   = S.t_dense(:)';                      % [1 x N] row
    xy  = S.Xpo(:, 1:2)';                     % [2 x N]
    th  = unwrap(S.Xpo(:, 3))';               % [1 x N] (unwrap once, here)
    S.pp_xy = pchip(t, xy);
    S.pp_th = pchip(t, th);
end

% -------------------------------------------------------------------------

function alpha = local_find_alpha_cont(S, xy)
% Continuous projection of 2-D position xy onto the PO spline.
% Uses a coarse discrete search over Xpo knots to identify the correct basin,
% then refines with fminbnd over a ±5-knot window around the discrete seed.
% This gives sub-knot accuracy regardless of orbit period.
    % Coarse seed (nearest of N_po discrete knots)
    d        = hypot(S.Xpo(:,1) - xy(1), S.Xpo(:,2) - xy(2));
    [~, idx] = min(d);
    t0       = S.t_dense(idx);
    % Refine over ±5-knot window (narrow: orbit is smooth, correct basin found above)
    dt  = S.Tf_PO / (numel(S.t_dense) - 1);
    lo  = max(0,        t0 - 5*dt);
    hi  = min(S.Tf_PO,  t0 + 5*dt);
    xy  = xy(:);    % ensure 2×1 column — matches ppval(S.pp_xy, scalar) output
    alpha = fminbnd(@(a) sum((ppval(S.pp_xy, a) - xy).^2), lo, hi);
end

% -------------------------------------------------------------------------

function alpha = local_find_alpha(S, xy)
% Discrete nearest-knot search (kept as internal fallback).
    d = hypot(S.Xpo(:,1) - xy(1), S.Xpo(:,2) - xy(2));
    [~, idx] = min(d);
    alpha = S.t_dense(idx);
end

% -------------------------------------------------------------------------

function s = local_conv_str(converged)
    if converged, s = 'CONVERGED'; else, s = 'NOT CONVERGED'; end
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
