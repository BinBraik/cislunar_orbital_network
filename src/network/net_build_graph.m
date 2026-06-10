function [A, W, D_sym, T_sym, edges_kept, DVcap_true_mps, Tmax_true_days, skip] = ...
    net_build_graph(DV_raw, TOF_raw, DVcap_nd, Tmax_nd, VU_mps, TU_days, budget_factor)
%NET_BUILD_GRAPH  Steps 1-2: sanitise one snapshot and apply feasibility filter.
%
% Inputs
%   DV_raw       [N×N double]  raw DV matrix from DVmatrix_sweep{di,dj} (m/s)
%   TOF_raw      [N×N double]  raw TOF matrix from TOFmatrix_sweep{di,dj} (days)
%   DVcap_nd     scalar        nondimensional DV-cap for this snapshot
%   Tmax_nd      scalar        nondimensional Tmax for this snapshot
%   VU_mps       scalar        velocity unit (m/s), from cfg.units.VU_mps
%   TU_days      scalar        time unit (days),    from cfg.units.TU_days
%   budget_factor scalar       multiplier applied to convert nd budgets to physical
%                              (= 2 for two-manoeuvre / two-leg transfers)
%
% Outputs
%   A              [N×N double]  unweighted feasibility adjacency (0/1, symmetric)
%   W              [N×N double]  DV-weighted adjacency; W(i,j) = D_sym(i,j) if
%                                feasible, Inf otherwise; W(i,i) = 0
%   D_sym          [N×N double]  symmetrised DV matrix; NaN replaced by Inf (m/s)
%   T_sym          [N×N double]  symmetrised TOF matrix; NaN replaced by Inf (days)
%   edges_kept     scalar        number of feasible directed edges (upper+lower)
%   DVcap_true_mps scalar        physical DV budget = budget_factor*DVcap_nd*VU_mps
%   Tmax_true_days scalar        physical TOF budget = budget_factor*Tmax_nd*TU_days
%   skip           logical       true → snapshot is all-NaN/empty; caller should skip

% ── Early-exit check ────────────────────────────────────────────────────────
skip = isempty(DV_raw) || all(isnan(DV_raw(:)));
if skip
    A = []; W = []; D_sym = []; T_sym = [];
    edges_kept = 0;
    DVcap_true_mps  = budget_factor * DVcap_nd * VU_mps;
    Tmax_true_days  = budget_factor * Tmax_nd   * TU_days;
    return
end

N = size(DV_raw, 1);

% ── Step 1: sanitise (NaN → Inf) ────────────────────────────────────────────
D = double(DV_raw);
T = double(TOF_raw);
D(isnan(D)) = Inf;
T(isnan(T)) = Inf;

% Enforce symmetry
%   DV : take the cheaper option at each pair
D = min(D, D');
%   TOF: average the two directions
T = (T + T') / 2;

% Zero diagonal
D(1:N+1:end) = 0;
T(1:N+1:end) = 0;

D_sym = D;
T_sym = T;

% ── Step 2: feasibility filter ───────────────────────────────────────────────
DVcap_true_mps  = budget_factor * DVcap_nd * VU_mps;
Tmax_true_days  = budget_factor * Tmax_nd   * TU_days;

% Edge (i,j) feasible iff ALL four conditions hold:
%   (a) D(i,j) is finite       (b) D(i,j) <= DVcap_true_mps
%   (c) T(i,j) is finite       (d) T(i,j) <= Tmax_true_days
feas = isfinite(D) & (D <= DVcap_true_mps) & ...
       isfinite(T) & (T <= Tmax_true_days);
feas(1:N+1:end) = false;          % no self-loops

A = double(feas);

% DV-weighted adjacency (Inf for infeasible, 0 on diagonal)
W            = Inf(N);
W(feas)      = D(feas);
W(1:N+1:end) = 0;

edges_kept = sum(feas(:));

end
