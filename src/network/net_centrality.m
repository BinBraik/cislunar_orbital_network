function metrics = net_centrality(A, W, D_sym, dist, DVcap_true_mps)
%NET_CENTRALITY  Compute node centrality metrics for one snapshot (Steps 5-9).
%
%   All metrics are computed on the feasible-edge topology established by
%   net_build_graph; infeasible edges are completely excluded from all
%   path computations.
%
% Inputs
%   A              [N×N double]  unweighted feasibility adjacency (0/1)
%   W              [N×N double]  DV-weighted adjacency (Inf infeasible, 0 diag)
%   D_sym          [N×N double]  symmetrised DV matrix (Inf for NaN/missing, 0 diag)
%   dist           [N×N double]  all-pairs shortest paths from net_floyd_warshall
%   DVcap_true_mps scalar        physical DV budget (m/s)
%
% Outputs — struct with fields:
%   R_budget           [N×1]  fraction of budget-reachable destinations
%   harmonic_closeness [N×1]  H(A) = 1/(N-1) * sum_{B!=A} 1/d(A,B)
%   betweenness        [N×1]  B(A) = 2/((N-1)(N-2)) * sum_{S<T,S,T!=A} sigma_ST(A)/sigma_ST
%   strength           [N×1]  S(A) = 1/(N-1) * sum_{B in Gamma(A)} 1/w(A,B)
%   is_articulation   [N×1 logical]  articulation-point flag
%   budget_reach      [N×N logical]  pair-level budget-reachability matrix

N       = size(A, 1);
tol_rel = 1e-6;             % relative tolerance for path-cost equality

% ── Step 5: Budgeted reachability per node ───────────────────────────────────
% Pair (i,j) is budget-reachable iff shortest-path cost <= DVcap AND finite.
budget_reach = isfinite(dist) & (dist <= DVcap_true_mps);
budget_reach(1:N+1:end) = false;          % exclude self

R_budget = sum(budget_reach, 2) / (N - 1);

% ── Step 6: Harmonic closeness ───────────────────────────────────────────────
% H(A) = 1/(N-1) * sum_{B != A} 1/d(A,B)
% Pairs with no admissible path have d = inf, contributing 1/inf = 0.
harmonic_closeness = zeros(N, 1);
for i = 1:N
    dests = find(budget_reach(i, :));
    if isempty(dests), continue; end
    harmonic_closeness(i) = sum(1 ./ dist(i, dests)) / (N - 1);
end

% ── Step 7: Betweenness ──────────────────────────────────────────────────────
% B(A) = 2/((N-1)(N-2)) * sum_{S<T, S,T != A} sigma_ST(A)/sigma_ST
%
% Accumulated using ordered pairs (equivalent since network is undirected;
% ordered sum = 2 × unordered sum, and the factor 2 cancels with the
% 1/(N-1)/(N-2) normalization to give B(A) = ordered_sum / ((N-1)*(N-2))).
%
% Normalisation denominator is FIXED at (N-1)*(N-2), independent of budget.
% Pairs with no admissible path contribute 0.
%
% Pre-compute sigma(s,v) = number of DV-shortest paths from s to v.
sigma = zeros(N, N);
for s = 1:N
    sigma(s, s) = 1;
    [~, order] = sort(dist(s, :));   % process in increasing dist from s
    for vi = 1:N
        v = order(vi);
        if v == s,                  continue; end
        if ~isfinite(dist(s, v)),   continue; end

        tol_v = tol_rel * (1 + dist(s, v));
        for u = 1:N
            if u == v,                continue; end
            if ~isfinite(W(u, v)),    continue; end   % must be direct edge
            if ~isfinite(dist(s, u)), continue; end
            if abs(dist(s, u) + W(u, v) - dist(s, v)) <= tol_v
                sigma(s, v) = sigma(s, v) + sigma(s, u);
            end
        end
    end
end

% Accumulate betweenness via ordered-pair loop
betweenness = zeros(N, 1);

for s = 1:N
    for t = 1:N
        if s == t || ~budget_reach(s, t), continue; end
        if sigma(s, t) == 0, continue; end

        d_st    = dist(s, t);
        tol_st  = tol_rel * (1 + d_st);

        for i = 1:N
            if i == s || i == t,       continue; end
            if ~isfinite(dist(s, i)) || ~isfinite(dist(i, t)), continue; end
            if abs(dist(s, i) + dist(i, t) - d_st) <= tol_st
                betweenness(i) = betweenness(i) + ...
                    (sigma(s, i) * sigma(i, t)) / sigma(s, t);
            end
        end
    end
end

% Fixed normalisation: ordered_sum / ((N-1)*(N-2))  [equivalent to manuscript]
if N >= 3
    betweenness = betweenness / ((N - 1) * (N - 2));
end

% ── Step 8: Strength ─────────────────────────────────────────────────────────
% S(A) = 1/(N-1) * sum_{B in Gamma(A)} 1/w(A,B)
% Uses raw (symmetrised) DV for direct feasible edges; normalised by (N-1).
strength = zeros(N, 1);
for i = 1:N
    nbrs = find(A(i, :));
    if isempty(nbrs), continue; end
    strength(i) = sum(1 ./ D_sym(i, nbrs)) / (N - 1);
end

% ── Step 9: Articulation points ──────────────────────────────────────────────
is_articulation = net_articulation(A);

% ── Pack output ──────────────────────────────────────────────────────────────
metrics.R_budget           = R_budget;
metrics.harmonic_closeness = harmonic_closeness;
metrics.betweenness        = betweenness;
metrics.strength           = strength;
metrics.is_articulation    = is_articulation;
metrics.budget_reach       = budget_reach;

end
