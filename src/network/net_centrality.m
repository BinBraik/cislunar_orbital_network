function metrics = net_centrality(A, W, D_sym, dist, DVcap_true_mps, VU_mps, Rmin_guard)
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
%   VU_mps         scalar        velocity unit (m/s)
%   Rmin_guard     scalar        min R_budget to qualify for reach winner
%                                (stored in output; caller applies it in tie logic)
%
% Outputs — struct with fields:
%   R_budget       [N×1]  fraction of budget-reachable destinations
%   reach          [N×1]  normalised sum of inverse shortest-path distances
%   gateway        [N×1]  DV-betweenness over budget-reachable OD pairs
%   hub            [N×1]  sum of inverse direct-edge DV costs
%   is_articulation [N×1 logical]  articulation-point flag
%   budget_reach   [N×N logical]  pair-level budget-reachability matrix
%   Rmin_guard     scalar  (passed through for downstream use)

N       = size(A, 1);
eps_dv  = 1e-3 * VU_mps;   % ~1 m/s floor; prevents hub division by zero
tol_rel = 1e-6;             % relative tolerance for path-cost equality

% ── Step 5: Budgeted reachability per node ───────────────────────────────────
% Pair (i,j) is budget-reachable iff shortest-path cost <= DVcap AND finite.
budget_reach = isfinite(dist) & (dist <= DVcap_true_mps);
budget_reach(1:N+1:end) = false;          % exclude self

R_budget = sum(budget_reach, 2) / (N - 1);

% ── Step 6: REACH metric ─────────────────────────────────────────────────────
% reach(i) = sum_{j budget-reachable} [ 1/dist(i,j) ] / (N-1)
% Rewards nodes that sit close (in DV) to many others.
reach = zeros(N, 1);
for i = 1:N
    dests = find(budget_reach(i, :));
    if isempty(dests), continue; end
    reach(i) = sum(1 ./ dist(i, dests)) / (N - 1);
end

% ── Step 7: GATEWAY metric ───────────────────────────────────────────────────
% DV-weighted betweenness centrality over budget-reachable OD pairs only.
% gateway(i) = (sum over valid (s,t) of fraction of shortest paths through i)
%              / max(1, number_of_valid_OD_pairs)
%
% First pre-compute sigma(s,v) = number of DV-shortest paths from s to v.
% Process nodes in increasing distance from s; predecessor u satisfies
%   isfinite(W(u,v))  [direct edge]  AND  dist(s,u)+W(u,v) ≈ dist(s,v).

sigma = zeros(N, N);          % sigma(s, v)
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

% Accumulate betweenness
gateway     = zeros(N, 1);
valid_pairs = 0;

for s = 1:N
    for t = 1:N
        if s == t || ~budget_reach(s, t), continue; end
        valid_pairs = valid_pairs + 1;

        if sigma(s, t) == 0, continue; end

        d_st    = dist(s, t);
        tol_st  = tol_rel * (1 + d_st);

        for i = 1:N
            if i == s || i == t,       continue; end
            if ~isfinite(dist(s, i)) || ~isfinite(dist(i, t)), continue; end
            if abs(dist(s, i) + dist(i, t) - d_st) <= tol_st
                % Fraction of s→t shortest paths passing through i:
                %   sigma(s,i) * sigma(i,t) / sigma(s,t)
                gateway(i) = gateway(i) + ...
                    (sigma(s, i) * sigma(i, t)) / sigma(s, t);
            end
        end
    end
end

if valid_pairs > 0
    gateway = gateway / valid_pairs;
end

% ── Step 8: HUB metric ───────────────────────────────────────────────────────
% hub(i) = sum_{j: direct feasible edge} [ 1 / (D(i,j) + eps) ]
% Uses raw (symmetrised) DV, not shortest-path distances.
hub = zeros(N, 1);
for i = 1:N
    nbrs = find(A(i, :));
    if isempty(nbrs), continue; end
    hub(i) = sum(1 ./ (D_sym(i, nbrs) + eps_dv));
end

% ── Step 9: Articulation points ──────────────────────────────────────────────
is_articulation = net_articulation(A);

% ── Pack output ──────────────────────────────────────────────────────────────
metrics.R_budget        = R_budget;
metrics.reach           = reach;
metrics.gateway         = gateway;
metrics.hub             = hub;
metrics.is_articulation = is_articulation;
metrics.budget_reach    = budget_reach;
metrics.Rmin_guard      = Rmin_guard;

end
