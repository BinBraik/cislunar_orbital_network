function [dist, pred] = net_floyd_warshall(W)
%NET_FLOYD_WARSHALL  All-pairs shortest-path distances (Floyd-Warshall).
%
%   Step 3 of the network centrality algorithm.
%
%   CRITICAL: Only FEASIBLE edges appear in W (infeasible edges are Inf).
%   An infeasible edge CANNOT be used as an intermediate hop, so the
%   multi-hop path cost is the sum of individual feasible-edge DV costs.
%
% Inputs
%   W    [N×N double]  DV-weighted adjacency matrix.
%                      W(i,j) = DV cost (m/s) for direct feasible edge.
%                      W(i,j) = Inf  for infeasible / missing edge.
%                      W(i,i) = 0.
%
% Outputs
%   dist [N×N double]  All-pairs minimum DV cost.
%                      dist(i,j) = Inf means no feasible-edge path exists.
%                      dist(i,i) = 0.
%   pred [N×N uint8]   Next-hop matrix for path reconstruction (optional).
%                      pred(i,j) = first node after i on optimal path to j.
%                      pred(i,j) = 0 if j is unreachable from i or i==j.
%                      Pass to net_path_reconstruct(pred, src, dst) to get
%                      the full hop sequence [src, h1, ..., dst].

N    = size(W, 1);
dist = W;

% ── Initialise predecessor matrix (only if caller requests it) ───────────────
compute_pred = (nargout > 1);
if compute_pred
    % pred(i,j) = j  when a direct feasible edge i→j exists
    % pred(i,j) = 0  when no edge (unreachable or diagonal)
    pred = uint8(repmat(uint8(1:N), N, 1));  % pred(i,j) = j for all j
    pred(~isfinite(W)) = 0;                  % no direct edge → unreachable
    pred(1:N+1:end)    = 0;                  % diagonal → no hop
end

% ── Floyd-Warshall ────────────────────────────────────────────────────────────
for k = 1:N
    % Vectorised distance update
    new_dist = dist(:,k) + dist(k,:);   % [N×N] via broadcast
    improved = new_dist < dist;

    dist(improved) = new_dist(improved);

    % Predecessor update: route i→j via k means first hop from i is pred(i,k)
    if compute_pred && any(improved(:))
        pred_via_k         = repmat(pred(:,k), 1, N);  % (i,j) entry = pred(i,k)
        pred(improved)     = pred_via_k(improved);
    end
end

% Guarantee exact zeros on diagonal (floating-point safety)
dist(1:N+1:end) = 0;

end
