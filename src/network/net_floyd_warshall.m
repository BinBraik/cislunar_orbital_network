function dist = net_floyd_warshall(W)
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

N    = size(W, 1);
dist = W;

for k = 1:N
    % Vectorised update: dist(i,j) = min( dist(i,j),  dist(i,k)+dist(k,j) )
    % dist(:,k) is [N×1], dist(k,:) is [1×N] → broadcast to [N×N]
    dist = min(dist, dist(:,k) + dist(k,:));
end

% Guarantee exact zeros on diagonal (floating-point safety)
dist(1:N+1:end) = 0;

end
