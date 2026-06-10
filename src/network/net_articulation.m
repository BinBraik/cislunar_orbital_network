function is_ap = net_articulation(A)
%NET_ARTICULATION  Identify articulation points by sequential node removal.
%
%   Step 9 of the network centrality algorithm.
%
%   Physical meaning: an articulation point is an orbit whose removal
%   disconnects the feasible-transfer network — some orbit pairs lose ALL
%   feasible routes (not just more expensive ones, but zero routes).
%
%   Implementation: For small N (= 13) the simplest correct approach is
%   O(N²) node-removal + BFS.  A node v is an articulation point iff
%   removing it strictly increases the number of connected components.
%
% Inputs
%   A     [N×N double]  unweighted undirected adjacency matrix (0/1)
%
% Outputs
%   is_ap [N×1 logical] true if node i is an articulation point

N = size(A, 1);
is_ap = false(N, 1);

if N <= 1
    return
end

% Baseline component count on the full graph
n_comp_base = i_count_components(A, N);

for rm = 1:N
    % Build subgraph without node rm
    keep  = [1:rm-1, rm+1:N];
    A_sub = A(keep, keep);
    n_sub = numel(keep);

    n_comp_sub = i_count_components(A_sub, n_sub);

    % Articulation point: removal increases component count
    if n_comp_sub > n_comp_base
        is_ap(rm) = true;
    end
end

end

% ── Local helper: BFS component counter ─────────────────────────────────────
function nc = i_count_components(A, n)
%I_COUNT_COMPONENTS  Count connected components via BFS.

if n == 0
    nc = 0;
    return
end

visited = false(n, 1);
nc      = 0;

for s = 1:n
    if visited(s), continue; end

    nc           = nc + 1;
    visited(s)   = true;
    queue        = s;
    qi           = 1;

    while qi <= numel(queue)
        u  = queue(qi);
        qi = qi + 1;

        nbrs = find(A(u,:) | A(:,u)');
        for v = nbrs
            if ~visited(v)
                visited(v)   = true;
                queue(end+1) = v;   %#ok<AGROW>
            end
        end
    end
end

end
