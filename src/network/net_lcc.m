function [lcc_size, lcc_full, comp_labels] = net_lcc(A)
%NET_LCC  Largest Connected Component (LCC) of an undirected graph.
%
%   Step 4 of the network centrality algorithm.
%   Uses the UNWEIGHTED adjacency A for connectivity analysis.
%
% Inputs
%   A           [N×N double]  unweighted adjacency matrix (symmetric, 0/1)
%
% Outputs
%   lcc_size    scalar        number of nodes in the largest component
%   lcc_full    scalar        1 if lcc_size == N (fully connected), else 0
%   comp_labels [N×1 double]  component index per node (1-based integers)

N = size(A, 1);

if N == 0
    lcc_size    = 0;
    lcc_full    = 0;
    comp_labels = zeros(0,1);
    return
end

comp_labels = zeros(N, 1);
comp_id     = 0;

for start_node = 1:N
    if comp_labels(start_node) ~= 0
        continue                % already visited
    end

    comp_id               = comp_id + 1;
    comp_labels(start_node) = comp_id;

    % BFS queue (pre-allocated with room to grow)
    queue = start_node;
    qi    = 1;

    while qi <= numel(queue)
        u = queue(qi);
        qi = qi + 1;

        % Undirected: check both row and column (A is symmetric but guard anyway)
        nbrs = find(A(u,:) | A(:,u)');

        for v = nbrs
            if comp_labels(v) == 0
                comp_labels(v) = comp_id;
                queue(end+1)   = v;        %#ok<AGROW>
            end
        end
    end
end

% Tally component sizes
if comp_id == 0
    lcc_size = 0;
else
    comp_sizes = accumarray(comp_labels, ones(N,1), [comp_id, 1]);
    lcc_size   = max(comp_sizes);
end

lcc_full = double(lcc_size == N);

end
