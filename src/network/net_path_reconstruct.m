function path = net_path_reconstruct(pred, src, dst)
%NET_PATH_RECONSTRUCT  Recover hop sequence from Floyd-Warshall predecessor matrix.
%
%   path = net_path_reconstruct(pred, src, dst)
%
%   Returns the full sequence of node indices  [src, h1, h2, ..., dst]
%   representing the FW-shortest path from src to dst.
%   Returns [] if no path exists (pred(src,dst) == 0 or src == dst).
%
% Inputs
%   pred  [N×N uint8]  Next-hop matrix from net_floyd_warshall (2nd output).
%                      pred(i,j) = first node to visit after i on way to j.
%                      pred(i,j) = 0 means j is unreachable from i.
%   src   scalar       Source node index (1-based).
%   dst   scalar       Destination node index (1-based).
%
% Output
%   path  [1×K int32]  Node sequence [src, ..., dst].  K >= 2 for a valid path.
%                      [] if unreachable or src == dst.

if src == dst || pred(src, dst) == 0
    path = [];
    return
end

N       = size(pred, 1);
path    = int32(src);
current = src;

for iter = 1:N   % a cycle-free path visits each node at most once
    nxt = double(pred(current, dst));
    if nxt == 0
        % Path broke — unreachable from current to dst
        path = [];
        return
    end
    path    = [path, int32(nxt)];  %#ok<AGROW>
    current = nxt;
    if current == dst
        return
    end
end

% If we reach here the path exceeded N hops — cycle guard triggered
warning('net_path_reconstruct: path from %d to %d exceeded %d hops (possible cycle).', ...
    src, dst, N);
path = [];

end
