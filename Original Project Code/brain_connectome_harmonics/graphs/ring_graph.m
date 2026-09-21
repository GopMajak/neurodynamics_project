function A = ring_graph( N, k )
% Builds a ring network of N nodes, where each node connects to its k
% nearest neighbors on either side.
%
% Inputs:
%   N - number of nodes
%   k - number of neighbors connected on each side
%
% Output:
%   A - adjacency matrix (N x N, unweighted)

A = false( N, N );
m = -k:k; m( m == 0 ) = [];

for ii = 0:(N-1)
    for jj = m
        A( ii+1, mod(ii+jj,N)+1 ) = 1;
    end
end
