function A = distance_dependent_graph(N, alpha)
% Builds a ring of N nodes where connection strength falls off with
% distance around the ring, following a power law (1/distance^alpha).
% Weights are normalized so each node's connections sum to the same total.
%
% Inputs:
%   N     - number of nodes
%   alpha - power-law exponent (bigger = strength drops off faster)
%
% Output:
%   A - weighted adjacency matrix

A = zeros( N, N );
d = nan( N );
for ii = 1:N, for jj = 1:N
        if (ii==jj), continue; end; dist = abs( ii - jj ); d(ii,jj) = min( dist, N - dist );
end; end
eta = nansum( 1 ./ d(1,:).^alpha ); A = (1.0/eta) * (1 ./ d.^alpha);
for ii = 1:N, A(ii,ii) = 0; end
