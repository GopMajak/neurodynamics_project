function [v,d] = circulant_eigensystem( a )
% Computes the eigenvectors and eigenvalues of a circulant matrix (like a
% ring network's adjacency matrix) using the known closed-form solution
% (Fourier modes), instead of a general-purpose eigensolver.
%
% Input:
%   a - adjacency matrix (must actually be circulant - not checked here)
%
% Outputs:
%   v - eigenvectors (N x N)
%   d - eigenvalues (N x N diagonal matrix)

N = size( a, 1 );

v = zeros( N ); d = zeros( N );
for ii = 1:N
    for jj = 1:N
        v(ii,jj) = (1/sqrt(N)) * exp( -2*pi*1i/N * (ii-1) * (jj-1) );
    end
end

for ii = 1:N
    d(ii,ii) = sum( a(1,:) .* exp( -2*pi*1i/N * ( ii - 1 ) * ( (1:N) - 1 ) ) );
end
