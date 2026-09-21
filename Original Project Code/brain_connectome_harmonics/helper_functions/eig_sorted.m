function [vs,ds] = eig_sorted( a )
% Computes the eigenvectors/eigenvalues of a matrix and sorts them from
% largest to smallest eigenvalue (MATLAB's eig() does not guarantee order).
%
% Input:
%   a - square matrix
%
% Outputs:
%   vs - eigenvectors, sorted to match ds
%   ds - eigenvalues, sorted descending

[v,d] = eig( double(a) );
[ds,ii_d] = sort( real(diag(d)), 'descend' );
vs = nan( size(a,1) );
for ii = 1:size(a,1)
    vs(:,ii) = (v(:,ii_d(ii)));
end
