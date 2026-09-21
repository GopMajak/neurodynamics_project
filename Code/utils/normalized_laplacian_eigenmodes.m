function [eigenvectors, eigenvalues] = normalized_laplacian_eigenmodes(W, num_modes)
%NORMALIZED_LAPLACIAN_EIGENMODES Solve the normalized graph Laplacian
%   eigenproblem (Vohryzek et al. 2025 eq. for L^norm; Pang et al. 2023
%   eq. 6-7): L = D - A, L_norm = D^{-1/2} L D^{-1/2}, solve
%   L_norm * psi_k = lambda_k * psi_k for the NUM_MODES smallest
%   eigenvalues/eigenvectors.
%
%   W : N x N (sparse) weighted adjacency matrix. Symmetrized internally
%       as (W + W')/2 to match Pang et al. 2023 eq. 7.

W = (W + W') / 2;
d = full(sum(W, 2));
d(d == 0) = eps; % avoid divide-by-zero for isolated nodes

n = numel(d);
D_inv_sqrt = spdiags(1 ./ sqrt(d), 0, n, n);
L = spdiags(d, 0, n, n) - W;
L_norm = D_inv_sqrt * L * D_inv_sqrt;
L_norm = (L_norm + L_norm') / 2;

opts.tol = 1e-10;
opts.maxit = 5000;
[eigenvectors, eigenvalues_diag] = eigs(L_norm, num_modes, 'smallestabs', opts);

eigenvalues = diag(eigenvalues_diag);
[eigenvalues, order] = sort(eigenvalues, 'ascend');
eigenvectors = eigenvectors(:, order);
end
