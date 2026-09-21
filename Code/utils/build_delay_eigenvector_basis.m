function [psi_real, eigvals] = build_delay_eigenvector_basis(W, tau, omega, kappa, dt)
%BUILD_DELAY_EIGENVECTOR_BASIS Real orthonormal basis from the complex
%   delay-rotated coupling operator's eigenvectors.
%
%   W_delay(i,j) = kappa * exp(-1i*omega(i)*dt*tau(i,j)) * W(i,j)
%   (matches km_dynamics_hcp.m's construction). Its eigenvectors are
%   complex and, since W_delay is not Hermitian in general, not
%   orthogonal -- so they cannot be used directly as a projection basis
%   for real-valued reconstruction targets (this is exactly the
%   obstacle the thesis proposal's Materials/Methods flags).
%
%   This builds a real, ORTHONORMAL basis instead by taking
%   [real(v_1), imag(v_1), real(v_2), imag(v_2), ...] for the eigenvectors
%   v_k sorted by descending real part of eigenvalue (matching
%   eig_sorted.m's convention -- most dominant/theoretically important
%   mode first), then orthonormalizing via economy QR. QR naturally caps
%   the result at N columns and (without column pivoting) preserves the
%   input order, so mode 1 of the returned basis is built from the
%   dominant eigenvector, mode 2 refines it, etc. -- an ordered,
%   real-valued analogue of the (non-orthogonal) complex eigenvectors.
%
%   No simulation needed -- this is pure linear algebra on W_delay.
%
%   Inputs: W (N x N real weighted graph, zero diagonal), tau (N x N
%   delay matrix, timesteps), omega (N x 1, rad/s), kappa (scalar), dt (s)
%
%   Outputs: psi_real (N x N real orthonormal basis, ordered),
%            eigvals (N x 1 complex eigenvalues of W_delay, same order
%            used to build the candidate directions)

N = size(W, 1);
W_delay = kappa .* (exp((-1i * omega * dt) .* tau) .* W);

[V, D] = eig(double(W_delay));
eigvals_all = diag(D);
[eigvals, order] = sort(real(eigvals_all), 'descend');
V = V(:, order);
eigvals = eigvals_all(order); % keep the actual complex eigenvalues, sorted by real part

candidates = zeros(N, 2*N);
candidates(:, 1:2:end) = real(V);
candidates(:, 2:2:end) = imag(V);

[Q, ~] = qr(candidates, 0); % economy QR, no column pivoting -> preserves input order
psi_real = Q;
end
