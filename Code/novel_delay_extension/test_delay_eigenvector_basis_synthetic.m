%% Smoke test for build_delay_eigenvector_basis.m (synthetic, small scale)
clear; clc;
this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(fullfile(code_dir, 'utils'));

rng(4);
n = 30;
theta = acos(2*rand(n,1) - 1); phi = 2*pi*rand(n,1); R = 70 + 3*randn(n,1);
vertices = [R.*sin(theta).*cos(phi), R.*sin(theta).*sin(phi), R.*cos(theta)];
D = squareform(pdist(vertices));
W = exp(-0.02*D); W(1:n+1:end) = 0;

tau = floor(D/1000/5/0.001); tau(tau<1 & D>0) = 1; % crude delay matrix, speed=5m/s, dt=0.001
omega = 2*pi*10*ones(n,1);
kappa = 6; dt = 0.001;

[psi, eigvals] = build_delay_eigenvector_basis(W, tau, omega, kappa, dt);

fprintf('psi size: %d x %d\n', size(psi,1), size(psi,2));
ortho_err = max(max(abs(psi'*psi - eye(size(psi,2)))));
fprintf('Max deviation from orthonormality (psi''*psi vs I): %.2e\n', ortho_err);
fprintf('First 5 eigenvalues (real part): %s\n', mat2str(round(real(eigvals(1:5)),3)'));
assert(ortho_err < 1e-8, 'basis is not orthonormal');
fprintf('PASS: basis is orthonormal, correct shape, no crash.\n');
