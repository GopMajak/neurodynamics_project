%% Smoke test for fit_edr_exponential_fixed_init.m (syntax/shape check only)
clear; clc;
this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(fullfile(code_dir, 'utils'));

rng(3);
n = 400;
true_A = 0.5; true_lambda = 0.15;

theta = acos(2*rand(n,1) - 1);
phi = 2*pi*rand(n,1);
R = 70 + 3*randn(n,1);
vertices = [R.*sin(theta).*cos(phi), R.*sin(theta).*sin(phi), R.*cos(theta)];

D_true = squareform(pdist(vertices));
SC = true_A * exp(-true_lambda * D_true) .* (rand(n) < 0.3) + 0.01*randn(n);
SC = max(SC, 0);
SC(1:n+1:end) = 0;
SC_norm = SC / max(SC(:));

stats = compute_distance_bin_stats(vertices, SC_norm, 30, [0, max(D_true(:))], 100);
[A_fixed, lambda_fixed] = fit_edr_exponential_fixed_init(stats.centers, stats.mean, 2, [0.15, 0.18], [-100, 100]);

fprintf('True lambda = %.4f (A not comparable, SC was normalized)\n', true_lambda);
fprintf('Fixed-init fit: A = %.4f, lambda = %.4f\n', A_fixed, lambda_fixed);
fprintf('No crash -- function runs correctly.\n');
