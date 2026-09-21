%% Synthetic smoke test for build_anatomical_graphs_blocked_v2.m
clear; clc;
this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(fullfile(code_dir, 'utils'));

rng(5);
n = 400;
true_A = 0.5; true_lambda = 0.05;

theta = acos(2*rand(n,1) - 1);
phi = 2*pi*rand(n,1);
R = 70 + 3*randn(n,1);
vertices = [R.*sin(theta).*cos(phi), R.*sin(theta).*sin(phi), R.*cos(theta)];

D_true = squareform(pdist(vertices));
detect_prob = exp(-0.02 * D_true);
SC = (true_A * exp(-true_lambda * D_true)) .* (rand(n) < detect_prob);
SC(1:n+1:end) = 0;

% inject deliberate LR exceptions
n_outliers = 10;
far_i = randi(n, n_outliers, 1); far_j = randi(n, n_outliers, 1);
for k = 1:n_outliers
    i = far_i(k); j = far_j(k);
    if D_true(i,j) > 60
        SC(i,j) = 3; SC(j,i) = 3;
    end
end

max_sc = max(SC(:));
C_norm = SC / max_sc;

n_bins = 30;
max_dist = max(D_true(:));
stats = compute_distance_bin_stats(vertices, C_norm, n_bins, [0, max_dist], 100);
[A_fit, lambda_fit] = fit_edr_exponential_fixed_init(stats.centers, stats.mean, 2, [0.15, 0.18], [-100, 100]);
fprintf('Fitted (normalized data): A=%.4f, lambda=%.4f\n', A_fit, lambda_fit);

[A_bin, A_cont, A_lr, report] = build_anatomical_graphs_blocked_v2( ...
    vertices, C_norm, 0.05, A_fit, lambda_fit, stats.edges, stats.mean, stats.std, ...
    2, n_bins-1, 3, 60, 1e-6, 100, 42);

assert(issymmetric(A_bin), 'EDR binary must be symmetric');
assert(issymmetric(A_cont), 'EDR continuous must be symmetric');
assert(issymmetric(A_lr), 'EDR+LR must be symmetric');
assert(report.n_lr_exceptions > 0, 'expected to detect injected LR outliers');
fprintf('LR exceptions detected: %d (%.2f%% of all pairs, injected %d candidates)\n', ...
    report.n_lr_exceptions, report.pct_lr_exceptions, n_outliers);
fprintf('Densities -- binary: %.2f%%, continuous: %.2f%%, EDR+LR: %.2f%%\n', ...
    report.density_binary*100, report.density_continuous*100, report.density_edrlr*100);
fprintf('No crash -- all assertions passed.\n');
