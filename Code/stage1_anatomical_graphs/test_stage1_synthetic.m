%% Synthetic smoke test for Stage 1 utility functions.
% Uses a small synthetic point cloud + synthetic connectome (NOT the real
% 29,696-vertex HCP data) to check that every function runs, returns
% correctly shaped/sane outputs, before running the full-scale pipeline
% on the server.

clear; clc;
this_file = mfilename('fullpath');
code_dir = fileparts(fileparts(this_file));
addpath(fullfile(code_dir, 'utils'));

rng(1);
n = 400;

%% Synthetic cortical-like point cloud (points on a bumpy sphere, ~70mm radius)
theta = acos(2*rand(n,1) - 1);
phi = 2*pi*rand(n,1);
R = 70 + 3*randn(n,1);
vertices = [R.*sin(theta).*cos(phi), R.*sin(theta).*sin(phi), R.*cos(theta)];

%% Synthetic structural connectome: EDR-like decay + a few strong LR outliers
D_true = squareform(pdist(vertices));
SC = 5*exp(-0.05*D_true) + 0.05*randn(n);
SC = max(SC, 0);
SC(1:n+1:end) = 0; % zero diagonal

% inject deliberate long-range exceptions: strong weight at large distance
n_outliers = 15;
far_pairs_i = randi(n, n_outliers, 1);
far_pairs_j = randi(n, n_outliers, 1);
for k = 1:n_outliers
    i = far_pairs_i(k); j = far_pairs_j(k);
    if D_true(i,j) > 60
        SC(i,j) = 10; SC(j,i) = 10; % force a strong outlier weight
    end
end

fprintf('=== Test 1: compute_distance_bin_stats ===\n');
dist_range = [0, max(D_true(:))*1.01];
stats = compute_distance_bin_stats(vertices, SC, 20, dist_range, 100);
assert(numel(stats.centers) == 20, 'expected 20 bin centers');
assert(all(stats.count >= 0), 'bin counts must be non-negative');
fprintf('  OK -- %d bins, total pairs counted = %d (expected %d)\n', ...
    numel(stats.centers), sum(stats.count), n*(n-1)/2);
assert(sum(stats.count) == n*(n-1)/2, 'pair count mismatch (upper triangle)');

fprintf('=== Test 2: fit_edr_exponential ===\n');
fit_stats = compute_distance_bin_stats(vertices, SC, 30, [0, max(D_true(:))], 100);
[A_fit, lambda_fit] = fit_edr_exponential(fit_stats.centers, fit_stats.mean, 2);
fprintf('  Fitted A = %.4f, lambda = %.4f (true-ish A~5, lambda~0.05)\n', A_fit, lambda_fit);
assert(A_fit > 0 && lambda_fit > 0, 'fitted parameters must be positive');

fprintf('=== Test 3: build_anatomical_graphs_blocked ===\n');
lr_stats = compute_distance_bin_stats(vertices, SC, 20, [0, max(D_true(:))*1.01], 100);
[A_bin, A_cont, A_lr, report] = build_anatomical_graphs_blocked( ...
    vertices, SC, 0.05, lambda_fit, lr_stats, 3, 60, 1e-6, 100, 42);

assert(issymmetric(A_bin), 'EDR binary must be symmetric');
assert(issymmetric(A_cont), 'EDR continuous must be symmetric');
assert(issymmetric(A_lr), 'EDR+LR must be symmetric');
assert(nnz(A_bin) > 0, 'EDR binary graph is empty');
assert(nnz(A_cont) > 0, 'EDR continuous graph is empty');
assert(nnz(A_lr) >= nnz(A_cont), 'EDR+LR should have at least as many edges as its continuous base');
assert(report.n_lr_exceptions > 0, 'expected to detect the injected LR outliers');
fprintf('  OK -- LR exceptions detected: %d (injected %d candidates)\n', report.n_lr_exceptions, n_outliers);
fprintf('  Densities -- binary: %.2f%%, continuous: %.2f%%, EDR+LR: %.2f%%\n', ...
    report.density_binary*100, report.density_continuous*100, report.density_edrlr*100);

fprintf('=== Test 4: normalized_laplacian_eigenmodes ===\n');
num_modes = 20;
[eig_vec, eig_val] = normalized_laplacian_eigenmodes(A_cont, num_modes);
assert(isequal(size(eig_vec), [n, num_modes]), 'eigenvector matrix shape mismatch');
assert(numel(eig_val) == num_modes, 'eigenvalue count mismatch');
assert(all(diff(eig_val) >= -1e-8), 'eigenvalues must be sorted ascending');
assert(abs(eig_val(1)) < 1e-6, 'first eigenvalue should be ~0 for a connected graph');
fprintf('  OK -- eigenvalues[1:5] = %s\n', mat2str(eig_val(1:5)', 4));

fprintf('\nAll synthetic smoke tests passed.\n');
