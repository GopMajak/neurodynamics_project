%% Synthetic test: does distance-dependent zero-inflation distort the EDR fit?
% Small synthetic point cloud (NOT real data) with a KNOWN true decay rate,
% where zero-inflation probability increases with distance (mimicking real
% tractography sparsity). Checks that compute_distance_bin_stats_v2
% correctly recovers the true lambda from nonzero-only means, while the
% all-pairs (zeros included) mean is distorted -- the hypothesized bug.

clear; clc;
this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(fullfile(code_dir, 'utils'));

rng(7);
n = 500;
true_A = 5.0;
true_lambda = 0.05;

theta = acos(2*rand(n,1) - 1);
phi = 2*pi*rand(n,1);
R = 70 + 3*randn(n,1);
vertices = [R.*sin(theta).*cos(phi), R.*sin(theta).*sin(phi), R.*cos(theta)];

D_true = squareform(pdist(vertices));
true_weight = true_A * exp(-true_lambda * D_true);

% distance-dependent zero-inflation: detection probability decays with
% distance (mimics tractography failing to detect long-range fibers)
detect_prob = exp(-0.02 * D_true);
detected = rand(n) < detect_prob;
SC = true_weight .* detected;
SC(1:n+1:end) = 0;

stats = compute_distance_bin_stats_v2(vertices, SC, 30, [0, max(D_true(:))], 100);

[A_all, lambda_all] = fit_edr_exponential(stats.centers, stats.mean_all, 2);
[A_nz, lambda_nz] = fit_edr_exponential(stats.centers, stats.mean_nonzero, 2);

fprintf('True params:              A = %.4f, lambda = %.4f\n', true_A, true_lambda);
fprintf('All-pairs fit (zeros in): A = %.4f, lambda = %.4f\n', A_all, lambda_all);
fprintf('Nonzero-only fit:         A = %.4f, lambda = %.4f\n', A_nz, lambda_nz);

err_all = abs(lambda_all - true_lambda) / true_lambda;
err_nz = abs(lambda_nz - true_lambda) / true_lambda;
fprintf('\nRelative lambda error -- all-pairs: %.1f%%, nonzero-only: %.1f%%\n', ...
    err_all*100, err_nz*100);

if err_nz < err_all
    fprintf('CONFIRMED: nonzero-only fit recovers the true decay rate much better\n');
    fprintf('           than the all-pairs (zero-inflated) fit -- matches the\n');
    fprintf('           hypothesized cause of Stage 1''s bad lambda.\n');
else
    fprintf('UNEXPECTED: nonzero-only fit did not improve on all-pairs fit --\n');
    fprintf('            zero-inflation may not be the (sole) explanation.\n');
end
