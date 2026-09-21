%% Synthetic test: do rare strong long-range outliers bias a per-bin MEAN
%  fit, and does per-bin MEDIAN recover the true decay rate instead?

clear; clc;
this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(fullfile(code_dir, 'utils'));

rng(11);
n = 500;
true_A = 5.0;
true_lambda = 0.05;

theta = acos(2*rand(n,1) - 1);
phi = 2*pi*rand(n,1);
R = 70 + 3*randn(n,1);
vertices = [R.*sin(theta).*cos(phi), R.*sin(theta).*sin(phi), R.*cos(theta)];

D_true = squareform(pdist(vertices));

% background: sparse detection, weight follows true EDR decay
detect_prob = exp(-0.02 * D_true);
detected = rand(n) < detect_prob;
SC = (true_A * exp(-true_lambda * D_true)) .* detected;

% inject rare "exception" connections: 0.5% of long-distance (>60mm) pairs
% get an anomalously strong weight unrelated to distance decay
long_range_mask = triu(D_true > 60, 1);
candidates = find(long_range_mask);
n_exceptions = round(0.02 * numel(candidates));
exception_idx = candidates(randperm(numel(candidates), n_exceptions));
[ei, ej] = ind2sub(size(SC), exception_idx);
for k = 1:numel(ei)
    strong = true_A * 8; % ~8x stronger than the local baseline amplitude
    SC(ei(k), ej(k)) = strong;
    SC(ej(k), ei(k)) = strong;
end
SC(1:n+1:end) = 0;

n_bins = 30;
stats_mean = compute_distance_bin_stats_v2(vertices, SC, n_bins, [0, max(D_true(:))], 100);
stats_robust = compute_distance_bin_median(vertices, SC, n_bins, [0, max(D_true(:))], 100);

[A_mean, lambda_mean] = fit_edr_exponential(stats_mean.centers, stats_mean.mean_nonzero, 1);
[A_med, lambda_med] = fit_edr_exponential(stats_robust.centers, stats_robust.median_nonzero, 1);
[A_trim, lambda_trim] = fit_edr_exponential(stats_robust.centers, stats_robust.trimmed_mean_nonzero, 1);

fprintf('True params:        A = %.4f, lambda = %.4f\n', true_A, true_lambda);
fprintf('Nonzero MEAN fit:   A = %.4f, lambda = %.4f (%.1f%% error)\n', ...
    A_mean, lambda_mean, 100*abs(lambda_mean-true_lambda)/true_lambda);
fprintf('Nonzero MEDIAN fit: A = %.4f, lambda = %.4f (%.1f%% error)\n', ...
    A_med, lambda_med, 100*abs(lambda_med-true_lambda)/true_lambda);
fprintf('Trimmed mean fit:   A = %.4f, lambda = %.4f (%.1f%% error)\n', ...
    A_trim, lambda_trim, 100*abs(lambda_trim-true_lambda)/true_lambda);
