%% Diagnostic v2: does a robust (median) per-bin statistic fix the EDR fit?
%
% Follow-up to diagnose_edr_fit.m. That script ruled out zero-inflation
% as the cause of Stage 1's bad lambda (nonzero-only MEAN made the fit
% worse, not better). A synthetic test (test_median_synthetic.m) then
% showed that rare, anomalously strong long-range "exception" connections
% -- exactly what Vohryzek et al. 2025 define as LR exceptions -- can
% catastrophically bias a per-bin MEAN (even a trimmed one), while a
% per-bin MEDIAN recovers the true decay rate almost exactly. This script
% checks whether that holds on the real structural connectome.
%
% NEW FILE -- does not modify diagnose_edr_fit.m, run_stage1.m, or any
% other existing file. All diagnostic attempts are preserved side by side.

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
project_root = fileparts(code_dir);
addpath(fullfile(code_dir, 'utils'));

data_dir = fullfile(project_root, 'data');
surface_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_midthickness-lh.vtk');
cortex_mask_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt');
sc_file = fullfile(data_dir, 'empirical', 'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat');

n_bins_fit = 400;
fit_dist_range = [10, 170];
exclude_bins = 25;
block_size = 2000;

fprintf('=== EDR fit diagnostic v2: robust (median) per-bin statistic ===\n');

fprintf('Loading surface + cortex mask...\n');
[vertices, ~] = read_vtk_surface(surface_file);
cortex = dlmread(cortex_mask_file); %#ok<DLMRD>
cortex_ind = find(cortex);
vertices_cortex = vertices(cortex_ind, :);

fprintf('Loading structural connectome...\n');
sc_data = load(sc_file);
SC = sc_data.avgSC_L;
clear sc_data;

fprintf('Computing per-bin median / trimmed-mean / mean of nonzero weights...\n');
fprintf('(this collects raw nonzero weights per bin -- expect higher memory use than the previous diagnostic)\n');
stats = compute_distance_bin_median(vertices_cortex, SC, n_bins_fit, fit_dist_range, block_size);
clear SC;

[A_mean, lambda_mean] = fit_edr_exponential(stats.centers, stats.mean_nonzero, exclude_bins);
[A_med, lambda_med] = fit_edr_exponential(stats.centers, stats.median_nonzero, exclude_bins);
[A_trim, lambda_trim] = fit_edr_exponential(stats.centers, stats.trimmed_mean_nonzero, exclude_bins);

fprintf('\n--- EDR fit comparison ---\n');
fprintf('  Paper (Vohryzek et al. 2025):  A = 0.066,  lambda = 0.162 mm^-1\n');
fprintf('  Nonzero MEAN fit:              A = %.4f, lambda = %.4f mm^-1\n', A_mean, lambda_mean);
fprintf('  Nonzero MEDIAN fit:            A = %.4f, lambda = %.4f mm^-1\n', A_med, lambda_med);
fprintf('  Nonzero TRIMMED-MEAN (95pct):  A = %.4f, lambda = %.4f mm^-1\n', A_trim, lambda_trim);

fprintf('\n--- Sample of per-bin values (first 30 fit-range bins) ---\n');
fit_bins = (exclude_bins+1):(exclude_bins+30);
for b = fit_bins
    fprintf('  d=%.1fmm: mean=%.4f median=%.4f trimmed=%.4f (n=%d)\n', ...
        stats.centers(b), stats.mean_nonzero(b), stats.median_nonzero(b), ...
        stats.trimmed_mean_nonzero(b), stats.count_nonzero(b));
end

out_dir = fullfile(data_dir, 'processed', 'eigenmodes_novel_v2');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
report.A_mean = A_mean; report.lambda_mean = lambda_mean;
report.A_med = A_med; report.lambda_med = lambda_med;
report.A_trim = A_trim; report.lambda_trim = lambda_trim;
save(fullfile(out_dir, 'edr_fit_diagnostic_median.mat'), 'report', 'stats');
fprintf('\nSaved -> %s\n', fullfile(out_dir, 'edr_fit_diagnostic_median.mat'));
fprintf('=== Diagnostic complete ===\n');
