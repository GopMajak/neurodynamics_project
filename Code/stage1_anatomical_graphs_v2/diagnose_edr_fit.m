%% Diagnostic: why did the Stage 1 EDR fit come out wrong?
%
% Stage 1 (Code/stage1_anatomical_graphs/run_stage1.m) fit the EDR
% exponential using the mean structural-connectome weight per distance
% bin, INCLUDING zero-weight pairs (most vertex pairs in a tractography-
% derived connectome have zero detected streamlines). This script checks
% whether that zero-inflation -- if it varies systematically with
% distance -- distorted the fitted decay rate (lambda came out ~0.0185
% instead of the paper's reported 0.162, and flagged an implausibly high
% 1.4% of all pairs as "long-range exceptions").
%
% This is a NEW, standalone diagnostic. It does not modify or overwrite
% run_stage1.m, build_anatomical_graphs_blocked.m, or any other Stage 1
% file -- those results are preserved for side-by-side comparison against
% whatever corrected pipeline follows from this diagnosis.

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

fprintf('=== EDR fit diagnostic ===\n');

fprintf('Loading surface + cortex mask...\n');
[vertices, ~] = read_vtk_surface(surface_file);
cortex = dlmread(cortex_mask_file); %#ok<DLMRD>
cortex_ind = find(cortex);
vertices_cortex = vertices(cortex_ind, :);

fprintf('Loading structural connectome...\n');
sc_data = load(sc_file);
SC = sc_data.avgSC_L;
clear sc_data;

%% Global sparsity / distribution check (cheap, no distance binning needed)
fprintf('\n--- Global connectome weight distribution ---\n');
sc_vec = SC(:);
n_total = numel(sc_vec);
n_zero = sum(sc_vec == 0);
fprintf('  Total entries: %d\n', n_total);
fprintf('  Zero entries:  %d (%.2f%%)\n', n_zero, 100*n_zero/n_total);
nz = sc_vec(sc_vec > 0);
fprintf('  Nonzero entries: %d (%.2f%%)\n', numel(nz), 100*numel(nz)/n_total);
fprintf('  Nonzero weight percentiles: p50=%.4g  p90=%.4g  p99=%.4g  p99.9=%.4g  max=%.4g\n', ...
    prctile(nz, 50), prctile(nz, 90), prctile(nz, 99), prctile(nz, 99.9), max(nz));
clear sc_vec nz;

%% Distance-binned stats: all-pairs (zeros included) vs. nonzero-only
fprintf('\nComputing distance-binned statistics (all-pairs vs. nonzero-only)...\n');
stats = compute_distance_bin_stats_v2(vertices_cortex, SC, n_bins_fit, fit_dist_range, block_size);

fprintf('\n--- Fraction of pairs with nonzero weight, by distance ---\n');
sample_bins = round(linspace(1, n_bins_fit, 10));
for b = sample_bins
    fprintf('  distance ~%.1f mm: frac_nonzero = %.4f (n_nonzero=%d / n_total=%d)\n', ...
        stats.centers(b), stats.frac_nonzero(b), stats.count_nonzero(b), stats.count_all(b));
end

%% Refit EDR exponential both ways
[A_all, lambda_all] = fit_edr_exponential(stats.centers, stats.mean_all, exclude_bins);
[A_nz, lambda_nz] = fit_edr_exponential(stats.centers, stats.mean_nonzero, exclude_bins);

fprintf('\n--- EDR fit comparison ---\n');
fprintf('  Paper (Vohryzek et al. 2025):        A = 0.066,  lambda = 0.162 mm^-1\n');
fprintf('  Stage 1 original (all pairs, zeros included): A = %.4f, lambda = %.4f mm^-1\n', A_all, lambda_all);
fprintf('  This diagnostic (nonzero pairs only):          A = %.4f, lambda = %.4f mm^-1\n', A_nz, lambda_nz);

%% Save
report.n_zero_pct = 100*n_zero/n_total;
report.A_all = A_all; report.lambda_all = lambda_all;
report.A_nz = A_nz; report.lambda_nz = lambda_nz;
out_dir = fullfile(data_dir, 'processed', 'eigenmodes_novel_v2');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
save(fullfile(out_dir, 'edr_fit_diagnostic.mat'), 'report', 'stats');
fprintf('\nSaved -> %s\n', fullfile(out_dir, 'edr_fit_diagnostic.mat'));
fprintf('=== Diagnostic complete ===\n');
