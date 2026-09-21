%% Diagnostic v3: faithful replication of Vohryzek's exact fitting recipe
%
% Prior diagnostics (diagnose_edr_fit.m, diagnose_edr_fit_median.m) ruled
% out zero-inflation and confirmed the underlying statistic Vohryzek's own
% longrange_derivation_project_laplacian_v2.m uses IS the plain all-pairs
% mean (matching Stage 1's original choice) -- so the bug isn't the
% statistic. Reading that script line-by-line surfaced three concrete
% procedural differences to test instead:
%   1. Connectome normalized by its global max BEFORE binning (C/max(C(:)))
%   2. Bins span the full [0, max(distance)] range in 400 equal bins, not
%      a hardcoded [10,170]mm range
%   3. lsqcurvefit uses a FIXED initial guess A0=[0.15, 0.18], not one
%      derived from a log-linear regression -- plausible source of a bad
%      local-minimum convergence on noisy real data
%
% NEW FILE -- does not modify run_stage1.m or any prior diagnostic.

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
project_root = fileparts(code_dir);
addpath(fullfile(code_dir, 'utils'));

data_dir = fullfile(project_root, 'data');
surface_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_midthickness-lh.vtk');
cortex_mask_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt');
sc_file = fullfile(data_dir, 'empirical', 'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat');

n_bins = 400;
fit_start_bin = 25;
block_size = 2000;

fprintf('=== EDR fit diagnostic v3: faithful replication ===\n');

fprintf('Loading surface + cortex mask...\n');
[vertices, ~] = read_vtk_surface(surface_file);
cortex = dlmread(cortex_mask_file); %#ok<DLMRD>
cortex_ind = find(cortex);
vertices_cortex = vertices(cortex_ind, :);

fprintf('Loading structural connectome...\n');
sc_data = load(sc_file);
SC = sc_data.avgSC_L;
clear sc_data;

%% Step 1: normalize by global max
max_sc = max(SC(:));
fprintf('Global max connectome weight: %.6f\n', max_sc);
SC_norm = SC / max_sc;
clear SC;

%% Step 2: find true max distance, bin over [0, max_distance]
fprintf('Finding true max pairwise distance (blocked)...\n');
n = size(vertices_cortex, 1);
sq_norms = sum(vertices_cortex.^2, 2);
max_dist = 0;
for start_row = 1:block_size:n
    end_row = min(start_row + block_size - 1, n);
    rows = (start_row:end_row)';
    bv = vertices_cortex(rows, :);
    d2 = sq_norms(rows) - 2*(bv * vertices_cortex') + sq_norms';
    d2(d2 < 0) = 0;
    max_dist = max(max_dist, sqrt(max(d2(:))));
end
fprintf('  max distance = %.4f mm\n', max_dist);

fprintf('Binning (all-pairs mean, matching Vohryzek''s own script) over [0, %.2f] mm, %d bins...\n', max_dist, n_bins);
stats = compute_distance_bin_stats(vertices_cortex, SC_norm, n_bins, [0, max_dist], block_size);

%% Step 3: fit with their exact initial guess/bounds
[A_fixed, lambda_fixed] = fit_edr_exponential_fixed_init(stats.centers, stats.mean, fit_start_bin, [0.15, 0.18], [-100, 100]);

%% For comparison: same data, our original log-linear-init fit
[A_loglinear, lambda_loglinear] = fit_edr_exponential(stats.centers, stats.mean, fit_start_bin - 1);

fprintf('\n--- EDR fit comparison ---\n');
fprintf('  Paper (Vohryzek et al. 2025):              A = 0.066,  lambda = 0.162 mm^-1\n');
fprintf('  Faithful replication (fixed A0=[0.15,0.18]): A = %.4f, lambda = %.4f mm^-1\n', A_fixed, lambda_fixed);
fprintf('  Same data, log-linear-init fit instead:      A = %.4f, lambda = %.4f mm^-1\n', A_loglinear, lambda_loglinear);

out_dir = fullfile(data_dir, 'processed', 'eigenmodes_novel_v2');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
report.max_sc = max_sc;
report.max_dist = max_dist;
report.A_fixed = A_fixed; report.lambda_fixed = lambda_fixed;
report.A_loglinear = A_loglinear; report.lambda_loglinear = lambda_loglinear;
save(fullfile(out_dir, 'edr_fit_diagnostic_faithful.mat'), 'report', 'stats');
fprintf('\nSaved -> %s\n', fullfile(out_dir, 'edr_fit_diagnostic_faithful.mat'));
fprintf('=== Diagnostic complete ===\n');
