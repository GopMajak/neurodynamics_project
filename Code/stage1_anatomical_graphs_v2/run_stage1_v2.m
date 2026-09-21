%% Stage 1 v2: corrected anatomical graph construction + eigenmode derivation
%
% Corrected version of Code/stage1_anatomical_graphs/run_stage1.m, using
% the faithful EDR fitting recipe validated in diagnose_edr_fit_faithful.m
% (A=0.0658, lambda=0.1616 -- matches Vohryzek et al. 2025's reported
% 0.066/0.162) and the exact LR-exception procedure read from their own
% longrange_derivation_project_laplacian_v2.m (see
% build_anatomical_graphs_blocked_v2.m for the full rationale).
%
% NEW FILE -- does not modify or overwrite run_stage1.m or any of its
% outputs (data/processed/eigenmodes_novel/). This script writes to a
% separate output directory (data/processed/eigenmodes_novel_v2/) so the
% original (bugged-fit) and corrected results exist side by side for
% comparison, alongside Vohryzek's own published/reproduced results.

clear; clc;

%% Paths
this_file = mfilename('fullpath');
code_dir = fileparts(fileparts(this_file));
project_root = fileparts(code_dir);
addpath(fullfile(code_dir, 'utils'));

data_dir = fullfile(project_root, 'data');
out_dir = fullfile(data_dir, 'processed', 'eigenmodes_novel_v2');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

surface_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_midthickness-lh.vtk');
cortex_mask_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt');
sc_file = fullfile(data_dir, 'empirical', 'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat');

%% Parameters
num_modes = 200;
alpha_short = 0.12;      % EDR binary/continuous decay (Pang et al. 2023, fixed constant)
n_bins = 400;
fit_start_bin = 25;      % EDR parameter fit: bins 25:400 (paper: "excluding first 25 bins")
nr_ini = 20;             % LR exception detection: bins 20:380 (Vohryzek's own script)
nr_fin = 380;
nstd = 3;
min_lr_distance = 40;    % mm
weight_floor = 1e-6;
block_size = 2000;
rng_seed = 42;

fprintf('=== Stage 1 v2: corrected anatomical graph construction ===\n');

%% Load surface + cortex mask
fprintf('Loading surface mesh...\n');
[vertices, ~] = read_vtk_surface(surface_file);
cortex = dlmread(cortex_mask_file); %#ok<DLMRD>
cortex_ind = find(cortex);
vertices_cortex = vertices(cortex_ind, :);
fprintf('  %d cortical vertices\n', numel(cortex_ind));

%% Load + normalize structural connectome
fprintf('Loading structural connectome...\n');
sc_data = load(sc_file);
SC = sc_data.avgSC_L;
clear sc_data;
max_sc = max(SC(:));
C_norm = SC / max_sc;
clear SC;
fprintf('  normalized by global max = %.6f\n', max_sc);

%% True max distance
fprintf('Finding true max pairwise distance...\n');
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

%% Distance-binned stats (400 bins over [0, max_dist]) on normalized connectome
fprintf('Binning connectome weights by distance...\n');
stats = compute_distance_bin_stats(vertices_cortex, C_norm, n_bins, [0, max_dist], block_size);

%% Fit EDR exponential (faithful recipe)
[A_fit, lambda_fit] = fit_edr_exponential_fixed_init(stats.centers, stats.mean, fit_start_bin, [0.15, 0.18], [-100, 100]);
fprintf('Fitted EDR params: A = %.4f, lambda = %.4f mm^-1 (paper: A=0.066, lambda=0.162)\n', A_fit, lambda_fit);

%% Build graphs
fprintf('Building EDR binary, EDR continuous, and EDR+LR adjacency matrices...\n');
[A_edr_binary, A_edr_continuous, A_edr_lr, report] = build_anatomical_graphs_blocked_v2( ...
    vertices_cortex, C_norm, alpha_short, A_fit, lambda_fit, stats.edges, stats.mean, stats.std, ...
    nr_ini, nr_fin, nstd, min_lr_distance, weight_floor, block_size, rng_seed);

fprintf('  LR exceptions detected: %d (%.4f%% of all pairs)\n', report.n_lr_exceptions, report.pct_lr_exceptions);
fprintf('  Density -- binary: %.4f%%, continuous: %.4f%%, EDR+LR: %.4f%%\n', ...
    report.density_binary*100, report.density_continuous*100, report.density_edrlr*100);

clear C_norm;

%% Eigendecomposition
graphs = {A_edr_binary, A_edr_continuous, A_edr_lr};
names = {'EDRbinary', 'EDRcontinuous', 'EDRLR'};

for g = 1:numel(graphs)
    fprintf('Solving normalized Laplacian eigenproblem for %s...\n', names{g});
    tic;
    [eig_vec, eig_val] = normalized_laplacian_eigenmodes(graphs{g}, num_modes);
    fprintf('  done in %.1f s\n', toc);

    save_file = fullfile(out_dir, sprintf('novel_v2_%s_eigenmodes_lh_200.mat', names{g}));
    save(save_file, 'eig_vec', 'eig_val', '-v7.3');
    fprintf('  saved -> %s\n', save_file);
end

%% Save diagnostic report
report.A_fit = A_fit;
report.lambda_fit = lambda_fit;
report.max_sc = max_sc;
report.max_dist = max_dist;
save(fullfile(out_dir, 'stage1_v2_report.mat'), 'report', 'stats');
fprintf('Report saved -> %s\n', fullfile(out_dir, 'stage1_v2_report.mat'));
fprintf('=== Stage 1 v2 complete ===\n');
