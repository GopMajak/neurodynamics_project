%% Stage 1: Anatomical graph construction + eigenmode derivation
%
% Independent reimplementation of the EDR / EDR+LR graph construction and
% Laplacian eigenmode derivation described in:
%   Vohryzek J, Sanz-Perl Y, Kringelbach ML, Deco G. "Human brain dynamics
%   are shaped by rare long-range connections over and above cortical
%   geometry." PNAS 122(1):e2415102122 (2025). See Materials and Methods:
%   "EDR", "EDR binary", "EDR continuous", "EDR+LR", "Laplacian
%   decomposition".
%   Pang JC et al. "Geometric constraints on human brain function."
%   Nature 618:566-574 (2023). See Methods: "Derivation of connectome
%   eigenmodes", "Derivation of EDR eigenmodes".
%
% This script builds EDR binary, EDR continuous, and EDR+LR harmonic
% modes from the structural connectome and cortical surface. Geometric
% (LBO) modes are NOT rederived here -- they are taken as given input
% data (Pang et al.'s precomputed fsLR_32k_midthickness-lh_emode_200.txt),
% exactly as both papers treat the HCP connectome/task data as given
% input data. All graph-construction and eigendecomposition code here is
% original, independent of the Pang/Vohryzek toolboxes.
%
% Outputs -> data/processed/eigenmodes_novel/novel_<name>_eigenmodes_lh_200.mat
%   each containing: eig_vec (num_vertices x 200), eig_val (200 x 1)

clear; clc;

%% Paths
this_file = mfilename('fullpath');
code_dir = fileparts(fileparts(this_file)); % Code/
project_root = fileparts(code_dir);
addpath(fullfile(code_dir, 'utils'));

data_dir = fullfile(project_root, 'data');
out_dir = fullfile(data_dir, 'processed', 'eigenmodes_novel');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

surface_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_midthickness-lh.vtk');
cortex_mask_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt');
sc_file = fullfile(data_dir, 'empirical', 'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat');

%% Parameters (from papers' Methods)
num_modes = 200;
alpha_short = 0.12;      % EDR binary/continuous decay (Pang et al. 2023)
n_bins_fit = 400;        % bins for EDR parameter fit
fit_dist_range = [10, 170];
exclude_bins = 25;       % excluded from fit (shortest-distance bins)
n_bins_lr = 400;         % bins for LR exception detection
sd_thresh = 3;           % LR exception: > mean + 3*SD
min_lr_distance = 40;    % mm
weight_floor = 1e-6;     % sparsification cutoff for continuous weights
block_size = 2000;
rng_seed = 42;

fprintf('=== Stage 1: Anatomical graph construction ===\n');

%% Load surface + cortex mask
fprintf('Loading surface mesh...\n');
[vertices, ~] = read_vtk_surface(surface_file);
cortex = dlmread(cortex_mask_file); %#ok<DLMRD>
cortex_ind = find(cortex);
vertices_cortex = vertices(cortex_ind, :);
fprintf('  %d total vertices, %d cortical (non-medial-wall)\n', numel(cortex), numel(cortex_ind));

%% Load structural connectome
fprintf('Loading structural connectome (this is a large file)...\n');
sc_data = load(sc_file);
SC = sc_data.avgSC_L;
clear sc_data;
fprintf('  connectome size: %d x %d\n', size(SC, 1), size(SC, 2));

if size(SC, 1) ~= numel(cortex_ind)
    error('run_stage1:sizeMismatch', ...
        'Connectome size (%d) does not match cortex vertex count (%d).', ...
        size(SC, 1), numel(cortex_ind));
end

%% Pass 1a: bin SC vs distance for EDR parameter fit
fprintf('Binning connectome weights by distance (EDR parameter fit)...\n');
fit_stats = compute_distance_bin_stats(vertices_cortex, SC, n_bins_fit, fit_dist_range, block_size);
[A_fit, lambda_fit] = fit_edr_exponential(fit_stats.centers, fit_stats.mean, exclude_bins);
fprintf('  Fitted EDR params: A = %.4f, lambda = %.4f mm^-1\n', A_fit, lambda_fit);
fprintf('  (Vohryzek et al. 2025 report A = 0.066, lambda = 0.162 mm^-1)\n');

%% Pass 1b: bin SC vs distance for LR exception detection (full range)
fprintf('Binning connectome weights by distance (LR exception detection)...\n');
% Determine the observed distance range cheaply from a random vertex subsample
sub_idx = randperm(size(vertices_cortex, 1), min(2000, size(vertices_cortex, 1)));
sub_d = pdist(vertices_cortex(sub_idx, :));
lr_dist_range = [0, max(sub_d) * 1.02];
lr_stats = compute_distance_bin_stats(vertices_cortex, SC, n_bins_lr, lr_dist_range, block_size);

%% Pass 2: build EDR binary / continuous / EDR+LR adjacency matrices
fprintf('Building EDR binary, EDR continuous, and EDR+LR adjacency matrices...\n');
[A_edr_binary, A_edr_continuous, A_edr_lr, report] = build_anatomical_graphs_blocked( ...
    vertices_cortex, SC, alpha_short, lambda_fit, lr_stats, sd_thresh, ...
    min_lr_distance, weight_floor, block_size, rng_seed);

fprintf('  LR exceptions detected: %d\n', report.n_lr_exceptions);
fprintf('  Density -- binary: %.4f%%, continuous: %.4f%%, EDR+LR: %.4f%%\n', ...
    report.density_binary*100, report.density_continuous*100, report.density_edrlr*100);

clear SC; % free ~7 GB before eigendecomposition

%% Eigendecomposition
graphs = {A_edr_binary, A_edr_continuous, A_edr_lr};
names = {'EDRbinary', 'EDRcontinuous', 'EDRLR'};

for g = 1:numel(graphs)
    fprintf('Solving normalized Laplacian eigenproblem for %s...\n', names{g});
    tic;
    [eig_vec, eig_val] = normalized_laplacian_eigenmodes(graphs{g}, num_modes);
    fprintf('  done in %.1f s\n', toc);

    save_file = fullfile(out_dir, sprintf('novel_%s_eigenmodes_lh_200.mat', names{g}));
    save(save_file, 'eig_vec', 'eig_val', '-v7.3');
    fprintf('  saved -> %s\n', save_file);
end

%% Save diagnostic report
report.A_fit = A_fit;
report.lambda_fit = lambda_fit;
report.alpha_short = alpha_short;
report_file = fullfile(out_dir, 'stage1_report.mat');
save(report_file, 'report', 'fit_stats', 'lr_stats');
fprintf('Diagnostic report saved -> %s\n', report_file);
fprintf('=== Stage 1 complete ===\n');
