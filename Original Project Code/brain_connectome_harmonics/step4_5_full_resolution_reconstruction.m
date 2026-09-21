% step4_5_full_resolution_reconstruction.m
%
% Full replication of Vohryzek et al. (2025)'s core result at the real
% 29,696-vertex resolution (step5_static_reconstruction.m does the same
% thing but at the faster, lower-resolution Glasser360 pilot scale):
% reconstruction accuracy vs. number of modes, for all four connectome
% types (Geometry, EDR binary, EDR continuous, EDR+LR), across all 47 HCP
% task contrasts.
%
% Reuses the checkpointed EDR+LR connectome and EDR fit saved by Stage 2
% of connectome_harmonics.m, instead of re-deriving it from scratch (that
% derivation is slow, ~170s, and doesn't need to be repeated once saved).
%
% EDR-continuous and EDR-binary are rebuilt using Pang et al.'s fixed
% decay rate (lambda = 0.12/mm, per Vohryzek et al. 2025's Methods) - NOT
% the checkpoint's own Afit/lambda, which is fit fresh to this dataset and
% is only used for the EDR+LR baseline. Using the same lambda for both
% would make EDR-continuous and EDR+LR share the same baseline, which
% would hide the very difference this analysis is trying to measure.
%
% Eigenmodes are computed at full vertex resolution, but reconstruction
% error is measured after parcellating both the reconstruction and the
% real data down to Glasser360 (180 parcels), matching the paper's Methods
% exactly. An earlier version of this script compared raw vertex-level
% maps instead, but that comparison is dominated by generic large-scale
% smoothness that's similar across any basis, and can hide a real but
% more modest effect like the long-range splice.

%% Setup

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Code from Vohryzek et al 2024', 'vohryzek2024_EDRLR');

addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

if ~exist('results', 'dir'); mkdir('results'); end

hemisphere = 'lh';
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';
num_modes = 200; % matches Vohryzek/Pang's published value

%% Load surface geometry and cortex mask (needed for the vertex distance matrix)

[vertices, faces] = read_vtk(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
num_vertices = length(cortex);
vertices_cortex = surface_midthickness.vertices(cortex_ind, :);

fprintf('Loaded surface: %d total vertices, %d cortex vertices\n', num_vertices, numel(cortex_ind));

% Glasser360 parcellation, used only to evaluate reconstruction error
% below (not for computing the eigenmodes themselves, which stay at full
% vertex resolution) - matches the paper's Methods, which measure
% reconstruction error at this parcel resolution rather than per-vertex.
parc_name = 'Glasser360';
parc = dlmread(fullfile(edrlr_data_dir, 'Data', 'parcellations', ...
    sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
fprintf('Loaded %s parcellation for evaluation: %d parcels\n', ...
    parc_name, numel(unique(parc_cortex(parc_cortex>0))));

%% Load the checkpointed EDR+LR connectome (from Stage 2 of connectome_harmonics.m)

checkpoint_file = fullfile('results', 'EDR_LR_connectome_full_resolution.mat');
if ~exist(checkpoint_file, 'file')
    error(['Checkpoint file not found: %s\nRun Stage 2 of connectome_harmonics.m ' ...
        'first (the connectome-derivation part only -- the eigendecomposition ' ...
        'that used to crash there is superseded by this script).'], checkpoint_file);
end
fprintf('Loading checkpointed EDR+LR connectome and EDR fit...\n');
load(checkpoint_file, 'EDR_LRE', 'Afit', 'lambda');
fprintf('Loaded EDR fit: A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

%% Rebuild EDR-continuous and EDR-binary using Pang et al.'s fixed decay rate
%
% Uses lambda_pang = 0.12, not the checkpoint's own Afit/lambda (which is
% fit fresh to this dataset and used only for the EDR+LR baseline above).
% This matches Vohryzek's own construction of these two graphs: a
% synthetic distance-probability surface, not a fit to the empirical
% connectome. Reusing Afit/lambda here (as an earlier version of this
% script did) would make EDR-continuous share EDR+LR's own baseline,
% which defeats the comparison this whole analysis is trying to make.
lambda_pang = 0.12;

tic
rr = squareform(pdist(single(vertices_cortex)));
fprintf('Computed %dx%d vertex distance matrix in %.1f s\n', size(rr,1), size(rr,2), toc);

Pspace = single(exp(-lambda_pang * double(rr)));
N = size(Pspace, 1);
Pspace(1:(N+1):end) = 0;
Pspace = Pspace / max(Pspace(:));

EDR_conn_full = Pspace;

rng(1); % matches Vohryzek's own EDR-binary construction seed
rand_prob = rand(size(Pspace), 'single');
rand_prob = triu(rand_prob,1) + triu(rand_prob,1)';
EDR_binary_full = single(rand_prob < Pspace);
clear rand_prob Pspace rr

EDR_conn_full(1:(N+1):end) = 0;
EDR_binary_full(1:(N+1):end) = 0;

%% Compute eigenmodes for the three network-derived graphs (Geometry is precomputed)

fprintf('\nEigendecomposing EDR continuous (%dx%d, %d modes)...\n', N, N, num_modes);
tic; [eig_vec_EDR, eig_val_EDR] = calc_network_eigenmode_lowmem(EDR_conn_full, num_modes); toc
clear EDR_conn_full

fprintf('\nEigendecomposing EDR binary (%dx%d, %d modes)...\n', N, N, num_modes);
tic; [eig_vec_EDRbin, eig_val_EDRbin] = calc_network_eigenmode_lowmem(EDR_binary_full, num_modes); toc
clear EDR_binary_full

fprintf('\nEigendecomposing EDR+LR (%dx%d, %d modes)...\n', N, N, num_modes);
tic; [eig_vec_EDRLR, eig_val_EDRLR] = calc_network_eigenmode_lowmem(EDR_LRE, num_modes); toc
clear EDR_LRE

fprintf('\nLoading precomputed Geometry eigenmodes...\n');
geom_evec_full = dlmread(fullfile(eigenmode_toolbox_dir, 'data', 'results', ...
    sprintf('basis_geometric_%s-%s_evec_200.txt', mesh_interest, hemisphere)));
eig_vec_geom = single(geom_evec_full(cortex_ind, 1:num_modes));
clear geom_evec_full

graph_names = {'Geometry', 'EDR binary', 'EDR continuous', 'EDR+LR'};
graph_colors = {'m', 'b', 'g', 'k'};
graph_eigvecs = {eig_vec_geom, eig_vec_EDRbin, eig_vec_EDR, eig_vec_EDRLR};
num_graphs = numel(graph_names);

%% Load task fMRI data and average across subjects, restricted to cortex vertices

fprintf('\nLoading task fMRI data (47 contrasts, 255-subject group average)...\n');
data = load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
task_names = fieldnames(data.zstat);
num_tasks = numel(task_names);

task_maps_cortex = zeros(numel(cortex_ind), num_tasks, 'single');
for t = 1:num_tasks
    group_avg = nanmean(data.zstat.(task_names{t}), 2); % [32492 x 1]
    task_maps_cortex(:, t) = single(group_avg(cortex_ind));
end
clear data

%% Reconstruction accuracy vs. number of modes, all four graphs
%
% The projection (beta) is computed using the full vertex-resolution data
% and eigenmodes. Reconstruction error is then measured after
% parcellating both the reconstruction and the real data to Glasser360,
% matching the paper's Methods. Parcellating is just an averaging
% operation, so it "commutes" with the reconstruction: parcellating the
% eigenmode basis once per graph (cheap) gives the same result as
% reconstructing and then parcellating for every mode count, but is much
% faster since it avoids redoing that step 47*4*200 times.

fprintf('\nParcellating eigenmode bases to %s for reconstruction-error evaluation...\n', parc_name);
graph_eigvecs_parc = cellfun(@(v) calc_parcellate(parc_cortex, double(v)), graph_eigvecs, 'UniformOutput', false);

fprintf('\nRunning reconstruction (accuracy vs. N modes, %d tasks x %d graphs, %d modes)...\n', ...
    num_tasks, num_graphs, num_modes);

corr_all = cell(1, num_graphs);
mse_all  = cell(1, num_graphs);
for g = 1:num_graphs
    corr_all{g} = nan(num_tasks, num_modes);
    mse_all{g}  = nan(num_tasks, num_modes);
end

tic
for t = 1:num_tasks
    y = double(task_maps_cortex(:, t));
    y_parc = calc_parcellate(parc_cortex, y);

    for g = 1:num_graphs
        basis = double(graph_eigvecs{g});
        basis_parc = graph_eigvecs_parc{g};
        beta = calc_eigendecomposition(y, basis, 'matrix');
        for N_modes = 1:num_modes
            y_hat_parc = basis_parc(:,1:N_modes) * beta(1:N_modes);
            corr_all{g}(t,N_modes) = corr(y_parc, y_hat_parc);
            mse_all{g}(t,N_modes)  = immse(y_parc, y_hat_parc);
        end
    end
    if mod(t,10) == 0
        fprintf('  ...task %d/%d done (%.1f s elapsed)\n', t, num_tasks, toc);
    end
end
fprintf('Done in %.1f s.\n', toc);

mean_corr = cellfun(@(x) mean(x,1), corr_all, 'UniformOutput', false);
mean_mse  = cellfun(@(x) mean(x,1), mse_all, 'UniformOutput', false);

%% Statistics: paired comparisons at N=20 (Vohryzek's reported peak-advantage range)

N_test = min(20, num_modes);
idx_EDR   = find(strcmp(graph_names, 'EDR continuous'));
idx_EDRLR = find(strcmp(graph_names, 'EDR+LR'));
idx_geom  = find(strcmp(graph_names, 'Geometry'));

fprintf('\n=== Paired comparisons at N=%d modes, across %d tasks ===\n', N_test, num_tasks);
for g = 1:num_graphs
    fprintf('%-16s mean corr = %.4f, mean MSE = %.4f\n', graph_names{g}, ...
        mean_corr{g}(N_test), mean_mse{g}(N_test));
end

[~, p_corr_edrlr_edr, ~, s_corr_edrlr_edr] = ttest(corr_all{idx_EDRLR}(:,N_test), corr_all{idx_EDR}(:,N_test));
[~, p_mse_edrlr_edr,  ~, s_mse_edrlr_edr]  = ttest(mse_all{idx_EDRLR}(:,N_test),  mse_all{idx_EDR}(:,N_test));
fprintf('\nPrimary comparison, EDR+LR vs. EDR-only:\n');
fprintf('  Correlation: paired t(%d) = %.3f, p = %.4g\n', s_corr_edrlr_edr.df, s_corr_edrlr_edr.tstat, p_corr_edrlr_edr);
fprintf('  MSE:         paired t(%d) = %.3f, p = %.4g\n', s_mse_edrlr_edr.df, s_mse_edrlr_edr.tstat, p_mse_edrlr_edr);

[~, p_corr_edrlr_geom, ~, s_corr_edrlr_geom] = ttest(corr_all{idx_EDRLR}(:,N_test), corr_all{idx_geom}(:,N_test));
[~, p_mse_edrlr_geom,  ~, s_mse_edrlr_geom]  = ttest(mse_all{idx_EDRLR}(:,N_test),  mse_all{idx_geom}(:,N_test));
fprintf('\nReplication check, EDR+LR vs. Geometry (Vohryzek''s own headline comparison):\n');
fprintf('  Correlation: paired t(%d) = %.3f, p = %.4g\n', s_corr_edrlr_geom.df, s_corr_edrlr_geom.tstat, p_corr_edrlr_geom);
fprintf('  MSE:         paired t(%d) = %.3f, p = %.4g\n', s_mse_edrlr_geom.df, s_mse_edrlr_geom.tstat, p_mse_edrlr_geom);

%% Save results

save(fullfile('results', 'step4_5_full_resolution_reconstruction.mat'), ...
    'task_names', 'graph_names', 'corr_all', 'mse_all', ...
    'eig_vec_geom', 'eig_vec_EDRbin', 'eig_val_EDRbin', ...
    'eig_vec_EDR', 'eig_val_EDR', 'eig_vec_EDRLR', 'eig_val_EDRLR', ...
    'Afit', 'lambda', '-v7.3');
fprintf('\nSaved results to results/step4_5_full_resolution_reconstruction.mat\n');

%% Figure

figure('Name', 'Full-resolution static reconstruction, all four graphs');
subplot(1,2,1);
hold on
for g = 1:num_graphs
    plot(1:num_modes, mean_corr{g}, [graph_colors{g} '-'], 'LineWidth', 2);
end
xlabel('Number of modes (N)'); ylabel('Correlation with empirical activity');
legend(graph_names, 'Location', 'southeast'); grid on
title('Reconstruction accuracy (correlation)'); xlim([1 num_modes]);

subplot(1,2,2);
hold on
for g = 1:num_graphs
    plot(1:num_modes, mean_mse{g}, [graph_colors{g} '-'], 'LineWidth', 2);
end
xlabel('Number of modes (N)'); ylabel('MSE vs. empirical activity');
legend(graph_names, 'Location', 'northeast'); grid on
title('Reconstruction error (MSE)'); xlim([1 num_modes]);

sgtitle(sprintf('Full vertex resolution (%d cortex vertices), %d-task average', numel(cortex_ind), num_tasks));
saveas(gcf, fullfile('results', 'step4_5_full_resolution_reconstruction.png'));
fprintf('Saved figure to results/step4_5_full_resolution_reconstruction.png\n');
