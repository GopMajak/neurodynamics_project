% figure_2_fMRI_reconstruction.m
%
% Reproduces Vohryzek et al. 2024 PNAS Figure 2: how well resting-state
% fMRI functional connectivity can be reconstructed from four different
% connectome eigenmode bases (Geometry, EDR binary, EDR continuous, EDR+LR),
% measured by mean squared error (MSE).
%
% Adapted from the original PNAS_Figure_2_main.m. The original script can
% recompute MSE directly from 255 subjects' raw fMRI data, but that raw
% data isn't included in the code release, so this loads the authors' own
% precomputed per-subject MSE results instead (Results/long_3T/*.mat).
%
% Other changes from the original: file paths fixed to this project's
% layout, and the daviolinplot function (not included in the release)
% replaced with MATLAB's built-in boxchart.

%% Setup

edrlr_data_dir = fullfile('Code from Vohryzek et al 2024', 'vohryzek2024_EDRLR');

%% Load precomputed per-subject reconstruction MSE (255 subjects x 200 modes)

condition_names = {'Geometry', 'EDR binary', 'EDR continuous', 'EDR+LR'};
result_files = {
    'fMRI_parc_reconstruction_all_subjects_Geometry.mat'
    'fMRI_parc_reconstruction_all_subjects_EDRbinary.mat'
    'fMRI_parc_reconstruction_all_subjects_EDRcontinuous.mat'
    'fMRI_parc_reconstruction_all_subjects_EDRLR.mat'
};

recon_mse_parc_version_all = [];
recon_mse_parc_SC_exceptions_version_all = [];
for v = 1:numel(result_files)
    s = load(fullfile(edrlr_data_dir, 'Results', 'long_3T', result_files{v}));
    recon_mse_parc_version_all(v,:,:) = s.recon_mse_parc_version; %#ok<SAGROW>
    recon_mse_parc_SC_exceptions_version_all(v,:,:,:) = s.recon_mse_parc_SC_exceptions_version; %#ok<SAGROW>
end
num_sbj = size(recon_mse_parc_version_all, 2);
num_modes = size(recon_mse_parc_version_all, 3);
fprintf('Loaded reconstruction MSE for %d subjects x %d modes x %d connectome types\n', ...
    num_sbj, num_modes, numel(result_files));

%% Figure 2: MSE at 200 modes, restricted to long-range functional connections (>40mm)

data2mse = squeeze(recon_mse_parc_SC_exceptions_version_all(:,:,1,200))'; % [num_sbj x 4]
group_inx = repelem(1:numel(condition_names), num_sbj);

figure('Name', 'Figure 2 - MSE at 200 modes (functional long-range exceptions)');
boxchart(categorical(group_inx, 1:numel(condition_names), condition_names), data2mse(:));
hold on
for v = 1:numel(condition_names)
    xj = v + (rand(num_sbj,1)-0.5)*0.15;
    plot(xj, data2mse(:,v), 'k.', 'markersize', 4)
end
ylabel('MSE'); grid on; set(gca, 'fontsize', 10)
title('Reconstruction MSE, long-range functional connections (>40mm), 200 modes')

%% Figure 2: MSE vs. number of modes, same long-range restriction

figure('Name', 'Figure 2 - MSE vs number of modes');
plot(squeeze(mean(recon_mse_parc_SC_exceptions_version_all(:,:,1,:), 2))', 'linewidth', 2)
axis square; grid on; ylim([0, 0.05])
legend(condition_names, 'Location', 'northeast')
ylabel('MSE'); xlabel('Modes');
set(gcf, 'Color', [1 1 1]);
title('Reconstruction MSE vs number of modes (mean across 255 subjects)')

%% Statistics: paired t-tests at 200 modes, each condition vs. EDR+LR

stats_modes_th = squeeze(recon_mse_parc_SC_exceptions_version_all(:,:,1,200)); % [4 x num_sbj]
edrlr_ind = find(strcmp(condition_names, 'EDR+LR'));
fprintf('\nPaired t-tests vs EDR+LR at 200 modes (Bonferroni-corrected p, n=3 comparisons):\n');
for v = 1:numel(condition_names)
    if v == edrlr_ind, continue; end
    [~, p] = ttest(stats_modes_th(v,:), stats_modes_th(edrlr_ind,:));
    fprintf('  %s vs EDR+LR: p = %.4g (corrected p = %.4g)\n', condition_names{v}, p, min(p*3,1));
end
