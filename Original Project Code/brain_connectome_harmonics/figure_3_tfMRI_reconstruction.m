% figure_3_tfMRI_reconstruction.m
%
% Reproduces Vohryzek et al. 2024 PNAS Figure 3 A/B/C: task-fMRI
% activation-map reconstruction accuracy (MSE) across 47 HCP tasks and
% 255 subjects, comparing the four connectome eigenmode bases.
%
% Adapted from the original PNAS_Figure_3_main_part_ABC.m. Recomputing MSE
% from raw activation maps for all 255 subjects x 47 tasks x 200 modes
% would be extremely slow, so this uses the authors' own precomputed
% per-subject MSE results instead (Results/long_3T_task/*.mat).
%
% Figure 3D (surface renderings of example reconstructions) is handled
% separately in figure_3D_surface_reconstruction.m, since it needs an
% actual reconstruction computed from eigenmodes rather than these
% precomputed MSE summaries.

%% Setup

edrlr_data_dir = fullfile('Code from Vohryzek et al 2024', 'vohryzek2024_EDRLR');
addpath(genpath(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', 'colormaps_add')));

condition_names = {'Geometry', 'EDR binary', 'EDR continuous', 'EDR+LR'};
representTask = [3, 11, 19, 30, 41, 44, 47]; % the 7 representative tasks used in the paper

%% Load task field names (used for labeling plots only)

data = load(fullfile(edrlr_data_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'), 'zstat');
fieldNames = fieldnames(data.zstat);
clear data

%% Load precomputed reconstruction MSE (47 tasks x 255 subjects x 200 modes)

all_files = {
    'tMRI_parcellated_reconstruction_all_subjects_Geometry.mat'
    'tMRI_parcellated_reconstruction_all_subjects_EDRbinary.mat'
    'tMRI_parcellated_reconstruction_all_subjects_EDRcontinuous.mat'
    'tMRI_parcellated_reconstruction_all_subjects_EDRLR.mat'
};
tk7_files = {
    'tMRI_parcellated_reconstruction_7tk_subjects_Geometry.mat'
    'tMRI_parcellated_reconstruction_7tk_subjects_EDRbinary.mat'
    'tMRI_parcellated_reconstruction_7tk_subjects_EDRcontinuous.mat'
    'tMRI_parcellated_reconstruction_7tk_subjects_EDRLR.mat'
};

recon_mse_parc_version = cell(1,4);      % {v} = [num_tk x num_modes], mean over subjects, all 47 tasks
recon_mse_parc_7tk_version = cell(1,4);  % {v} = [num_tk x num_modes], mean over subjects, only the 7 representative tasks
for v = 1:4
    s = load(fullfile(edrlr_data_dir, 'Results', 'long_3T_task', all_files{v}));
    recon_mse_parc_version{v} = squeeze(mean(s.recon_mse_parc_tk_sbj_version, 2, 'omitnan'))';

    s7 = load(fullfile(edrlr_data_dir, 'Results', 'long_3T_task', tk7_files{v}));
    recon_mse_parc_7tk_version{v} = squeeze(mean(s7.recon_mse_parc_tk_sbj_version, 2, 'omitnan'))';
end
num_modes = size(recon_mse_parc_version{1}, 1);
fprintf('Loaded task reconstruction MSE: %d modes x %d tasks x %d connectome types\n', ...
    num_modes, numel(fieldNames), numel(condition_names));

%% Figure 3A: normalized reconstruction MSE and its rate of change, 7 representative tasks

colors = {'m-', 'b-', 'g-', 'k-'};
nMod = 200;
figure('Name', 'Figure 3A - task reconstruction, 7 representative tasks');
for v = 1:4
    subplot(2,4,v)
    y = recon_mse_parc_7tk_version{v}(1:nMod, representTask);
    plot(1:nMod, y./max(y), colors{v}, 'linewidth', 2); grid on
    title(condition_names{v}); xlim([1,nMod]); ylim([0,1]); axis square

    subplot(2,4,v+4)
    y2 = recon_mse_parc_7tk_version{v}(2:nMod, representTask);
    plot(2:nMod, abs(diff([zeros(1,numel(representTask)); y2])), colors{v}, 'linewidth', 2); grid on
    title(condition_names{v}); xlim([4,nMod]); ylim([0,0.3]); axis square
end

%% Figure 3B: reconstruction MSE difference vs. Geometry, across all 47 tasks

th_mode = 20;
geom_ind = find(strcmp(condition_names, 'Geometry'));
diff_names = setdiff(1:4, geom_ind);

figure('Name', 'Figure 3B - task reconstruction difference vs Geometry');
for k = 1:numel(diff_names)
    v = diff_names(k);
    d = (recon_mse_parc_version{v}(1:th_mode,:) - recon_mse_parc_version{geom_ind}(1:th_mode,:))';

    subplot(2,3,k)
    imagesc(d, [-0.25, 0.25]); title(condition_names{v})
    xlabel('mode'); ylabel('task')

    subplot(2,3,k+3)
    bar(mean(d,1), 'k', 'linewidth', 2); title(condition_names{v}); ylim([-0.25, 0.25])
end
colormap(redblue)

%% Figure 3C: per-mode MSE contribution for one representative task, EDR+LR

edrlr_ind = find(strcmp(condition_names, 'EDR+LR'));
task_i = representTask(end); % last of the 7 representative tasks, as in the original script
figure('Name', 'Figure 3C - mode contribution to reconstruction (EDR+LR)');
bar(abs(diff([0; recon_mse_parc_version{edrlr_ind}(2:th_mode,task_i)])), 'k', 'linewidth', 2)
grid on; ylim([0, 0.05]); xlim([0.5, th_mode-0.5])
title(sprintf('Mode contribution to reconstruction for task %s', fieldNames{task_i}), 'Interpreter', 'none')
ylabel('FC MSE contribution'); xlabel('Modes')
