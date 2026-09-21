%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3 ABC
%%% Original: James Pang, Monash University,
%%% updated: Jakub Vohryzek, Universitat Pompeu Fabra
%%%
%%% PATCHED for Neurodynamics_Project (Goal A reproduction):
%%%  - original hardcoded repo_dir/addpath lines commented out, replaced by
%%%    a path-setup block (see PNAS_Figure_2_main.m for the same pattern).
%%%  No other changes -- this script has no third-party plotting
%%%  dependencies and its default (parfor_run = 0) branch only loads
%%%  precomputed Results/*.mat files we already have.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- PATCH: path setup ---
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir));
pang_dir = fullfile(project_root, 'Reference Code Material', 'pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
vohryzek_dir = fullfile(project_root, 'Reference Code Material', 'Vohryzeketal2024', 'vohryzek2024_EDRLR');
addpath(genpath(pang_dir));
addpath(genpath(vohryzek_dir));
cd(vohryzek_dir);
fig_out_dir = fullfile(project_root, 'Figures', 'reproduction');
if ~exist(fig_out_dir, 'dir'); mkdir(fig_out_dir); end
%% --- END PATCH ---

%% loading data
% repo_dir = '/BrainEigenmodes_EDRLR_PNAS-main';  % ORIGINAL -- superseded by PATCH above
% addpath(genpath(repo_dir))                       % ORIGINAL -- superseded by PATCH above

%% Load surface files for visualization
hemisphere = 'lh';
num_modes = 200;
num_sbj = 255;
num_tk = 47;
MRI_type = 'tMRI' % 'fMRI' or 'tMRI'
typeRecon = 'original' % 'parcellated' or 'original'
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(sprintf('Data/template_surfaces/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';
% Load cortex mask
cortex = dlmread(sprintf('Data/template_surfaces/%s_cortex-%s_mask.txt', surface_interest, hemisphere));
cortex_ind = find(cortex);

% Calculate number of vertices
num_vertices = length(cortex);
%% parcellated template
% At parcellated level
parc_name = 'Glasser360';
parc = dlmread(sprintf('Data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere));
num_parcels = length(unique(parc(parc>0)));
%%
data = load('Data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat');
fieldNames = fieldnames(data.zstat);
%% subject list

% Read the subject IDs from the text file
subjectList = importdata('Data/empirical/subject_list_HCP.txt');

%%
parfor_run = 0;
if parfor_run == 1
parfor i = [1,2,3,4]
    data_connectome_file = '';

    %% extract connectomes
    connectome_version = i;
    if connectome_version == 1
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRbinary_eigenmodes-%s-%s_%i.mat',surface_interest, hemisphere, num_modes);
    elseif connectome_version == 2
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRcontinuous_eigenmodes-%s-%s_%i.mat',surface_interest, hemisphere, num_modes);
    elseif connectome_version == 3
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_Geometry_eigenmodes_%s-%s_%i.mat',surface_interest, hemisphere, num_modes);
    elseif connectome_version == 4
        data_connectome_file = sprintf('Connectome_derivation/Connectomes/synthetic_EDRLR_version2_eigenmodes_%s-%s_%i.mat',surface_interest, hemisphere, num_modes);
    end
    data_connectome = load(data_connectome_file);
    eigenmodes = data_connectome.eig_vec;
    %% pre-allocation

    recon_mse_parc_tk_sbj_version1 = zeros(num_tk,num_sbj, num_modes);
    recon_mse_parc_tk_sbj_version2 = zeros(num_tk,num_sbj, num_modes);
    recon_mse_parc_tk_sbj_version3 = zeros(num_tk,num_sbj, num_modes);
    recon_mse_parc_tk_sbj_version4 = zeros(num_tk,num_sbj, num_modes);
    % Load example single-subject tfMRI z-stat data

    data = load('/project_laplacian/BrainEigenmodes-main/data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat');
    fieldNames = fieldnames(data.zstat);
    for tk = [3, 11, 19, 30, 41, 44, 47] % the 7 representative tasks

        for sbj = 1:255
            fieldName = fieldNames{tk};
            disp(['connectome v',num2str(i),' task ',fieldName,' subject #',num2str(sbj)])
            data_to_reconstruct = data.zstat.(fieldName)(:,sbj);

            recon_beta = zeros(num_modes, num_modes);
            for mode = 1:num_modes % 20 (for 47tasks) and num_modes (7 representative tasks)
                basis = eigenmodes(cortex_ind, 1:mode);
                recon_beta(1:mode,mode) = calc_eigendecomposition(data_to_reconstruct(cortex_ind), basis, 'matrix');
            end

            % At parcellated level
            parc_name = 'Glasser360';
            parc = dlmread(sprintf('Data/parcellations/fsLR_32k_%s-%s.txt', parc_name, hemisphere));

            recon_mse_parc = zeros(1, num_modes);
            for mode = 1:num_modes % 20 (for 47tasks) and num_modes (7 representative tasks)
                recon_temp = eigenmodes(:, 1:mode)*recon_beta(1:mode,mode);
                recon_mse_parc(mode)     = immse(calc_parcellate(parc, data_to_reconstruct), calc_parcellate(parc, recon_temp));

            end

            if connectome_version == 1
                recon_mse_parc_tk_sbj_version1(tk,sbj,:) = recon_mse_parc;

            elseif connectome_version == 2
                recon_mse_parc_tk_sbj_version2(tk,sbj,:) = recon_mse_parc;

            elseif connectome_version == 3
                recon_mse_parc_tk_sbj_version3(tk,sbj,:) = recon_mse_parc;

            elseif connectome_version == 4
                recon_mse_parc_tk_sbj_version4(tk,sbj,:) = recon_mse_parc;
            end
        end
    end

%% for 200 modes 7 tasks
     if connectome_version == 1
        saveFileName = 'Results/tMRI_parcellated_reconstruction_7tk_subjects_EDRbinary';
        Rev2_parsave(saveFileName,...
        recon_mse_parc_tk_sbj_version1)
    elseif connectome_version == 2
        saveFileName = 'Results/tMRI_parcellated_reconstruction_7tk_subjects_EDRcontinuous';
        Rev2_parsave(saveFileName,...
        recon_mse_parc_tk_sbj_version2)
    elseif connectome_version == 3
        saveFileName = 'Results/tMRI_parcellated_reconstruction_7tk_subjects_Geometry';
        Rev2_parsave(saveFileName,...
        recon_mse_parc_tk_sbj_version3)
    elseif connectome_version == 4
        saveFileName = 'Results/tMRI_parcellated_reconstruction_7tk_subjects_EDRLR';
        Rev2_parsave(saveFileName,...
        recon_mse_parc_tk_sbj_version4)
    end
end
elseif parfor_run ==0
    representTask = [3, 11, 19, 30, 41, 44, 47]; % the 7 representative tasks

    load('Results/long_3T_task/tMRI_parcellated_reconstruction_all_subjects_EDRbinary.mat')
    recon_mse_parc_version1   = squeeze(nanmean(recon_mse_parc_tk_sbj_version,2))';

    load('Results/long_3T_task/tMRI_parcellated_reconstruction_all_subjects_EDRcontinuous.mat')
    recon_mse_parc_version2   = squeeze(nanmean(recon_mse_parc_tk_sbj_version,2))';

    load('Results/long_3T_task/tMRI_parcellated_reconstruction_all_subjects_Geometry.mat')
    recon_mse_parc_version3   = squeeze(nanmean(recon_mse_parc_tk_sbj_version,2))';

    load('Results/long_3T_task/tMRI_parcellated_reconstruction_all_subjects_EDRLR.mat')
    recon_mse_parc_version4   = squeeze(nanmean(recon_mse_parc_tk_sbj_version,2))';

    load('Results/long_3T_task/tMRI_parcellated_reconstruction_7tk_subjects_EDRbinary.mat')
    recon_mse_parc_7tk_version1   = squeeze(nanmean(recon_mse_parc_tk_sbj_version,2))';

    load('Results/long_3T_task/tMRI_parcellated_reconstruction_7tk_subjects_EDRcontinuous.mat')
    recon_mse_parc_7tk_version2   = squeeze(nanmean(recon_mse_parc_tk_sbj_version,2))';

    load('Results/long_3T_task/tMRI_parcellated_reconstruction_7tk_subjects_Geometry.mat')
    recon_mse_parc_7tk_version3   = squeeze(nanmean(recon_mse_parc_tk_sbj_version,2))';

    load('Results/long_3T_task/tMRI_parcellated_reconstruction_7tk_subjects_EDRLR.mat')
    recon_mse_parc_7tk_version4   = squeeze(nanmean(recon_mse_parc_tk_sbj_version,2))';

end

%% FIGURE A
%% representative tasks
representTask = [3, 11, 19, 30, 41, 44, 47]; % the 7 representative tasks
nMod = 200;
figure
subplot(2,4,2)
plot(1:nMod, recon_mse_parc_7tk_version1(1:nMod,representTask)./max(recon_mse_parc_7tk_version1(1:nMod,representTask)), 'b-', 'linewidth', 2);grid on; title('EDR binary');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,3);axis square
plot(1:nMod, recon_mse_parc_7tk_version2(1:nMod,representTask)./max(recon_mse_parc_7tk_version2(1:nMod,representTask)), 'g-', 'linewidth', 2');grid on; title('EDR weighted');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,1);axis square
plot(1:nMod, recon_mse_parc_7tk_version3(1:nMod,representTask)./max(recon_mse_parc_7tk_version3(1:nMod,representTask)), 'm-', 'linewidth', 2);grid on;ylim([0,1]); title('Geometry');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,4);
plot(1:nMod, recon_mse_parc_7tk_version4(1:nMod,representTask)./max(recon_mse_parc_7tk_version4(1:nMod,representTask)), 'k-', 'linewidth', 2);grid on; title('EDR+LR');xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,6);
plot((2:nMod), abs(diff([zeros(1,7); recon_mse_parc_7tk_version1(2:nMod,representTask)])), 'b-', 'linewidth', 2);grid on; title('EDR binary');%ylim([0, .1]);
xlim([4,nMod]);ylim([0, 0.3])
subplot(2,4,7);
plot((2:nMod), abs(diff([zeros(1,7); recon_mse_parc_7tk_version2(2:nMod,representTask)])), 'g-', 'linewidth', 2);grid on; title('EDR weighted');%ylim([0, .1]);
xlim([4,nMod]);ylim([0, 0.3])
subplot(2,4,5);
plot((2:nMod), abs(diff([zeros(1,7); recon_mse_parc_7tk_version3(2:nMod,representTask)])), 'm-', 'linewidth', 2);grid on; title('Geometry');%ylim([0, .1]);
xlim([4,nMod]);ylim([0, 0.3])
subplot(2,4,8);
plot((2:nMod), abs(diff([zeros(1,7); recon_mse_parc_7tk_version4(2:nMod,representTask)])), 'k-', 'linewidth', 2);grid on; title('EDR+LR');%ylim([0, .1]);
xlim([4,nMod]);ylim([0, 0.3])
exportgraphics(gcf, fullfile(fig_out_dir, 'Figure3A_TAreconstruction_7tasks.png')); % PATCH: save for headless run

%% FIGURE B
figure('Name', 'tfMRI reconstruction - accuracy difference');
th_mode = 20;
customYTicks = 1:47; % [1, 3, 5, 7, 10];
customYTickLabels = fieldNames;
% Set the custom Y-axis ticks and labels
set(gca, 'YTick', customYTicks, 'YTickLabel', customYTickLabels,'TickLabelInterpreter','none');
subplot(2,3,1)
imagesc((recon_mse_parc_version1(1:th_mode,:)-recon_mse_parc_version3(1:th_mode,:))',[-0.25 0.25]);title('EDR binary')
subplot(2,3,2)
imagesc((recon_mse_parc_version2(1:th_mode,:)-recon_mse_parc_version3(1:th_mode,:))',[-0.25 0.25]);title('EDR weighted')
subplot(2,3,3)
imagesc((recon_mse_parc_version4(1:th_mode,:)-recon_mse_parc_version3(1:th_mode,:))',[-0.25 0.25]);title('EDR+LR')

colormap(redblue)
subplot(2,3,4)
bar(mean((recon_mse_parc_version1(1:th_mode,:)-recon_mse_parc_version3(1:th_mode,:))'),'k','LineWidth',2);title('EDR binary');ylim([-0.25 0.25]);
subplot(2,3,5)
bar(mean((recon_mse_parc_version2(1:th_mode,:)-recon_mse_parc_version3(1:th_mode,:))'),'k','LineWidth',2);title('EDR weighted');ylim([-0.25 0.25]);
subplot(2,3,6)
bar(mean((recon_mse_parc_version4(1:th_mode,:)-recon_mse_parc_version3(1:th_mode,:))'),'k','LineWidth',2);title('EDR+LR');ylim([-0.25 0.25]);
exportgraphics(gcf, fullfile(fig_out_dir, 'Figure3B_TAreconstruction_47tasks_diff.png')); % PATCH: save for headless run
%% FIGURE C
% for one task differences
figure
i=7;
bar(abs(diff([0; recon_mse_parc_version4(2:th_mode,representTask(i))])), 'k', 'linewidth', 2);grid on; title('EDR+LR');ylim([0, .05]);xlim([0.5,th_mode-0.5])
title(['Mode contribution to reconstruction for task ', customYTickLabels{representTask(i)}],'Interpreter', 'none');ylabel('FC mse contribution'); xlabel('Modes')
exportgraphics(gcf, fullfile(fig_out_dir, 'Figure3C_mode_contribution.png')); % PATCH: save for headless run
%% info for FIGURE D
customYTickLabels{14}
recon_mse_parc_version3([5,10,15,20],14);
recon_mse_parc_version4([5,10,15,20],14);
