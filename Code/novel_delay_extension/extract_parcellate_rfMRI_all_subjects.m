%% Extract + parcellate resting-state rfMRI timeseries for all 255 subjects
%
% Reads raw REST1_LR CIFTI files, extracts LEFT_CORTEX grayordinates
% (dropping the medial-wall NaN padding -- validated vertex-for-vertex
% against our cortex mask in test_cifti_read_single_subject.m), parcellates
% to Glasser360, and saves a per-subject parcellated timeseries plus an
% aggregate empirical FC array.
%
% Checkpointed: each subject's parcellated timeseries is saved to its own
% .mat file, and already-completed subjects are skipped on rerun, so this
% can be safely re-launched if interrupted.
%
% MUST run on the server (255 x ~438MB raw CIFTI reads).

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir));

pang_dir = fullfile(project_root, 'Reference Code Material', 'pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
vohryzek_dir = fullfile(project_root, 'Reference Code Material', 'Vohryzeketal2024', 'vohryzek2024_EDRLR');

addpath(genpath(fullfile(pang_dir, 'functions_matlab', 'cifti-matlab-master')));
addpath(fullfile(project_root, 'Code', 'utils'));

raw_dir = fullfile(project_root, 'data', 'empirical', 'rfMRI_raw');
out_dir = fullfile(project_root, 'data', 'processed', 'rfMRI_parcellated');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% Canonical subject order + cortex/parcellation setup (matches existing pipeline)
% subject_list_HCP.txt is single-column (one subject ID per line, no index column)
subject_ids = dlmread(fullfile(project_root, 'data', 'empirical', 'subject_list_HCP.txt')); %#ok<DLMRD>
n_subj = numel(subject_ids);
fprintf('Loaded %d subject IDs. First 3: %d, %d, %d. Last: %d\n', n_subj, subject_ids(1), subject_ids(2), subject_ids(3), subject_ids(end));
if n_subj ~= 255 || subject_ids(1) ~= 100206 || subject_ids(end) ~= 877269
    error('Subject list did not parse as expected (n=%d, first=%d, last=%d). Aborting before wasting the run.', n_subj, subject_ids(1), subject_ids(end));
end

cortex = dlmread(fullfile(vohryzek_dir, 'Data', 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt')); %#ok<DLMRD>
cortex_ind = find(cortex);
parc = dlmread(fullfile(vohryzek_dir, 'Data', 'parcellations', 'fsLR_32k_Glasser360-lh.txt')); %#ok<DLMRD>
parc_cortex = parc(cortex_ind);
num_parcels = numel(unique(parc_cortex(parc_cortex>0)));
fprintf('Cortex vertices: %d, parcels: %d\n', numel(cortex_ind), num_parcels);

%% Extract + parcellate, per subject, with checkpointing
t_start = tic;
n_done = 0; n_skipped = 0; n_missing = 0; n_failed = 0;

for si = 1:n_subj
    sid = subject_ids(si);
    out_file = fullfile(out_dir, sprintf('parc_rest_%d.mat', sid));

    if exist(out_file, 'file')
        n_skipped = n_skipped + 1;
        continue;
    end

    raw_file = fullfile(raw_dir, sprintf('%d_rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii', sid));
    if ~exist(raw_file, 'file')
        fprintf('  [%d/%d] subject %d: RAW FILE MISSING, skipping.\n', si, n_subj, sid);
        n_missing = n_missing + 1;
        continue;
    end

    try
        source = ft_read_cifti(raw_file, 'readsurface', false);
        left_label_idx = find(strcmpi(source.brainstructurelabel, 'CORTEX_LEFT'));
        left_mask = (source.brainstructure == left_label_idx);
        left_timeseries = source.dtseries(left_mask, :);

        if size(left_timeseries, 1) == numel(cortex_ind)
            left_timeseries_valid = left_timeseries;
        elseif size(left_timeseries, 1) == numel(cortex)
            valid_row = ~isnan(left_timeseries(:,1));
            if ~(sum(valid_row) == numel(cortex_ind) && isequal(find(valid_row), cortex_ind(:)))
                error('Vertex alignment check failed for subject %d.', sid);
            end
            left_timeseries_valid = left_timeseries(valid_row, :);
        else
            error('Unexpected LEFT_CORTEX row count (%d) for subject %d.', size(left_timeseries,1), sid);
        end

        parc_timeseries = parcellate_average(left_timeseries_valid, parc_cortex); % [180 x T]
        T = size(parc_timeseries, 2);

        save(out_file, 'parc_timeseries', 'sid', 'T', '-v7.3');
        n_done = n_done + 1;
    catch ME
        fprintf('  [%d/%d] subject %d: FAILED (%s)\n', si, n_subj, sid, ME.message);
        n_failed = n_failed + 1;
        continue;
    end

    if mod(n_done, 10) == 0 || si == n_subj
        elapsed = toc(t_start);
        fprintf('  [%d/%d] done=%d skipped=%d missing=%d failed=%d, elapsed=%.1fs\n', ...
            si, n_subj, n_done, n_skipped, n_missing, n_failed, elapsed);
    end
end

fprintf('\nExtraction pass complete: done=%d, skipped(already existed)=%d, missing=%d, failed=%d\n', ...
    n_done, n_skipped, n_missing, n_failed);

%% Aggregate: build empirical FC for every successfully extracted subject
fprintf('\nAggregating per-subject FC...\n');
FC_all = nan(num_parcels, num_parcels, n_subj);
T_all = nan(n_subj, 1);
valid_subj = false(n_subj, 1);

for si = 1:n_subj
    sid = subject_ids(si);
    out_file = fullfile(out_dir, sprintf('parc_rest_%d.mat', sid));
    if ~exist(out_file, 'file'); continue; end
    d = load(out_file, 'parc_timeseries', 'T');
    FC_all(:,:,si) = corr(d.parc_timeseries');
    T_all(si) = d.T;
    valid_subj(si) = true;
end

fprintf('Aggregated FC for %d/%d subjects.\n', sum(valid_subj), n_subj);
save(fullfile(out_dir, 'FC_all_subjects.mat'), 'FC_all', 'T_all', 'valid_subj', 'subject_ids', '-v7.3');
fprintf('Saved -> %s\n', fullfile(out_dir, 'FC_all_subjects.mat'));
fprintf('\n=== Full extraction + FC aggregation complete ===\n');
