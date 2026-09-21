%% Feasibility test: read ONE subject's raw rfMRI CIFTI file, extract+parcellate
%
% Before committing to a 255-subject extraction pipeline, this tests on a
% single subject: how long does ft_read_cifti actually take on a ~438MB
% file, does the LEFT_CORTEX grayordinate count match our cortex mask
% (29,696, the key alignment sanity check), and does parcellation to
% Glasser360 produce a sane result.
%
% MUST run on the server (real 438MB file read + full pipeline test).

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir));

pang_dir = fullfile(project_root, 'Reference Code Material', 'pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
vohryzek_dir = fullfile(project_root, 'Reference Code Material', 'Vohryzeketal2024', 'vohryzek2024_EDRLR');

addpath(genpath(fullfile(pang_dir, 'functions_matlab', 'cifti-matlab-master')));
addpath(fullfile(project_root, 'Code', 'utils'));

%% Locate one subject's file
raw_dir = fullfile(project_root, 'data', 'empirical', 'rfMRI_raw');
files = dir(fullfile(raw_dir, '*.dtseries.nii'));
if isempty(files)
    error('No .dtseries.nii files found in %s', raw_dir);
end
test_file = fullfile(raw_dir, files(1).name);
fprintf('Testing on: %s\n', files(1).name);
fprintf('File size: %.1f MB\n', files(1).bytes / 1e6);

%% Read the CIFTI file, timed
fprintf('\nReading CIFTI file (ft_read_cifti, readsurface=false)...\n');
tic;
source = ft_read_cifti(test_file, 'readsurface', false);
t_read = toc;
fprintf('Read time: %.1f s\n', t_read);

fprintf('\nStructure fields: %s\n', strjoin(fieldnames(source), ', '));
if isfield(source, 'brainstructurelabel')
    fprintf('Brain structures found: %s\n', strjoin(source.brainstructurelabel, ', '));
end

%% Extract LEFT_CORTEX
left_label_idx = find(strcmpi(source.brainstructurelabel, 'CORTEX_LEFT'));
if isempty(left_label_idx)
    error('CORTEX_LEFT not found in brainstructurelabel');
end
left_mask = (source.brainstructure == left_label_idx);
fprintf('\nLEFT_CORTEX grayordinates: %d\n', sum(left_mask));

% figure out which data field holds the timeseries
if isfield(source, 'dtseries')
    data_field = 'dtseries';
elseif isfield(source, 'data')
    data_field = 'data';
else
    error('Could not find a timeseries data field (expected dtseries or data). Fields: %s', strjoin(fieldnames(source), ', '));
end
timeseries_full = source.(data_field);
fprintf('Full data field ''%s'' size: %d x %d\n', data_field, size(timeseries_full,1), size(timeseries_full,2));

left_timeseries = timeseries_full(left_mask, :);
fprintf('Left-cortex timeseries size: %d x %d\n', size(left_timeseries,1), size(left_timeseries,2));

%% Compare against our cortex mask count (key alignment check)
cortex = dlmread(fullfile(vohryzek_dir, 'Data', 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt')); %#ok<DLMRD>
cortex_ind = find(cortex);
fprintf('\nOur cortex mask (fsLR_32k_cortex-lh_mask.txt) non-medial-wall count: %d\n', numel(cortex_ind));

alignment_ok = false;
if sum(left_mask) == numel(cortex_ind)
    fprintf('MATCH: CIFTI LEFT_CORTEX count equals our cortex mask count -- consistent with standard fs_LR 32k grayordinate ordering.\n');
    alignment_ok = true;
    left_timeseries_valid = left_timeseries;
elseif sum(left_mask) == numel(cortex)
    % CIFTI padded CORTEX_LEFT out to the FULL surface (medial wall included,
    % filled with NaN placeholder rows). Recover the true grayordinate set by
    % dropping NaN rows, then check it lines up with our mask vertex-for-vertex.
    fprintf('CIFTI LEFT_CORTEX count (%d) equals FULL surface vertex count -- looks like medial-wall padding.\n', sum(left_mask));
    valid_row = ~isnan(left_timeseries(:,1));
    fprintf('Non-NaN rows within LEFT_CORTEX: %d (expect %d)\n', sum(valid_row), numel(cortex_ind));
    if sum(valid_row) == numel(cortex_ind) && isequal(find(valid_row), cortex_ind(:))
        fprintf('MATCH (vertex-for-vertex): non-NaN row indices exactly equal our cortex mask indices.\n');
        alignment_ok = true;
        left_timeseries_valid = left_timeseries(valid_row, :);
    else
        fprintf('MISMATCH -- non-NaN row set does not exactly equal our cortex mask index set. Do not proceed without resolving this.\n');
    end
else
    fprintf('MISMATCH -- CIFTI count (%d) matches neither our cortex mask count (%d) nor the full surface count (%d). Do not proceed without resolving this.\n', ...
        sum(left_mask), numel(cortex_ind), numel(cortex));
end

%% Parcellate to Glasser360 and sanity-check the result
parc = dlmread(fullfile(vohryzek_dir, 'Data', 'parcellations', 'fsLR_32k_Glasser360-lh.txt')); %#ok<DLMRD>
parc_cortex = parc(cortex_ind);

if alignment_ok
    parc_timeseries = parcellate_average(left_timeseries_valid, parc_cortex); % [180 x T]
    fprintf('\nParcellated timeseries size: %d x %d\n', size(parc_timeseries,1), size(parc_timeseries,2));
    fprintf('Sample stats -- mean: %.4f, std: %.4f, any NaN: %d\n', ...
        mean(parc_timeseries(:)), std(parc_timeseries(:)), any(isnan(parc_timeseries(:))));

    % quick FC sanity check: parcel 1 should correlate more with nearby/related
    % parcels than a generic distant one -- just check the FC matrix looks
    % like a real correlation matrix (values in [-1,1], diagonal ~1)
    FC = corr(parc_timeseries');
    fprintf('FC matrix: diagonal mean = %.4f (expect ~1), off-diag range = [%.3f, %.3f]\n', ...
        mean(diag(FC)), min(FC(~eye(size(FC)))), max(FC(~eye(size(FC)))));
else
    fprintf('\nSkipping parcellation due to count mismatch above.\n');
end

fprintf('\n=== Single-subject feasibility test complete ===\n');
