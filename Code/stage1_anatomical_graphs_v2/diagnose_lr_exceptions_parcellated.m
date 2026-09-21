%% Diagnostic v4: are LR exceptions "rare" at the parcel level?
%
% run_stage1_v2.m reproduced Vohryzek's exact fitting/exception-detection
% formula faithfully and still found 6.36M vertex-pair "exceptions"
% (1.44% of all vertex pairs) -- much higher than the word "rare" seems
% to imply. Hypothesis: the paper's characterization (and Fig. 2A's
% visualization) is at the PARCELLATED (Glasser360, 180x180) level, since
% all downstream reconstruction analysis works on parcellated data
% anyway. 6.36M vertex-pair exceptions could still collapse into a small,
% genuinely sparse number of PARCEL-pair blocks if they cluster spatially
% rather than spreading evenly across all 16,110 possible parcel pairs.
%
% This script re-runs only the exception-detection blocked pass (reusing
% the already-fitted A/lambda and bin stats from stage1_v2_report.mat --
% no need to refit), and instead of building the vertex-level adjacency
% matrix, accumulates a 180x180 parcel-pair exception count matrix.
%
% NEW FILE -- does not modify run_stage1_v2.m or any prior file.

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
project_root = fileparts(code_dir);
addpath(fullfile(code_dir, 'utils'));

data_dir = fullfile(project_root, 'data');
surface_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_midthickness-lh.vtk');
cortex_mask_file = fullfile(data_dir, 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt');
parc_file = fullfile(data_dir, 'parcellations', 'fsLR_32k_Glasser360-lh.txt');
sc_file = fullfile(data_dir, 'empirical', 'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat');
report_file = fullfile(data_dir, 'processed', 'eigenmodes_novel_v2', 'stage1_v2_report.mat');

nr_ini = 20; nr_fin = 380; nstd = 3; min_lr_distance = 40;
block_size = 2000;

fprintf('=== Diagnostic v4: parcel-level rarity of LR exceptions ===\n');

fprintf('Loading previous fit/bin-stats report...\n');
loaded = load(report_file, 'report', 'stats');
A_fit = loaded.report.A_fit; %#ok<NASGU> -- not needed here, kept for provenance
lambda_fit = loaded.report.lambda_fit; %#ok<NASGU>
max_sc = loaded.report.max_sc;
bin_edges = loaded.stats.edges;
bin_mean = loaded.stats.mean;
bin_std = loaded.stats.std;
fprintf('  reusing bin stats (400 bins) and max_sc=%.6f from prior run\n', max_sc);

fprintf('Loading surface, cortex mask, Glasser360 parcellation...\n');
[vertices, ~] = read_vtk_surface(surface_file);
cortex = dlmread(cortex_mask_file); %#ok<DLMRD>
cortex_ind = find(cortex);
vertices_cortex = vertices(cortex_ind, :);
parc_full = dlmread(parc_file); %#ok<DLMRD>
parc = parc_full(cortex_ind); % parcel label per cortical vertex, aligned with vertices_cortex
labels = unique(parc(parc > 0));
num_parcels = numel(labels);
label_to_idx = zeros(max(labels), 1);
label_to_idx(labels) = 1:num_parcels;
fprintf('  %d parcels\n', num_parcels);

fprintf('Loading + normalizing structural connectome...\n');
sc_data = load(sc_file);
C_norm = sc_data.avgSC_L / max_sc;
clear sc_data;

n = size(vertices_cortex, 1);
sq_norms = sum(vertices_cortex.^2, 2);
parcel_exception_count = zeros(num_parcels, num_parcels);
n_total_exceptions = 0;

for start_row = 1:block_size:n
    end_row = min(start_row + block_size - 1, n);
    rows = (start_row:end_row)';

    bv = vertices_cortex(rows, :);
    d2 = sq_norms(rows) - 2*(bv * vertices_cortex') + sq_norms';
    d2(d2 < 0) = 0;
    D_block = sqrt(d2);
    for r = 1:numel(rows)
        D_block(r, 1:rows(r)) = Inf;
    end

    C_block = C_norm(rows, :);
    bidx = discretize(D_block, bin_edges);
    valid_exc_bin = ~isnan(bidx) & bidx >= nr_ini & bidx <= nr_fin;
    threshold = inf(size(D_block));
    threshold(valid_exc_bin) = bin_mean(bidx(valid_exc_bin)) + nstd * bin_std(bidx(valid_exc_bin));
    is_exception = valid_exc_bin & (C_block > threshold) & (D_block > min_lr_distance);

    [rr, cc] = find(is_exception);
    if ~isempty(rr)
        gi = rows(rr); % global vertex row indices
        gj = cc;       % global vertex col indices
        pi_raw = parc(gi);
        pj_raw = parc(gj);
        % guard: a vertex inside the cortex mask can still have parcel
        % label 0 if the Glasser parcellation's medial-wall boundary
        % doesn't exactly match the cortex mask's -- indexing
        % label_to_idx(0) would otherwise error.
        has_label = pi_raw > 0 & pj_raw > 0;
        pi_idx = zeros(size(pi_raw));
        pj_idx = zeros(size(pj_raw));
        pi_idx(has_label) = label_to_idx(pi_raw(has_label));
        pj_idx(has_label) = label_to_idx(pj_raw(has_label));
        valid_parc = has_label & pi_idx > 0 & pj_idx > 0;
        pi_idx = pi_idx(valid_parc); pj_idx = pj_idx(valid_parc);
        pair_idx = sub2ind([num_parcels, num_parcels], min(pi_idx,pj_idx), max(pi_idx,pj_idx));
        counts = accumarray(pair_idx, 1, [num_parcels*num_parcels, 1]);
        parcel_exception_count = parcel_exception_count + reshape(counts, num_parcels, num_parcels);
    end
    n_total_exceptions = n_total_exceptions + nnz(is_exception);
    fprintf('  rows %d-%d done (cumulative exceptions: %d)\n', start_row, end_row, n_total_exceptions);
end

n_possible_parcel_pairs = num_parcels*(num_parcels-1)/2;
n_parcel_pairs_with_exception = nnz(triu(parcel_exception_count, 1) > 0);
pct_parcel_pairs = 100 * n_parcel_pairs_with_exception / n_possible_parcel_pairs;

fprintf('\n--- Parcel-level rarity check ---\n');
fprintf('  Total vertex-pair exceptions: %d\n', n_total_exceptions);
fprintf('  Possible parcel pairs (180 choose 2): %d\n', n_possible_parcel_pairs);
fprintf('  Parcel pairs containing >=1 exception: %d (%.2f%%)\n', n_parcel_pairs_with_exception, pct_parcel_pairs);

[sorted_counts, sorted_idx] = sort(parcel_exception_count(:), 'descend');
fprintf('\n  Top 10 parcel-pairs by exception count:\n');
[pi_top, pj_top] = ind2sub([num_parcels, num_parcels], sorted_idx(1:10));
for k = 1:10
    if pi_top(k) < pj_top(k)
        fprintf('    parcel %d - parcel %d: %d exceptions\n', pi_top(k), pj_top(k), sorted_counts(k));
    end
end

out_dir = fullfile(data_dir, 'processed', 'eigenmodes_novel_v2');
save(fullfile(out_dir, 'lr_exceptions_parcellated.mat'), 'parcel_exception_count', ...
    'n_total_exceptions', 'n_parcel_pairs_with_exception', 'pct_parcel_pairs');
fprintf('\nSaved -> %s\n', fullfile(out_dir, 'lr_exceptions_parcellated.mat'));
fprintf('=== Diagnostic complete ===\n');
