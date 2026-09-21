% diag_lr_magnitude.m
%
% Diagnostic script: how big are the spliced-in long-range (LR) connections
% actually, compared to (a) the smooth EDR baseline they replace, and
% (b) typical short-range connections elsewhere in the connectome?
%
% This was written because step4_5 found no real difference between EDR+LR
% and EDR-continuous. One possible explanation: the long-range connections
% are "exceptions" only in a statistical sense (rare, >3 SD above the
% local mean) but numerically tiny, so they'd barely register in the
% eigenmodes even though they pass the paper's detection rule.
%
% Re-derives which connections count as long-range exceptions by replaying
% connectome_harmonics.m Stage 2's binning/detection logic exactly, then
% reports how large those connections actually are (rather than just
% counting them).

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Code from Vohryzek et al 2024', 'vohryzek2024_EDRLR');
addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

hemisphere = 'lh';
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
vertices_cortex = surface_midthickness.vertices(cortex_ind, :);

fprintf('Loading empirical connectome...\n');
load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', ...
    'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');
C = avgSC_L / max(avgSC_L(:));
clear avgSC_L

fprintf('Loading checkpointed Afit/lambda...\n');
load(fullfile('results', 'EDR_LR_connectome_full_resolution.mat'), 'Afit', 'lambda');

tic
rr = squareform(pdist(single(vertices_cortex)));
fprintf('Computed distance matrix in %.1f s\n', toc);

% replay Stage 2's exact binning/exception-detection logic
NR = 400; NRini = 20; NRfin = 380; NSTD = 3; DistRange = 40;
range_dist = max(rr(:));
delta = range_dist / NR;
index = uint16(min(floor(double(rr)/delta) + 1, NR));

tic
bin_mean = accumarray(index(:), C(:), [NR,1], @mean, single(NaN));
bin_std  = accumarray(index(:), C(:), [NR,1], @std, single(NaN));
fprintf('Binned in %.1f s\n', toc);

local_mean = bin_mean(index);
local_std = bin_std(index);
is_exception = (index >= NRini) & (index <= NRfin) & (rr > DistRange) & (C > local_mean + NSTD*local_std);
fprintf('Exceptions found: %d (upper triangle)\n', nnz(triu(is_exception,1)));

%% Compare exception weight vs. the EDR baseline they replace

actual_at_exception   = double(C(is_exception));
baseline_at_exception = Afit(1) * exp(-Afit(2) * double(rr(is_exception)));

% typical short-range weight, for scale: nonzero connections within the
% same distance window used to define "long-range" (<=40mm)
short_range_mask = (rr <= DistRange) & (rr > 0);
typical_short_range = double(C(short_range_mask));

fprintf('\n=== Magnitude comparison ===\n');
fprintf('Baseline EDR value at exception positions:  median=%.6g  mean=%.6g  max=%.6g\n', ...
    median(baseline_at_exception), mean(baseline_at_exception), max(baseline_at_exception));
fprintf('Actual (spliced) value at exception positions: median=%.6g  mean=%.6g  max=%.6g\n', ...
    median(actual_at_exception), mean(actual_at_exception), max(actual_at_exception));
fprintf('Ratio (actual/baseline) at exception positions: median=%.3g  mean=%.3g\n', ...
    median(actual_at_exception./baseline_at_exception), mean(actual_at_exception./baseline_at_exception));
fprintf('\nTypical short-range (<=%dmm) weight in same connectome: median=%.6g  mean=%.6g  max=%.6g\n', ...
    DistRange, median(typical_short_range), mean(typical_short_range), max(typical_short_range));

fprintf('\n=== Total edge-weight mass comparison ===\n');
total_mass_all = sum(double(C(:)));
total_mass_exception_added = sum(actual_at_exception - baseline_at_exception);
fprintf('Total weight mass of whole connectome C: %.6g\n', total_mass_all);
fprintf('Total EXTRA weight mass added by LR splice (actual - baseline, summed): %.6g (%.4g%% of total)\n', ...
    total_mass_exception_added, 100*total_mass_exception_added/total_mass_all);

fprintf('\n=== Per-vertex degree comparison ===\n');
d_baseline_only = sum(Afit(1)*exp(-Afit(2)*double(rr)), 2) - Afit(1); % minus self-term
extra = zeros(size(rr,1),1,'double');
[ii,~] = find(is_exception);
contrib = actual_at_exception - baseline_at_exception;
extra = accumarray(ii, contrib, [size(rr,1),1]);
fprintf('Median vertex degree from smooth EDR baseline alone: %.6g\n', median(d_baseline_only));
fprintf('Median EXTRA degree contributed by LR splice: %.6g\n', median(extra));
fprintf('Median relative degree increase from LR splice: %.4g%%\n', ...
    100*median(extra ./ d_baseline_only));
