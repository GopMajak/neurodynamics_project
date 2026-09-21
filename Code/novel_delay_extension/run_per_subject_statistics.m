%% Per-subject reconstruction + paired statistics (matches Vohryzek's actual design)
%
% Everything so far (Stage 8, multi-speed sweep, eigenbasis sweep) used
% GROUP-AVERAGED task maps -- a first-pass simplification. This instead
% reconstructs each of the 255 subjects individually and runs the paired
% comparisons the thesis proposal specifies at N=20 (its stated primary
% comparison point):
%   1. Delayed EDR+LR    vs. static EDR+LR
%   2. Delayed EDR-only  vs. static EDR-only
%   3. Delayed EDR-only  vs. static EDR+LR
%   4. Delayed EDR+LR    vs. delayed EDR-only
%
% "Delayed" = the W_delay eigenvector basis (the strongest delay-informed
% method found so far -- see build_delay_eigenvector_basis.m), evaluated
% at the STANDARD reference point used throughout this session
% (f_mu=10Hz, speed=5 m/s, kappa=6 -- kappa is provably irrelevant to
% this basis, verified in sweep_eigbasis_fmu.m), NOT the best point found
% in that parameter sweep -- testing significance at a cherry-picked
% "best" point would be circular/selection-biased.
%
% Bonferroni-corrected paired t-tests, matching Vohryzek's own approach
% (PNAS_Figure_2_main.m multiplies raw p-values by the number of
% comparisons).

clear; clc;

%% Paths
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir));

pang_dir = fullfile(project_root, 'Reference Code Material', 'pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
vohryzek_dir = fullfile(project_root, 'Reference Code Material', 'Vohryzeketal2024', 'vohryzek2024_EDRLR');
recovered_dir = fullfile(project_root, 'Original Project Code', 'brain_connectome_harmonics');

addpath(genpath(fullfile(pang_dir, 'functions_matlab')));
addpath(fullfile(recovered_dir, 'functions'));
addpath(fullfile(project_root, 'Code', 'utils'));

fig_out_dir = fullfile(project_root, 'Figures', 'novel_delay_extension');
results_out_dir = fullfile(project_root, 'Results', 'novel_delay_extension');
if ~exist(fig_out_dir, 'dir'); mkdir(fig_out_dir); end
if ~exist(results_out_dir, 'dir'); mkdir(results_out_dir); end

%% Rebuild parcel-level EDR-only / EDR+LR graphs
hemisphere = 'lh'; surface_interest = 'fsLR_32k'; mesh_interest = 'midthickness'; parc_name = 'Glasser360';
[vertices, faces] = read_vtk(fullfile(vohryzek_dir, 'Data', 'template_surfaces', sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices'; surface_midthickness.faces = faces';
cortex = dlmread(fullfile(vohryzek_dir, 'Data', 'template_surfaces', sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
parc = dlmread(fullfile(vohryzek_dir, 'Data', 'parcellations', sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
parcels = unique(parc_cortex(parc_cortex>0));
num_parcels = length(parcels);
fprintf('Loaded %s parcellation: %d parcels\n', parc_name, num_parcels);

load(fullfile(pang_dir, 'data', 'empirical', 'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');
vertices_cortex = surface_midthickness.vertices(cortex_ind, :);
centroids = zeros(num_parcels, 3);
for p = 1:num_parcels
    centroids(p,:) = mean(vertices_cortex(parc_cortex==parcels(p), :), 1);
end
rr_parc = squareform(pdist(centroids));
connectome_parc = calc_parcellate_matrix(parc_cortex, avgSC_L);
C_parc = connectome_parc / max(connectome_parc(:));
clear avgSC_L

NR = 60; NSTD = 3; DistRange = 40;
range_dist = max(rr_parc(:)); delta = range_dist / NR;
xcoor = delta/2 + delta*(0:NR-1);
index_parc = floor(rr_parc/delta) + 1; index_parc(index_parc > NR) = NR;
sc_density = cell(1, NR); sc_density_i = cell(1, NR); sc_density_j = cell(1, NR);
ycoor2 = nan(1, NR);
for n = 1:NR
    [idx_i, idx_j] = find(index_parc == n);
    idx = find(index_parc == n);
    sc_density{n} = C_parc(idx); sc_density_i{n} = idx_i; sc_density_j{n} = idx_j;
    if ~isempty(idx); ycoor2(n) = mean(C_parc(idx)); end
end
fit_start = find(xcoor >= 10, 1);
fit_ind = fit_start:NR; fit_ind = fit_ind(~isnan(ycoor2(fit_ind)));
expfunc = @(A, x) (A(1)*exp(-A(2)*x));
sse = @(A) sum((expfunc(A, xcoor(fit_ind)) - ycoor2(fit_ind)).^2);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, 'Display', 'off');
Afit = fminsearch(sse, [0.15, 0.18], options);
Clong = zeros(num_parcels, num_parcels);
for i = fit_start:NR
    if isempty(sc_density{i}); continue; end
    mv = mean(sc_density{i}); st = std(sc_density{i});
    ind_exc = find(sc_density{i} > mv + NSTD*st);
    for n = 1:numel(ind_exc)
        ii = sc_density_i{i}(ind_exc(n)); jj = sc_density_j{i}(ind_exc(n));
        if rr_parc(ii,jj) > DistRange
            Clong(ii,jj) = sc_density{i}(ind_exc(n));
        end
    end
end
EDR_conn = Afit(1)*exp(-Afit(2)*rr_parc);
EDR_LRE_parc = EDR_conn; EDR_LRE_parc(Clong>0) = Clong(Clong>0);

N3 = num_parcels;
EDR_conn_km = EDR_conn;     EDR_conn_km(1:N3+1:end) = 0;
EDR_LRE_km  = EDR_LRE_parc; EDR_LRE_km(1:N3+1:end)  = 0;

%% Build the four bases
fprintf('\nComputing static (Laplacian) eigenmodes...\n');
[psi_static_edr, ~]   = normalized_laplacian_eigenmodes(EDR_conn_km, num_parcels);
[psi_static_edrlr, ~] = normalized_laplacian_eigenmodes(EDR_LRE_km, num_parcels);

fprintf('Computing W_delay eigenvector basis (f_mu=10Hz, speed=5m/s, kappa=6 -- standard reference point)...\n');
dt = 1e-3; f_mu = 10; omega = f_mu*2*pi*ones(N3,1); kappa = 6; speed_basis = 5;
tau = delay_matrix_from_distance(rr_parc, speed_basis, dt);
[psi_delay_edr, ~] = build_delay_eigenvector_basis(EDR_conn_km, tau, omega, kappa, dt);
[psi_delay_edrlr, ~] = build_delay_eigenvector_basis(EDR_LRE_km, tau, omega, kappa, dt);

conditions = struct( ...
    'name', {'static_EDR', 'static_EDRLR', 'delay_EDR', 'delay_EDRLR'}, ...
    'basis', {psi_static_edr, psi_static_edrlr, psi_delay_edr, psi_delay_edrlr});
n_cond = numel(conditions);

%% Per-subject reconstruction, all 7 representative tasks, N=1..180
fprintf('\nLoading task data and reconstructing per-subject (255 subjects x 7 tasks x %d conditions)...\n', n_cond);
data = load(fullfile(vohryzek_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
fieldNames = fieldnames(data.zstat);
representTask = [3, 11, 19, 30, 41, 44, 47];
n_tasks = numel(representTask);
max_N = num_parcels;
N_target = 20; % the proposal's primary comparison point

acc20 = nan(255, n_tasks, n_cond); % per-subject accuracy at N=20, per task, per condition

for ti = 1:n_tasks
    fieldName = fieldNames{representTask(ti)};
    Y = parcellate_average(data.zstat.(fieldName)(cortex_ind, :), parc_cortex); % [180 x 255]
    n_subj = size(Y, 2);

    for ci = 1:n_cond
        Psi = conditions(ci).basis;
        A = Psi' * Y; % [180 x 255] projection coefficients, all modes at once
        Y_hat_N = Psi(:, 1:N_target) * A(1:N_target, :); % reconstruction using first N_target modes
        acc20(1:n_subj, ti, ci) = corr_columns(Y, Y_hat_N)';
    end
    fprintf('  task %d/%d (%s) done\n', ti, n_tasks, fieldName);
end

% average across the 7 tasks -> one accuracy value per subject per condition
acc20_avg = squeeze(mean(acc20, 2)); % [255 x 4]

%% Paired statistics (Bonferroni-corrected, matching Vohryzek's own approach)
labels = {'static EDR', 'static EDR+LR', 'delay EDR', 'delay EDR+LR'};
comparisons = {
    'delay EDR+LR vs static EDR+LR', 4, 2;
    'delay EDR-only vs static EDR-only', 3, 1;
    'delay EDR-only vs static EDR+LR', 3, 2;
    'delay EDR+LR vs delay EDR-only', 4, 3;
};
n_comp = size(comparisons, 1);

n_nan_subjects = sum(any(isnan(acc20_avg), 2));
if n_nan_subjects > 0
    fprintf('\nNote: %d subject(s) have a NaN accuracy value in at least one condition (likely a zero-variance\n', n_nan_subjects);
    fprintf('activation map for that subject/task -- correlation is undefined). Excluded pairwise per comparison,\n');
    fprintf('matching how ttest itself handles them.\n');
end

fprintf('\n=== Paired statistics at N=%d (255 subjects, 7-task average), Bonferroni-corrected x%d ===\n', N_target, n_comp);
stats_results = cell(n_comp, 6);
for k = 1:n_comp
    name = comparisons{k,1};
    ia = comparisons{k,2}; ib = comparisons{k,3};
    x = acc20_avg(:, ia); y = acc20_avg(:, ib);
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid); y = y(valid);
    [~, p, ~, st] = ttest(x, y);
    p_corrected = min(p * n_comp, 1);
    d = mean(x - y) / std(x - y); % Cohen's d for paired samples
    fprintf('  %-38s: mean(%s)=%.4f, mean(%s)=%.4f, diff=%.4f, t(%d)=%.3f, p=%.2e, p_corr=%.2e, d=%.3f, n=%d\n', ...
        name, labels{ia}, mean(x), labels{ib}, mean(y), mean(x-y), st.df, st.tstat, p, p_corrected, d, numel(x));
    stats_results(k,:) = {name, mean(x), mean(y), st.tstat, p_corrected, d};
end

%% Figure: per-subject accuracy@N=20 boxplot, matching PNAS Figure 2C style
fg = figure('Name', 'Per-subject accuracy at N=20', 'Visible', 'off');
boxplot(acc20_avg, 'Labels', labels);
hold on
c = [0.3 0.3 0.9; 0.9 0.3 0.3; 0.3 0.7 0.3; 0.9 0.6 0.1];
for ci = 1:n_cond
    xi = ci + (rand(255,1)-0.5)*0.15;
    scatter(xi, acc20_avg(:,ci), 6, c(ci,:), 'filled', 'MarkerFaceAlpha', 0.4);
end
ylabel(sprintf('Correlation with empirical activity (N=%d)', N_target));
title(sprintf('Per-subject reconstruction accuracy at N=%d (255 subjects, 7-task average)', N_target));
grid on
exportgraphics(fg, fullfile(fig_out_dir, 'per_subject_accuracy_N20_boxplot.png'));

%% Save
save(fullfile(results_out_dir, 'per_subject_statistics_results.mat'), ...
    'acc20', 'acc20_avg', 'labels', 'stats_results', 'N_target', 'f_mu', 'speed_basis', 'kappa', '-v7.3');
fprintf('\nSaved figure -> %s\n', fullfile(fig_out_dir, 'per_subject_accuracy_N20_boxplot.png'));
fprintf('Saved results -> %s\n', fullfile(results_out_dir, 'per_subject_statistics_results.mat'));
fprintf('=== Per-subject statistics complete ===\n');
