%% Hybrid basis: does a fixed delay:static RATIO help across the whole N range?
%
% Follow-up to run_hybrid_basis_test.m, which fixed N_total=20 and swept
% the absolute number of delay modes -- found a modest, significant
% improvement at N_delay=3,6,7 (~15-35% delay fraction). This instead
% fixes several delay:static RATIOS and sweeps them across the full
% N_total=1..180 range, to check whether a modest delay fraction helps
% consistently (not just a coincidence of N_total=20), and whether the
% best ratio changes with the total mode budget.

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

%% Rebuild parcel-level EDR+LR graph
hemisphere = 'lh'; surface_interest = 'fsLR_32k'; mesh_interest = 'midthickness'; parc_name = 'Glasser360';
[vertices, faces] = read_vtk(fullfile(vohryzek_dir, 'Data', 'template_surfaces', sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices'; surface_midthickness.faces = faces';
cortex = dlmread(fullfile(vohryzek_dir, 'Data', 'template_surfaces', sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
parc = dlmread(fullfile(vohryzek_dir, 'Data', 'parcellations', sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
parcels = unique(parc_cortex(parc_cortex>0));
num_parcels = length(parcels);

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
EDR_LRE_km = EDR_LRE_parc; EDR_LRE_km(1:N3+1:end) = 0;

%% Build static + delay-informed bases
fprintf('Computing static (Laplacian) EDR+LR eigenmodes...\n');
[psi_static, ~] = normalized_laplacian_eigenmodes(EDR_LRE_km, num_parcels);

fprintf('Computing W_delay eigenvector basis (f_mu=10Hz, speed=5m/s, kappa=6)...\n');
dt = 1e-3; f_mu = 10; omega = f_mu*2*pi*ones(N3,1); kappa = 6; speed_basis = 5;
tau = delay_matrix_from_distance(rr_parc, speed_basis, dt);
[psi_delay, ~] = build_delay_eigenvector_basis(EDR_LRE_km, tau, omega, kappa, dt);

%% Load per-subject task data
fprintf('Loading task data...\n');
data = load(fullfile(vohryzek_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
fieldNames = fieldnames(data.zstat);
representTask = [3, 11, 19, 30, 41, 44, 47];
n_tasks = numel(representTask);

Y_tasks = cell(1, n_tasks);
for ti = 1:n_tasks
    Y_tasks{ti} = parcellate_average(data.zstat.(fieldNames{representTask(ti)})(cortex_ind, :), parc_cortex);
end
n_subj = size(Y_tasks{1}, 2);

%% Sweep: several fixed delay:static ratios, across N_total = 1..40 (fine) + a few larger anchors
N_total_values = [1:1:40, 60, 80, 100, 140, 180];
n_N = numel(N_total_values);
ratios = [0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.00];
n_ratios = numel(ratios);

mean_acc = nan(n_ratios, n_N);
sem_acc = nan(n_ratios, n_N);

fprintf('\n=== Ratio sweep: %d ratios x %d N_total values ===\n', n_ratios, n_N);
for ri = 1:n_ratios
    ratio = ratios(ri);
    for ni = 1:n_N
        N_total = N_total_values(ni);
        N_delay = round(ratio * N_total);
        N_static = N_total - N_delay;

        if N_static > 0 && N_delay > 0
            candidate = [psi_static(:,1:N_static), psi_delay(:,1:N_delay)];
            [Psi, ~] = qr(candidate, 0);
        elseif N_static > 0
            Psi = psi_static(:,1:N_static);
        else
            Psi = psi_delay(:,1:N_delay);
        end

        acc_this = nan(n_subj, n_tasks);
        for ti = 1:n_tasks
            Y = Y_tasks{ti};
            A = Psi' * Y;
            Y_hat = Psi * A;
            acc_this(:,ti) = corr_columns(Y, Y_hat)';
        end
        acc_avg = mean(acc_this, 2);
        mean_acc(ri, ni) = mean(acc_avg, 'omitnan');
        sem_acc(ri, ni) = std(acc_avg, 0, 'omitnan') / sqrt(sum(~isnan(acc_avg)));
    end
    fprintf('  ratio=%.2f done (mean accuracy at N=20: %.4f)\n', ratio, mean_acc(ri, N_total_values==20));
end

%% Figure: full curves for each ratio
fg = figure('Name', 'Delay:static ratio sweep - full curves', 'Visible', 'off');
colors = lines(n_ratios);
hold on
for ri = 1:n_ratios
    plot(N_total_values, mean_acc(ri,:), 'LineWidth', 1.5, 'Color', colors(ri,:));
end
grid on; xlabel('Total number of modes (N)'); ylabel('Mean per-subject accuracy (7-task average)');
legend(arrayfun(@(r) sprintf('%.0f%% delay', r*100), ratios, 'UniformOutput', false), 'Location', 'southeast');
title('Reconstruction accuracy vs. N, for fixed delay:static ratios (EDR+LR)');
exportgraphics(fg, fullfile(fig_out_dir, 'hybrid_ratio_sweep_full.png'));

%% Zoomed view, N<=40 (where the effect is most relevant)
fg2 = figure('Name', 'Delay:static ratio sweep - zoomed N<=40', 'Visible', 'off');
hold on
zoom_idx = N_total_values <= 40;
for ri = 1:n_ratios
    plot(N_total_values(zoom_idx), mean_acc(ri,zoom_idx), 'LineWidth', 1.5, 'Color', colors(ri,:));
end
grid on; xlabel('Total number of modes (N)'); ylabel('Mean per-subject accuracy (7-task average)');
legend(arrayfun(@(r) sprintf('%.0f%% delay', r*100), ratios, 'UniformOutput', false), 'Location', 'southeast');
title('Zoomed (N<=40): reconstruction accuracy vs. N, for fixed delay:static ratios');
exportgraphics(fg2, fullfile(fig_out_dir, 'hybrid_ratio_sweep_zoomed.png'));

%% Figure: improvement over pure static (ratio=0), per ratio, vs N
fg3 = figure('Name', 'Improvement over pure static', 'Visible', 'off');
hold on
for ri = 2:n_ratios
    plot(N_total_values(zoom_idx), mean_acc(ri,zoom_idx) - mean_acc(1,zoom_idx), 'LineWidth', 1.5, 'Color', colors(ri,:));
end
yline(0, 'k--');
grid on; xlabel('Total number of modes (N)'); ylabel('Accuracy minus pure-static accuracy');
legend(arrayfun(@(r) sprintf('%.0f%% delay', r*100), ratios(2:end), 'UniformOutput', false), 'Location', 'best');
title('Improvement (or cost) of adding delay modes, relative to pure static EDR+LR');
exportgraphics(fg3, fullfile(fig_out_dir, 'hybrid_ratio_improvement.png'));

%% Save
save(fullfile(results_out_dir, 'hybrid_ratio_sweep_results.mat'), ...
    'N_total_values', 'ratios', 'mean_acc', 'sem_acc', '-v7.3');
fprintf('\nSaved figures -> %s\n', fig_out_dir);
fprintf('Saved results -> %s\n', fullfile(results_out_dir, 'hybrid_ratio_sweep_results.mat'));
fprintf('=== Ratio sweep complete ===\n');
