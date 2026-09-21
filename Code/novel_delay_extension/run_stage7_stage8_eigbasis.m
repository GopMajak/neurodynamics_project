%% Stage 7 + 8, W_delay eigenvector basis: theoretically-grounded alternative
%
% Both prior delay-informed basis attempts (SVD of raw cos(theta),
% SVD of carrier-removed cos(phi)) required either accepting a degenerate
% ~2-dimensional basis (identical frequencies) or deviating from
% Budzinski's exact model (frequency heterogeneity) to get a
% multi-dimensional basis at all -- and even then, static eigenmodes
% won a robust head-to-head across the entire speed range.
%
% This tries a different construction entirely: instead of PCA of
% SIMULATED dynamics, build the basis directly from the EIGENVECTORS of
% the complex delay-rotated coupling operator W_delay (the same object
% already used for the gamma-prediction diagnostics and the ring-graph
% validation). No simulation needed -- pure linear algebra -- and no
% deviation from identical frequencies (sigma_f=0, Budzinski's exact
% setup). See build_delay_eigenvector_basis.m for how the complex,
% non-orthogonal eigenvectors are turned into a real orthonormal basis.
%
% Cheap enough (no simulation) to test across the full speed range
% directly, unlike the earlier multi-speed sweep which needed a full
% Kuramoto simulation per speed.

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

%% Static baseline
fprintf('\nComputing static (Laplacian) eigenmodes at parcel resolution...\n');
[psi_static_edr, ~]   = normalized_laplacian_eigenmodes(EDR_conn_km, num_parcels);
[psi_static_edrlr, ~] = normalized_laplacian_eigenmodes(EDR_LRE_km, num_parcels);

data = load(fullfile(vohryzek_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
fieldNames = fieldnames(data.zstat);
representTask = [3, 11, 19, 30, 41, 44, 47];
n_tasks = numel(representTask);
max_N = num_parcels;

y_tasks = zeros(num_parcels, n_tasks);
for ti = 1:n_tasks
    activation_map = nanmean(data.zstat.(fieldNames{representTask(ti)}), 2);
    y_tasks(:,ti) = parcellate_average(activation_map(cortex_ind), parc_cortex);
end

[~, corr_static_edr] = reconstruct_all(psi_static_edr, y_tasks, max_N);
[~, corr_static_edrlr] = reconstruct_all(psi_static_edrlr, y_tasks, max_N);
mean_corr_static_edr = mean(corr_static_edr, 1);
mean_corr_static_edrlr = mean(corr_static_edrlr, 1);
fprintf('Static baseline Accuracy@N=20 -- EDR: %.3f, EDR+LR: %.3f\n', mean_corr_static_edr(20), mean_corr_static_edrlr(20));

%% W_delay eigenvector basis, across the same speed grid as before
dt = 1e-3;
f_mu = 10; omega = f_mu*2*pi*ones(N3,1); % IDENTICAL frequencies -- no model deviation needed
kappa = 6;
speeds = linspace(1, 30, 20);
n_speeds = numel(speeds);

acc20_eig_edr = nan(1, n_speeds);
acc20_eig_edrlr = nan(1, n_speeds);

fprintf('\n=== W_delay eigenvector basis across speeds ===\n');
for si = 1:n_speeds
    speed = speeds(si);
    tau = delay_matrix_from_distance(rr_parc, speed, dt);

    [psi_edr, ~] = build_delay_eigenvector_basis(EDR_conn_km, tau, omega, kappa, dt);
    [psi_edrlr, ~] = build_delay_eigenvector_basis(EDR_LRE_km, tau, omega, kappa, dt);

    [~, corr_edr] = reconstruct_all(psi_edr, y_tasks, max_N);
    [~, corr_edrlr] = reconstruct_all(psi_edrlr, y_tasks, max_N);
    mc_edr = mean(corr_edr, 1);
    mc_edrlr = mean(corr_edrlr, 1);
    acc20_eig_edr(si) = mc_edr(20);
    acc20_eig_edrlr(si) = mc_edrlr(20);

    fprintf('  speed=%5.2f m/s | Acc@N=20 EDR=%.3f EDR+LR=%.3f\n', speed, acc20_eig_edr(si), acc20_eig_edrlr(si));
end

%% Figure
fg = figure('Name', 'W_delay eigenvector basis reconstruction accuracy', 'Visible', 'off');
plot(speeds, acc20_eig_edr, 'o-', 'LineWidth', 1.8); hold on
plot(speeds, acc20_eig_edrlr, 's-', 'LineWidth', 1.8);
yline(mean_corr_static_edr(20), '--', 'Color', [0.3 0.3 0.9], 'LineWidth', 1.5);
yline(mean_corr_static_edrlr(20), '--', 'Color', [0.9 0.3 0.3], 'LineWidth', 1.5);
grid on; xlabel('Conduction speed (m/s)'); ylabel('Correlation with empirical activity (N=20 modes)');
legend({'W_{delay} eigenbasis EDR', 'W_{delay} eigenbasis EDR+LR', 'Static EDR (reference)', 'Static EDR+LR (reference)'}, 'Location', 'best');
title('Reconstruction accuracy at N=20 vs. speed (W_{delay} eigenvector basis)');
exportgraphics(fg, fullfile(fig_out_dir, 'eigbasis_accuracy_N20.png'));

%% Full curves at the reference speed (5 m/s) for a direct like-for-like comparison
[~, ref_idx] = min(abs(speeds - 5));
tau_ref = delay_matrix_from_distance(rr_parc, speeds(ref_idx), dt);
[psi_edr_ref, ~] = build_delay_eigenvector_basis(EDR_conn_km, tau_ref, omega, kappa, dt);
[psi_edrlr_ref, ~] = build_delay_eigenvector_basis(EDR_LRE_km, tau_ref, omega, kappa, dt);
[~, corr_edr_ref] = reconstruct_all(psi_edr_ref, y_tasks, max_N);
[~, corr_edrlr_ref] = reconstruct_all(psi_edrlr_ref, y_tasks, max_N);
mean_corr_edr_ref = mean(corr_edr_ref, 1);
mean_corr_edrlr_ref = mean(corr_edrlr_ref, 1);

fg2 = figure('Name', 'Full accuracy curve at 5 m/s', 'Visible', 'off');
plot(1:max_N, mean_corr_static_edr, 'LineWidth', 1.8); hold on
plot(1:max_N, mean_corr_static_edrlr, 'LineWidth', 1.8);
plot(1:max_N, mean_corr_edr_ref, 'LineWidth', 1.8);
plot(1:max_N, mean_corr_edrlr_ref, 'LineWidth', 1.8);
grid on; xlabel('Number of modes (N)'); ylabel('Correlation with empirical activity');
legend({'Static EDR', 'Static EDR+LR', 'W_{delay} eigenbasis EDR', 'W_{delay} eigenbasis EDR+LR'}, 'Location', 'southeast');
title('Full reconstruction curve, speed=5 m/s (W_{delay} eigenvector basis)');
exportgraphics(fg2, fullfile(fig_out_dir, 'eigbasis_full_curve_5ms.png'));

%% Save
save(fullfile(results_out_dir, 'stage7_stage8_eigbasis_results.mat'), ...
    'speeds', 'acc20_eig_edr', 'acc20_eig_edrlr', 'mean_corr_static_edr', 'mean_corr_static_edrlr', ...
    'mean_corr_edr_ref', 'mean_corr_edrlr_ref', 'kappa', 'f_mu', '-v7.3');
fprintf('\nSaved figures -> %s\n', fig_out_dir);
fprintf('Saved results -> %s\n', fullfile(results_out_dir, 'stage7_stage8_eigbasis_results.mat'));
fprintf('=== W_delay eigenvector basis reconstruction complete ===\n');

%% Local functions (must appear at the end of a MATLAB script file)
function [mse_c, corr_c] = reconstruct_all(Psi, y_tasks, max_N)
    n_tasks = size(y_tasks, 2);
    mse_c = nan(n_tasks, max_N);
    corr_c = nan(n_tasks, max_N);
    for ti = 1:n_tasks
        y = y_tasks(:,ti);
        for N = 1:max_N
            y_hat = Psi(:,1:N) * (Psi(:,1:N)' * y);
            mse_c(ti,N) = mean((y - y_hat).^2);
            cc = corrcoef(y, y_hat);
            corr_c(ti,N) = cc(1,2);
        end
    end
end
