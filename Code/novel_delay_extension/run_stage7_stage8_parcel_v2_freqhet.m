%% Stage 7 + 8 v2 (parcel-level): with frequency heterogeneity
%
% Identical to run_stage7_stage8_parcel.m EXCEPT for one change: natural
% frequencies are drawn from a Gaussian spread (sigma_f = 0.5 Hz) around
% the shared 10 Hz mean, instead of being identical for every node.
%
% Rationale (see sweep_coupling_strength.m and sliding_window_complexity.m
% for the ruled-out alternatives): with identical frequencies,
% cos(theta(t)) is mathematically confined to ~2 effective dimensions
% (participation ratio never exceeded 2.3 across an entire coupling sweep
% AND an entire simulation timecourse), regardless of coupling strength or
% analysis window -- a structural consequence of every oscillator sharing
% one carrier frequency, not a tuning problem. sweep_frequency_heterogeneity.m
% found sigma_f=0.5 Hz escapes this (participation ratio ~2.9-3.5) while
% still preserving a clear EDR-vs-EDR+LR difference in synchrony (R=0.26
% vs 0.54) -- larger spreads wash out that difference entirely (both
% graphs converge to the same near-incoherent floor by sigma_f=1 Hz),
% which would mean the resulting basis reflects frequency noise rather
% than the connectome structure being studied.
%
% This file is kept alongside the original (identical-frequency) version
% rather than overwriting it, so both operating points remain available
% for comparison.

clear; clc;

%% Paths
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir));

pang_dir = fullfile(project_root, 'Reference Code Material', 'pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
vohryzek_dir = fullfile(project_root, 'Reference Code Material', 'Vohryzeketal2024', 'vohryzek2024_EDRLR');
recovered_dir = fullfile(project_root, 'Original Project Code', 'brain_connectome_harmonics');

addpath(genpath(fullfile(pang_dir, 'functions_matlab')));
addpath(fullfile(recovered_dir, 'functions'));
addpath(fullfile(recovered_dir, 'simulation'));
addpath(fullfile(project_root, 'Code', 'utils'));

fig_out_dir = fullfile(project_root, 'Figures', 'novel_delay_extension');
results_out_dir = fullfile(project_root, 'Results', 'novel_delay_extension');
if ~exist(fig_out_dir, 'dir'); mkdir(fig_out_dir); end
if ~exist(results_out_dir, 'dir'); mkdir(results_out_dir); end

%% Rebuild Stage 1's parcel-level EDR-only / EDR+LR graphs
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
fprintf('EDR fit (Glasser360): A = %.4f, lambda = %.4f /mm\n', Afit(1), Afit(2));
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
fprintf('Long-range exceptions: %d parcel pairs\n', nnz(triu(Clong,1)));

N3 = num_parcels;
EDR_conn_km = EDR_conn;     EDR_conn_km(1:N3+1:end) = 0;
EDR_LRE_km  = EDR_LRE_parc; EDR_LRE_km(1:N3+1:end)  = 0;

%% Stage 4 (static baseline, parcel-level): normalized Laplacian eigenmodes
fprintf('\nComputing static (Laplacian) eigenmodes at parcel resolution...\n');
[psi_static_edr, ~]   = normalized_laplacian_eigenmodes(EDR_conn_km, num_parcels);
[psi_static_edrlr, ~] = normalized_laplacian_eigenmodes(EDR_LRE_km, num_parcels);

%% Stage 7: delay-informed basis via SVD of simulated cos(theta(t))
fprintf('\n=== STAGE 7 (v2, frequency heterogeneity): delay-informed basis extraction ===\n');
dt = 1e-3;
T_transient = 5.0;
T_steady = 5.0;
t_full = 0:dt:(T_transient + T_steady);
n_transient_steps = round(T_transient/dt);

f_mu = 10;
sigma_f = 0.5; % Hz -- validated operating point (sweep_frequency_heterogeneity.m)
kappa = 6;
speed_basis = 5;
rng(2); omega = 2*pi*(f_mu + sigma_f*randn(N3,1)); % matches the sweep's frequency draw
rng(1); theta0 = 2*pi*rand(N3,1);
tau = delay_matrix_from_distance(rr_parc, speed_basis, dt);

fprintf('Simulating EDR-only (speed=%g m/s, sigma_f=%g Hz, T=%g s) for basis extraction...\n', speed_basis, sigma_f, t_full(end));
theta_edr = simulate_dKM_fast(EDR_conn_km, tau, omega, kappa, theta0, t_full, dt, 'euler');
fprintf('Simulating EDR+LR (speed=%g m/s, sigma_f=%g Hz, T=%g s) for basis extraction...\n', speed_basis, sigma_f, t_full(end));
theta_edrlr = simulate_dKM_fast(EDR_LRE_km, tau, omega, kappa, theta0, t_full, dt, 'euler');

X_edr = cos(theta_edr(n_transient_steps+1:end, :))';
X_edrlr = cos(theta_edrlr(n_transient_steps+1:end, :))';

[U_edr, S_edr, ~] = svd(X_edr, 'econ');
[U_edrlr, S_edrlr, ~] = svd(X_edrlr, 'econ');

sv_edr = diag(S_edr); sv_edrlr = diag(S_edrlr);
peff_edr = (sum(sv_edr.^2))^2 / sum(sv_edr.^4);
peff_edrlr = (sum(sv_edrlr.^2))^2 / sum(sv_edrlr.^4);
fprintf('Participation ratio -- EDR: %.2f, EDR+LR: %.2f\n', peff_edr, peff_edrlr);
fprintf('EDR-only:  top-5 singular values (%% variance): %s\n', ...
    mat2str(round(100*sv_edr(1:5).^2/sum(sv_edr.^2), 1)'));
fprintf('EDR+LR:    top-5 singular values (%% variance): %s\n', ...
    mat2str(round(100*sv_edrlr(1:5).^2/sum(sv_edrlr.^2), 1)'));

psi_delay_edr = U_edr;
psi_delay_edrlr = U_edrlr;

%% Stage 8: reconstruct group-averaged task-activation maps
fprintf('\n=== STAGE 8: task-activation reconstruction ===\n');
data = load(fullfile(vohryzek_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
fieldNames = fieldnames(data.zstat);
representTask = [3, 11, 19, 30, 41, 44, 47];

conditions = struct( ...
    'name', {'static_EDR', 'static_EDRLR', 'delay_EDR', 'delay_EDRLR'}, ...
    'basis', {psi_static_edr, psi_static_edrlr, psi_delay_edr, psi_delay_edrlr});

max_N = num_parcels;
mse_curves = nan(numel(conditions), numel(representTask), max_N);
corr_curves = nan(numel(conditions), numel(representTask), max_N);

for ti = 1:numel(representTask)
    tk = representTask(ti);
    fieldName = fieldNames{tk};
    activation_map = nanmean(data.zstat.(fieldName), 2);
    y = parcellate_average(activation_map(cortex_ind), parc_cortex);

    for ci = 1:numel(conditions)
        Psi = conditions(ci).basis;
        for N = 1:max_N
            y_hat = Psi(:,1:N) * (Psi(:,1:N)' * y);
            mse_curves(ci, ti, N) = mean((y - y_hat).^2);
            cc = corrcoef(y, y_hat);
            corr_curves(ci, ti, N) = cc(1,2);
        end
    end
    fprintf('  task %d/%d (%s) done\n', ti, numel(representTask), fieldName);
end

%% Figures
mean_mse = squeeze(mean(mse_curves, 2));
mean_corr = squeeze(mean(corr_curves, 2));

fg1 = figure('Name', 'Stage 8 v2 - reconstruction accuracy', 'Visible', 'off');
subplot(1,2,1)
plot(1:max_N, mean_corr', 'LineWidth', 1.8); grid on
xlabel('Number of modes (N)'); ylabel('Correlation with empirical activity')
legend({'Static EDR', 'Static EDR+LR', 'Delay-informed EDR', 'Delay-informed EDR+LR'}, 'Location', 'southeast')
title('Reconstruction accuracy (correlation)')
subplot(1,2,2)
plot(1:max_N, mean_mse', 'LineWidth', 1.8); grid on
xlabel('Number of modes (N)'); ylabel('MSE vs. empirical activity')
legend({'Static EDR', 'Static EDR+LR', 'Delay-informed EDR', 'Delay-informed EDR+LR'}, 'Location', 'northeast')
title('Reconstruction error (MSE)')
sgtitle(sprintf('Glasser360, 7-task average (sigma_f=%g Hz, freq. heterogeneity)', sigma_f))
exportgraphics(fg1, fullfile(fig_out_dir, 'stage8_v2_reconstruction_accuracy.png'));

%% Zoomed view (first 20 modes -- the proposal's primary comparison range)
fg2 = figure('Name', 'Stage 8 v2 - zoomed (N<=20)', 'Visible', 'off');
subplot(1,2,1)
plot(1:20, mean_corr(:,1:20)', 'LineWidth', 1.8); grid on
xlabel('Number of modes (N)'); ylabel('Correlation with empirical activity')
legend({'Static EDR', 'Static EDR+LR', 'Delay-informed EDR', 'Delay-informed EDR+LR'}, 'Location', 'southeast')
title('Reconstruction accuracy (N<=20)')
subplot(1,2,2)
plot(1:20, mean_mse(:,1:20)', 'LineWidth', 1.8); grid on
xlabel('Number of modes (N)'); ylabel('MSE vs. empirical activity')
legend({'Static EDR', 'Static EDR+LR', 'Delay-informed EDR', 'Delay-informed EDR+LR'}, 'Location', 'northeast')
title('Reconstruction error (N<=20)')
exportgraphics(fg2, fullfile(fig_out_dir, 'stage8_v2_zoomed_N20.png'));

%% Save
save(fullfile(results_out_dir, 'stage7_stage8_v2_results.mat'), ...
    'mse_curves', 'corr_curves', 'representTask', 'fieldNames', 'conditions', ...
    'speed_basis', 'kappa', 'f_mu', 'sigma_f', 'sv_edr', 'sv_edrlr', 'peff_edr', 'peff_edrlr', '-v7.3');
fprintf('\nSaved figures -> %s\n', fig_out_dir);
fprintf('Saved results -> %s\n', fullfile(results_out_dir, 'stage7_stage8_v2_results.mat'));
fprintf('=== Stage 7 + 8 v2 (parcel-level, freq. heterogeneity) complete ===\n');
