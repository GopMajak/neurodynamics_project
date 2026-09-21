%% Stage 1 + Stage 3 (parcel-level): EDR+LR construction + delayed-Kuramoto extension
%
% Combines Stage 1 and Stage 3 of the recovered
% "Original Project Code/brain_connectome_harmonics/connectome_harmonics.m"
% (your own earlier work on this thesis), skipping its Stage 2 (the
% expensive vertex-resolution eigendecomposition, which crashed there on
% memory limits and which we've since separately completed successfully
% via the server -- see Code/stage1_anatomical_graphs_v2/run_stage1_v2.m,
% whose EDR fit and exception count independently match this recovered
% script's Stage 2 log EXACTLY: A=0.0658/0.0659, lambda=0.1616/0.1617,
% 6,356,624 exceptions both times).
%
% This script:
%   Stage 1 (parcel-level, ~180 parcels): builds a quick EDR-only and
%     EDR+LR connectome directly at Glasser360 resolution (own EDR fit,
%     own long-range exceptions -- not aggregated from vertex level).
%   Stage 3: runs Budzinski et al. 2023's delayed Kuramoto model on both
%     graphs, single-run comparison + conduction-speed scan, using the
%     validated simulate_dKM_fast.m (confirmed numerically identical to
%     the reference simulate_dKM.m in Code/validation/validate_simulate_dKM_fast.m).
%
% All helper functions (calc_parcellate_matrix, delay_matrix_from_distance,
% simulate_dKM_fast, order_parameter) are reused directly from
% "Original Project Code/brain_connectome_harmonics/" via addpath -- not
% copied or reimplemented -- since this is the user's own prior code for
% this same thesis, not third-party material.

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
addpath(fullfile(recovered_dir, 'analysis'));

fig_out_dir = fullfile(project_root, 'Figures', 'novel_delay_extension');
results_out_dir = fullfile(project_root, 'Results', 'novel_delay_extension');
if ~exist(fig_out_dir, 'dir'); mkdir(fig_out_dir); end
if ~exist(results_out_dir, 'dir'); mkdir(results_out_dir); end

%% Stage 1: load data (Glasser360 parcellation, left hemisphere)
hemisphere = 'lh';
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';
parc_name = 'Glasser360';

[vertices, faces] = read_vtk(fullfile(vohryzek_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

cortex = dlmread(fullfile(vohryzek_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
num_vertices = length(cortex);

parc = dlmread(fullfile(vohryzek_dir, 'Data', 'parcellations', ...
    sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
parcels = unique(parc_cortex(parc_cortex>0));
num_parcels = length(parcels);
fprintf('Loaded %s parcellation: %d parcels\n', parc_name, num_parcels);

load(fullfile(pang_dir, 'data', 'empirical', ...
    'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');

%% Stage 1: parcel centroid distances and parcellated connectome
vertices_cortex = surface_midthickness.vertices(cortex_ind, :);
centroids = zeros(num_parcels, 3);
for p = 1:num_parcels
    centroids(p,:) = mean(vertices_cortex(parc_cortex==parcels(p), :), 1);
end
rr_parc = squareform(pdist(centroids));

connectome_parc = calc_parcellate_matrix(parc_cortex, avgSC_L);
C_parc = connectome_parc / max(connectome_parc(:));
clear avgSC_L

%% Stage 1: bin connectivity by distance and fit the exponential distance rule
NR = 60;
NSTD = 3;
DistRange = 40;

range_dist = max(rr_parc(:));
delta = range_dist / NR;
xcoor = delta/2 + delta*(0:NR-1);

index_parc = floor(rr_parc/delta) + 1;
index_parc(index_parc > NR) = NR;

sc_density = cell(1, NR);
sc_density_i = cell(1, NR);
sc_density_j = cell(1, NR);
ycoor2 = nan(1, NR);
for n = 1:NR
    [idx_i, idx_j] = find(index_parc == n);
    idx = find(index_parc == n);
    sc_density{n} = C_parc(idx);
    sc_density_i{n} = idx_i;
    sc_density_j{n} = idx_j;
    if ~isempty(idx)
        ycoor2(n) = mean(C_parc(idx));
    end
end

fit_start = find(xcoor >= 10, 1);
fit_ind = fit_start:NR;
fit_ind = fit_ind(~isnan(ycoor2(fit_ind)));

expfunc = @(A, x) (A(1)*exp(-A(2)*x));
sse = @(A) sum((expfunc(A, xcoor(fit_ind)) - ycoor2(fit_ind)).^2);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, 'Display', 'off');
A0 = [0.15, 0.18];
Afit = fminsearch(sse, A0, options);
lambda = Afit(2);
yl = Afit(1)*exp(-Afit(2)*xcoor);

fprintf('EDR fit (Glasser360): A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

%% Stage 1: detect long-range exceptions and build the EDR+LR connectome
Clong = zeros(num_parcels, num_parcels);
Clong_all = zeros(num_parcels, num_parcels);

for i = fit_start:NR
    if isempty(sc_density{i})
        continue
    end
    mv = mean(sc_density{i});
    st = std(sc_density{i});
    ind_exc = find(sc_density{i} > mv + NSTD*st);
    for n = 1:numel(ind_exc)
        ii = sc_density_i{i}(ind_exc(n));
        jj = sc_density_j{i}(ind_exc(n));
        Clong_all(ii,jj) = sc_density{i}(ind_exc(n));
        if rr_parc(ii,jj) > DistRange
            Clong(ii,jj) = sc_density{i}(ind_exc(n));
        end
    end
end

EDR_conn = Afit(1)*exp(-Afit(2)*rr_parc);
EDR_LRE_parc = EDR_conn;
EDR_LRE_parc(Clong>0) = Clong(Clong>0);

fprintf('Long-range exceptions found: %d parcel pairs (out of %d)\n', ...
    nnz(triu(Clong,1)), num_parcels*(num_parcels-1)/2);

fg1 = figure('Name', 'Stage 1 - EDR fit', 'Visible', 'off');
errorbar(xcoor, cellfun(@(x) mean(x,'omitnan'), sc_density), cellfun(@(x) std(x,'omitnan'), sc_density), 'o');
hold on
plot(xcoor, yl, 'r-', 'linewidth', 2)
xlabel('Distance (mm)'); ylabel('Connection strength (normalized)')
legend('binned SC (mean +/- SD)', 'EDR fit')
title(sprintf('Glasser360: lambda = %.4f /mm', lambda))
grid on
exportgraphics(fg1, fullfile(fig_out_dir, 'stage1_EDR_fit.png'));

fg2 = figure('Name', 'Stage 1 - Connectome comparison', 'Visible', 'off');
subplot(2,2,1); imagesc(C_parc); axis square; colorbar; title('Connectome (parcellated)')
subplot(2,2,2); imagesc(EDR_conn); axis square; colorbar; title('EDR fit')
subplot(2,2,3); imagesc(EDR_LRE_parc); axis square; colorbar; title('EDR+LR')
subplot(2,2,4); imagesc(Clong>0); axis square; colorbar; title(sprintf('Long-range exceptions (n=%d)', nnz(triu(Clong,1))))
colormap(flipud(bone))
exportgraphics(fg2, fullfile(fig_out_dir, 'stage1_connectome_comparison.png'));

%% Stage 3: time-delayed Kuramoto simulation, EDR vs EDR+LR connectomes
fprintf('\n=== STAGE 3: delayed Kuramoto on EDR vs EDR+LR (Glasser360 lh) ===\n');

N3 = size(EDR_conn, 1);
EDR_conn_km = EDR_conn;     EDR_conn_km(1:N3+1:end) = 0;
EDR_LRE_km  = EDR_LRE_parc; EDR_LRE_km(1:N3+1:end)  = 0;

%% Stage 3: single-run comparison (delayed vs non-delayed, EDR vs EDR+LR)
dt3 = 1e-3; T3 = 5.0; t3 = 0:dt3:T3;
f_mu3 = 10;
omega3 = f_mu3*2*pi*ones(N3,1);
kappa3 = 6;
speed_ref3 = 5;

rng(1); theta0_3 = 2*pi*rand(N3,1);
tau_ref3 = delay_matrix_from_distance(rr_parc, speed_ref3, dt3);

fprintf('Simulating single run (speed = %g m/s, T = %g s)...\n', speed_ref3, T3);
theta_edr_delay     = simulate_dKM_fast(EDR_conn_km, tau_ref3, omega3, kappa3, theta0_3, t3, dt3, 'euler');
theta_edrlr_delay   = simulate_dKM_fast(EDR_LRE_km,  tau_ref3, omega3, kappa3, theta0_3, t3, dt3, 'euler');
theta_edr_nodelay   = simulate_KM(EDR_conn_km, omega3, kappa3, theta0_3, t3, dt3, 'euler', 0);
theta_edrlr_nodelay = simulate_KM(EDR_LRE_km,  omega3, kappa3, theta0_3, t3, dt3, 'euler', 0);

R_edr_delay     = order_parameter(theta_edr_delay,     N3);
R_edrlr_delay   = order_parameter(theta_edrlr_delay,   N3);
R_edr_nodelay   = order_parameter(theta_edr_nodelay,   N3);
R_edrlr_nodelay = order_parameter(theta_edrlr_nodelay, N3);

n_tail3 = round(2/dt3);
fprintf('Mean R, last 2s (delayed):    EDR = %.4f   EDR+LR = %.4f\n', ...
    mean(R_edr_delay(end-n_tail3+1:end)), mean(R_edrlr_delay(end-n_tail3+1:end)));
fprintf('Mean R, last 2s (non-delayed): EDR = %.4f   EDR+LR = %.4f\n', ...
    mean(R_edr_nodelay(end-n_tail3+1:end)), mean(R_edrlr_nodelay(end-n_tail3+1:end)));

fg3 = figure('Name', 'Stage 3 - order parameter time series', 'Visible', 'off');
plot(t3, R_edr_delay, 'LineWidth', 1.5); hold on
plot(t3, R_edrlr_delay, 'LineWidth', 1.5);
plot(t3, R_edr_nodelay, '--', 'LineWidth', 1);
plot(t3, R_edrlr_nodelay, '--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Order parameter R(t)'); ylim([0 1]); grid on
legend('EDR (delayed)', 'EDR+LR (delayed)', 'EDR (no delay)', 'EDR+LR (no delay)', 'Location', 'best')
title(sprintf('Glasser360 lh: speed = %g m/s, \\kappa = %g', speed_ref3, kappa3))
exportgraphics(fg3, fullfile(fig_out_dir, 'stage3_order_parameter_timeseries.png'));

save(fullfile(results_out_dir, 'stage3_dKM_single_run.mat'), ...
    't3', 'R_edr_delay', 'R_edrlr_delay', 'R_edr_nodelay', 'R_edrlr_nodelay', ...
    'speed_ref3', 'kappa3', 'f_mu3', '-v7.3');

%% Stage 3: conduction-speed scan, EDR vs EDR+LR
n_speeds3 = 20;
speeds3 = linspace(1, 30, n_speeds3);
T_scan3 = 3.0; t_scan3 = 0:dt3:T_scan3;
n_half3 = floor(length(t_scan3)/2);

rng(42); theta0_scan3 = 2*pi*rand(N3,1);

R_edr_scan = nan(1, n_speeds3);
R_edrlr_scan = nan(1, n_speeds3);

fprintf('\nConduction speed scan (%d speeds, EDR vs EDR+LR)...\n', n_speeds3);
for si = 1:n_speeds3
    tau_s3 = delay_matrix_from_distance(rr_parc, speeds3(si), dt3);

    th_edr   = simulate_dKM_fast(EDR_conn_km, tau_s3, omega3, kappa3, theta0_scan3, t_scan3, dt3, 'euler');
    th_edrlr = simulate_dKM_fast(EDR_LRE_km,  tau_s3, omega3, kappa3, theta0_scan3, t_scan3, dt3, 'euler');

    R_edr_scan(si)   = mean(order_parameter(th_edr(n_half3:end,:),   N3));
    R_edrlr_scan(si) = mean(order_parameter(th_edrlr(n_half3:end,:), N3));

    fprintf('  speed = %5.2f m/s | R_EDR = %.3f | R_EDR+LR = %.3f\n', ...
        speeds3(si), R_edr_scan(si), R_edrlr_scan(si));
end

fg4 = figure('Name', 'Stage 3 - conduction speed scan', 'Visible', 'off');
plot(speeds3, R_edr_scan, 'o-', 'LineWidth', 1.5); hold on
plot(speeds3, R_edrlr_scan, 's-', 'LineWidth', 1.5);
xlabel('Conduction speed (m/s)'); ylabel('Steady-state order parameter R'); ylim([0 1]); grid on
legend('EDR', 'EDR+LR', 'Location', 'best')
title('Glasser360 lh: synchrony vs. conduction speed, EDR vs EDR+LR')
exportgraphics(fg4, fullfile(fig_out_dir, 'stage3_speed_scan.png'));

save(fullfile(results_out_dir, 'stage3_dKM_speed_scan.mat'), ...
    'speeds3', 'R_edr_scan', 'R_edrlr_scan', 'kappa3', 'f_mu3', 'T_scan3', '-v7.3');
fprintf('\nSaved figures -> %s\n', fig_out_dir);
fprintf('Saved results -> %s\n', results_out_dir);
fprintf('=== Stage 1 + Stage 3 (parcel-level) complete ===\n');
