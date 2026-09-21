%% Robustness check for the delay-vs-no-delay simulated FC result
%
% simulate_FC_delay_vs_nodelay.m found delay beating a no-delay control on
% both EDR-only and EDR+LR, peaking sharply around speed=3 m/s -- the
% first positive result for delay in the whole investigation, but based
% on a coarse 12-point speed grid and a single random seed (rng(7)).
%
% This checks whether that peak survives (a) finer speed resolution near
% where it was found, and (b) repetition across multiple random initial
% phase conditions -- before treating it as a robust finding rather than
% a lucky single-seed / coarse-grid artifact.
%
% Reuses the same graph-construction and simulation code as
% simulate_FC_delay_vs_nodelay.m.

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
addpath(fullfile(project_root, 'Code', 'utils'));

results_out_dir = fullfile(project_root, 'Results', 'novel_delay_extension');
fig_out_dir = fullfile(project_root, 'Figures', 'novel_delay_extension');

%% Load group-average empirical resting-state FC
rest_dir = fullfile(project_root, 'data', 'processed', 'rfMRI_parcellated');
d = load(fullfile(rest_dir, 'FC_all_subjects.mat'), 'FC_all', 'valid_subj');
FC_all = d.FC_all(:,:,d.valid_subj);
n_valid = size(FC_all, 3);
num_parcels = size(FC_all, 1);
Z = atanh(FC_all); Z(~isfinite(Z)) = NaN;
FC_empirical_group = tanh(mean(Z, 3, 'omitnan'));
FC_empirical_group(1:num_parcels+1:end) = 1;
fprintf('Loaded group-average empirical FC (%d subjects, %d parcels).\n', n_valid, num_parcels);

%% Rebuild EDR-only / EDR+LR parcel graphs (identical recipe to simulate_FC_delay_vs_nodelay.m)
hemisphere = 'lh'; surface_interest = 'fsLR_32k'; mesh_interest = 'midthickness'; parc_name = 'Glasser360';
[vertices, faces] = read_vtk(fullfile(vohryzek_dir, 'Data', 'template_surfaces', sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices'; surface_midthickness.faces = faces';
cortex = dlmread(fullfile(vohryzek_dir, 'Data', 'template_surfaces', sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere))); %#ok<DLMRD>
cortex_ind = find(cortex);
parc = dlmread(fullfile(vohryzek_dir, 'Data', 'parcellations', sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere))); %#ok<DLMRD>
parc_cortex = parc(cortex_ind);
parcels = unique(parc_cortex(parc_cortex>0));

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

N = num_parcels;
EDR_conn_km = EDR_conn;    EDR_conn_km(1:N+1:end) = 0;
EDR_LRE_km  = EDR_LRE_parc; EDR_LRE_km(1:N+1:end)  = 0;

graph_names = {'EDR-only', 'EDR+LR'};
graph_W = {EDR_conn_km, EDR_LRE_km};

%% Simulation parameters
dt = 1e-3;
T_total = 20.0; t = 0:dt:T_total;
transient_s = 5.0; n_transient = round(transient_s/dt);
f_mu = 10; omega = f_mu*2*pi*ones(N,1);
kappa = 6;

speeds_fine = [0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 7 8]; % finer grid around/beyond the speed=3 peak
seeds = [7 11 23 42 99]; % rng(7) matches the original run; 4 new seeds added for robustness

fc_upper_mask = triu(true(N), 1);
fc_empirical_vec = FC_empirical_group(fc_upper_mask);

n_graphs = numel(graph_names); n_speeds = numel(speeds_fine); n_seeds = numel(seeds);
r_delay = nan(n_graphs, n_speeds, n_seeds);
r_nodelay = nan(n_graphs, n_seeds);

for gi = 1:n_graphs
    W = graph_W{gi};
    fprintf('\n=== Graph: %s ===\n', graph_names{gi});

    for sdi = 1:n_seeds
        rng(seeds(sdi)); theta0 = 2*pi*rand(N,1);

        theta_nodelay = simulate_KM(W, omega, kappa, theta0, t, dt, 'euler', 0);
        sig_nodelay = cos(theta_nodelay(n_transient+1:end, :));
        FC_nodelay = corr(sig_nodelay);
        r_nodelay(gi,sdi) = corr(FC_nodelay(fc_upper_mask), fc_empirical_vec);
        fprintf('  seed=%3d | no-delay r = %.4f\n', seeds(sdi), r_nodelay(gi,sdi));

        for spi = 1:n_speeds
            tau = delay_matrix_from_distance(rr_parc, speeds_fine(spi), dt);
            theta_delay = simulate_dKM_fast(W, tau, omega, kappa, theta0, t, dt, 'euler');
            sig_delay = cos(theta_delay(n_transient+1:end, :));
            FC_delay = corr(sig_delay);
            r_delay(gi,spi,sdi) = corr(FC_delay(fc_upper_mask), fc_empirical_vec);
        end
        [best_r, best_idx] = max(r_delay(gi,:,sdi));
        fprintf('  seed=%3d | best delayed: speed=%.1f m/s, r=%.4f (vs no-delay %.4f)\n', ...
            seeds(sdi), speeds_fine(best_idx), best_r, r_nodelay(gi,sdi));
    end
end

%% Summary across seeds
fprintf('\n=== Summary: mean +/- SD across %d seeds ===\n', n_seeds);
for gi = 1:n_graphs
    fprintf('\n-- %s --\n', graph_names{gi});
    fprintf('  no-delay:  mean r = %.4f +/- %.4f\n', mean(r_nodelay(gi,:)), std(r_nodelay(gi,:)));
    for spi = 1:n_speeds
        m = mean(r_delay(gi,spi,:)); s = std(r_delay(gi,spi,:));
        fprintf('  speed=%4.1f m/s: mean r = %.4f +/- %.4f\n', speeds_fine(spi), m, s);
    end
end

save(fullfile(results_out_dir, 'simulate_FC_delay_robustness_results.mat'), ...
    'r_delay', 'r_nodelay', 'speeds_fine', 'seeds', 'graph_names', 'kappa', 'f_mu', 'T_total', 'transient_s', 'n_valid', '-v7.3');

%% Figure: mean +/- SD across seeds, per graph
fg = figure('Name', 'Delay-vs-no-delay robustness (multi-seed)', 'Visible', 'off');
colors = lines(n_graphs);
hold on
for gi = 1:n_graphs
    m = squeeze(mean(r_delay(gi,:,:), 3));
    s = squeeze(std(r_delay(gi,:,:), 0, 3));
    errorbar(speeds_fine, m, s, 'o-', 'LineWidth', 1.5, 'Color', colors(gi,:), 'DisplayName', sprintf('%s (delayed)', graph_names{gi}));
    nd_m = mean(r_nodelay(gi,:));
    plot([min(speeds_fine) max(speeds_fine)], [nd_m nd_m], '--', 'Color', colors(gi,:), 'DisplayName', sprintf('%s (no delay)', graph_names{gi}));
end
xlabel('Conduction speed (m/s)'); ylabel('Correlation with empirical group-average FC (mean +/- SD, %d seeds)');
legend('Location', 'best'); grid on
title(sprintf('Robustness check: %d seeds, finer speed grid near the speed=3 m/s peak', n_seeds));
exportgraphics(fg, fullfile(fig_out_dir, 'simulate_FC_delay_robustness.png'));

fprintf('\nSaved results -> %s\n', results_out_dir);
fprintf('Saved figure -> %s\n', fig_out_dir);
fprintf('=== Robustness check complete ===\n');
