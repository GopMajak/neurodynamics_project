%% Simulated FC (delay vs no-delay) vs empirical group-average resting-state FC
%
% Pivots from testing delay as a LINEAR RECONSTRUCTION BASIS (robustly
% falsified across six separate tests: hybrid basis, residual variance,
% generalization check, per-subject stats, fmu/kappa sweeps) to testing it
% as a GENERATIVE DYNAMICAL FEATURE instead: simulate the delayed Kuramoto
% model (dKM) directly on the EDR-only / EDR+LR structural graphs, and see
% whether the resulting simulated FC matches empirical resting-state FC
% better than a no-delay Kuramoto (KM) control on the SAME graph. This
% mirrors the logic of Cabral et al. 2011 / Budzinski et al. 2023: with
% identical oscillator frequencies, conduction delay + network topology
% alone can generate realistic spatial FC structure that a delay-free
% model on the same graph cannot.
%
% Target: group-averaged empirical FC across all successfully extracted
% HCP subjects (Fisher z-transform average, not per-subject comparison --
% matches the group-average structural connectome already used throughout
% this pipeline). See extract_parcellate_rfMRI_all_subjects.m.
%
% All simulation functions (simulate_dKM_fast, simulate_KM,
% delay_matrix_from_distance) are reused directly from
% "Original Project Code/brain_connectome_harmonics/" (validated earlier
% in this project), not reimplemented.

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

fig_out_dir = fullfile(project_root, 'Figures', 'novel_delay_extension');
results_out_dir = fullfile(project_root, 'Results', 'novel_delay_extension');
if ~exist(fig_out_dir, 'dir'); mkdir(fig_out_dir); end
if ~exist(results_out_dir, 'dir'); mkdir(results_out_dir); end

%% Load group-average empirical resting-state FC
rest_dir = fullfile(project_root, 'data', 'processed', 'rfMRI_parcellated');
fc_file = fullfile(rest_dir, 'FC_all_subjects.mat');
if ~exist(fc_file, 'file')
    error('FC_all_subjects.mat not found (%s) -- run extract_parcellate_rfMRI_all_subjects.m on the server first.', fc_file);
end
d = load(fc_file, 'FC_all', 'valid_subj');
FC_all = d.FC_all(:,:,d.valid_subj);
n_valid = size(FC_all, 3);
fprintf('Loaded empirical FC for %d valid subjects.\n', n_valid);
if n_valid < 2
    error('Fewer than 2 valid subjects in FC_all_subjects.mat -- extraction likely did not finish successfully.');
end

num_parcels = size(FC_all, 1);
Z = atanh(FC_all);
Z(~isfinite(Z)) = NaN; % guard diagonal (atanh(1) = Inf)
Z_mean = mean(Z, 3, 'omitnan');
FC_empirical_group = tanh(Z_mean);
FC_empirical_group(1:num_parcels+1:end) = 1; % restore exact diagonal
fprintf('Group-average empirical FC built (Fisher z-average, %d subjects, %d parcels).\n', n_valid, num_parcels);

%% Rebuild EDR-only / EDR+LR parcel graphs (same fitting recipe as run_stage1_stage3_parcel.m)
hemisphere = 'lh'; surface_interest = 'fsLR_32k'; mesh_interest = 'midthickness'; parc_name = 'Glasser360';
[vertices, faces] = read_vtk(fullfile(vohryzek_dir, 'Data', 'template_surfaces', sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices'; surface_midthickness.faces = faces';
cortex = dlmread(fullfile(vohryzek_dir, 'Data', 'template_surfaces', sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere))); %#ok<DLMRD>
cortex_ind = find(cortex);
parc = dlmread(fullfile(vohryzek_dir, 'Data', 'parcellations', sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere))); %#ok<DLMRD>
parc_cortex = parc(cortex_ind);
parcels = unique(parc_cortex(parc_cortex>0));
assert(numel(parcels) == num_parcels, 'Parcel count mismatch between empirical FC (%d) and structural parcellation (%d).', num_parcels, numel(parcels));

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

fprintf('EDR fit: A=%.4f, lambda=%.4f, LR exceptions=%d/%d parcel pairs\n', Afit(1), Afit(2), nnz(triu(Clong,1)), N*(N-1)/2);

graph_names = {'EDR-only', 'EDR+LR'};
graph_W = {EDR_conn_km, EDR_LRE_km};

%% Simulation parameters
dt = 1e-3;
T_total = 20.0;          % seconds of simulated model time
t = 0:dt:T_total;
transient_s = 5.0;       % discard as transient before computing FC
n_transient = round(transient_s/dt);
f_mu = 10;                % Hz, IDENTICAL across nodes -- delay + topology alone should generate FC structure
omega = f_mu*2*pi*ones(N,1);
kappa = 6;
speeds = [1 2 3 4 5 6 8 10 15 20 25 30]; % m/s, conduction speed scan

rng(7); theta0 = 2*pi*rand(N,1);

fc_upper_mask = triu(true(N), 1);
fc_empirical_vec = FC_empirical_group(fc_upper_mask);

results = repmat(struct('name', '', 'r_nodelay', nan, 'speeds', speeds, 'r_delay', nan(1,numel(speeds)), 'FC_nodelay', []), 1, numel(graph_names));

for gi = 1:numel(graph_names)
    W = graph_W{gi};
    gname = graph_names{gi};
    fprintf('\n=== Graph: %s ===\n', gname);

    fprintf('  Simulating no-delay KM control...\n');
    tic;
    theta_nodelay = simulate_KM(W, omega, kappa, theta0, t, dt, 'euler', 0);
    fprintf('  done (%.1f s)\n', toc);
    sig_nodelay = cos(theta_nodelay(n_transient+1:end, :));
    FC_nodelay = corr(sig_nodelay);
    r_nodelay = corr(FC_nodelay(fc_upper_mask), fc_empirical_vec);
    fprintf('  No-delay KM: FC-empirical correlation = %.4f\n', r_nodelay);

    r_delay = nan(1, numel(speeds));
    for si = 1:numel(speeds)
        tau = delay_matrix_from_distance(rr_parc, speeds(si), dt);
        theta_delay = simulate_dKM_fast(W, tau, omega, kappa, theta0, t, dt, 'euler');
        sig_delay = cos(theta_delay(n_transient+1:end, :));
        FC_delay = corr(sig_delay);
        r_delay(si) = corr(FC_delay(fc_upper_mask), fc_empirical_vec);
        fprintf('  speed = %5.2f m/s | FC-empirical correlation = %.4f\n', speeds(si), r_delay(si));
    end

    [best_r, best_idx] = max(r_delay);
    verdict = 'no-delay wins';
    if best_r > r_nodelay; verdict = 'DELAY WINS'; end
    fprintf('  Best delayed speed: %.2f m/s, r = %.4f (vs no-delay r = %.4f) -- %s\n', ...
        speeds(best_idx), best_r, r_nodelay, verdict);

    results(gi).name = gname;
    results(gi).r_nodelay = r_nodelay;
    results(gi).r_delay = r_delay;
    results(gi).FC_nodelay = FC_nodelay;
end

save(fullfile(results_out_dir, 'simulate_FC_delay_vs_nodelay_results.mat'), ...
    'results', 'speeds', 'kappa', 'f_mu', 'T_total', 'transient_s', 'FC_empirical_group', 'n_valid', '-v7.3');

%% Figure
fg = figure('Name', 'Simulated vs empirical FC: delay vs no-delay', 'Visible', 'off');
hold on
colors = lines(numel(graph_names));
for gi = 1:numel(graph_names)
    plot(results(gi).speeds, results(gi).r_delay, 'o-', 'LineWidth', 1.5, 'Color', colors(gi,:), ...
        'DisplayName', sprintf('%s (delayed)', results(gi).name));
    yl = results(gi).r_nodelay;
    plot([min(speeds) max(speeds)], [yl yl], '--', 'Color', colors(gi,:), ...
        'DisplayName', sprintf('%s (no delay)', results(gi).name));
end
xlabel('Conduction speed (m/s)'); ylabel('Correlation with empirical group-average FC');
legend('Location', 'best'); grid on
title(sprintf('Simulated FC accuracy: delayed dKM vs no-delay KM (n=%d subjects)', n_valid));
exportgraphics(fg, fullfile(fig_out_dir, 'simulate_FC_delay_vs_nodelay.png'));

fprintf('\nSaved results -> %s\n', results_out_dir);
fprintf('Saved figure -> %s\n', fig_out_dir);
fprintf('=== Simulated FC (delay vs no-delay) comparison complete ===\n');
