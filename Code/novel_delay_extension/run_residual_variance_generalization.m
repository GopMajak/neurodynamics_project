%% Does the residual-variance null result generalize across graph type and parameter point?
%
% run_residual_variance_test.m found delay-informed modes explain LESS of
% the static-N=20 residual than even random noise, for EDR+LR at the
% standard reference point (f_mu=10Hz, speed=5m/s). This checks whether
% that holds for EDR-only too, and whether the (f_mu=80Hz, speed=8.63m/s)
% point -- the single combination that showed a (small, unreplicated)
% positive accuracy margin in sweep_eigbasis_fmu.m -- behaves any
% differently. If delay's residual-explaining power is ever going to beat
% random noise, this is the most likely place to see it.

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

results_out_dir = fullfile(project_root, 'Results', 'novel_delay_extension');

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
EDR_conn_km = EDR_conn;    EDR_conn_km(1:N3+1:end) = 0;
EDR_LRE_km  = EDR_LRE_parc; EDR_LRE_km(1:N3+1:end)  = 0;

graphs = struct('name', {'EDR-only', 'EDR+LR'}, 'W', {EDR_conn_km, EDR_LRE_km});

%% Load per-subject task data
fprintf('Loading task data...\n');
data = load(fullfile(vohryzek_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
fieldNames = fieldnames(data.zstat);
representTask = [3, 11, 19, 30, 41, 44, 47];
n_tasks = numel(representTask);
N_base = 20;
K_values = [5, 10, 20];

Y_tasks = cell(1, n_tasks);
for ti = 1:n_tasks
    Y_tasks{ti} = parcellate_average(data.zstat.(fieldNames{representTask(ti)})(cortex_ind, :), parc_cortex);
end
n_subj = size(Y_tasks{1}, 2);

%% Random control bases (shared across all conditions for consistency)
n_random_draws = 5;
psi_random = cell(1, n_random_draws);
for d = 1:n_random_draws
    rng(1000 + d);
    [Q, ~] = qr(randn(num_parcels, num_parcels), 0);
    psi_random{d} = Q;
end

%% Parameter points to test
dt = 1e-3;
param_points = struct('f_mu', {10, 80}, 'speed', {5, 8.63}, 'label', {'standard (f_mu=10, speed=5)', 'best-margin point (f_mu=80, speed=8.63)'});

fprintf('\n=== Residual variance: does it generalize across graph and parameter point? ===\n');
summary_rows = {};

for gi = 1:numel(graphs)
    W = graphs(gi).W;
    [psi_static, ~] = normalized_laplacian_eigenmodes(W, num_parcels);
    Psi_base = psi_static(:,1:N_base);

    Residual_tasks = cell(1, n_tasks);
    for ti = 1:n_tasks
        Y = Y_tasks{ti};
        Y_hat_base = Psi_base * (Psi_base' * Y);
        Residual_tasks{ti} = Y - Y_hat_base;
    end

    for pi = 1:numel(param_points)
        f_mu = param_points(pi).f_mu;
        omega = f_mu*2*pi*ones(N3,1);
        speed = param_points(pi).speed;
        tau = delay_matrix_from_distance(rr_parc, speed, dt);
        [psi_delay, ~] = build_delay_eigenvector_basis(W, tau, omega, 6, dt);

        fprintf('\n-- %s, %s --\n', graphs(gi).name, param_points(pi).label);
        for K = K_values
            r2_static_ext_all = nan(n_subj, n_tasks);
            r2_delay_all = nan(n_subj, n_tasks);
            r2_random_all = nan(n_subj, n_tasks);

            for ti = 1:n_tasks
                R = Residual_tasks{ti};
                var_R = sum(R.^2, 1);

                Psi_ext = psi_static(:, N_base+1:N_base+K);
                R_hat = Psi_ext * (Psi_ext' * R);
                r2_static_ext_all(:,ti) = 1 - sum((R - R_hat).^2, 1)' ./ var_R';

                Psi_d = psi_delay(:, 1:K);
                R_hat = Psi_d * (Psi_d' * R);
                r2_delay_all(:,ti) = 1 - sum((R - R_hat).^2, 1)' ./ var_R';

                r2_rand_draws = nan(n_subj, n_random_draws);
                for d = 1:n_random_draws
                    Psi_r = psi_random{d}(:, 1:K);
                    R_hat = Psi_r * (Psi_r' * R);
                    r2_rand_draws(:,d) = 1 - sum((R - R_hat).^2, 1)' ./ var_R';
                end
                r2_random_all(:,ti) = mean(r2_rand_draws, 2);
            end

            m_ext = mean(mean(r2_static_ext_all, 2), 'omitnan');
            m_delay = mean(mean(r2_delay_all, 2), 'omitnan');
            m_rand = mean(mean(r2_random_all, 2), 'omitnan');
            fprintf('  K=%2d: static_ext=%.4f, delay=%.4f, random=%.4f  %s\n', K, m_ext, m_delay, m_rand, ...
                merge_flag(m_delay > m_rand));
            summary_rows(end+1,:) = {graphs(gi).name, param_points(pi).label, K, m_ext, m_delay, m_rand}; %#ok<AGROW>
        end
    end
end

save(fullfile(results_out_dir, 'residual_variance_generalization_results.mat'), 'summary_rows');
fprintf('\nSaved -> %s\n', fullfile(results_out_dir, 'residual_variance_generalization_results.mat'));
fprintf('=== Generalization check complete ===\n');

function s = merge_flag(cond)
    if cond, s = '<-- delay beats random here'; else, s = ''; end
end
