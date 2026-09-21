%% Residual variance test (Option B): does delay explain what static EDR+LR misses?
%
% Fixes the static EDR+LR reconstruction at N=20 (the proposal's primary
% comparison point) and takes its per-subject, per-task RESIDUAL (leftover
% error). Then asks: how much of that residual can be explained by K
% additional dimensions, comparing three sources for those K dimensions:
%   1. MORE STATIC modes (21, 22, ..., 20+K) -- the natural "just keep
%      extending the same model" control.
%   2. DELAY-INFORMED modes (1..K of the W_delay eigenvector basis) -- the
%      condition of interest.
%   3. RANDOM orthonormal directions (averaged over several draws) -- a
%      null-model control for "any generic K-dimensional subspace will
%      explain some residual variance just by dimension-counting".
%
% Unlike the hybrid test (Option A), this does NOT trade away any of the
% 20 static modes -- delay dimensions are purely additive on top of the
% fixed static-20 reconstruction, so it's a genuinely different question:
% not "should delay replace some static capacity" but "does delay carry
% information the static model's residual doesn't already contain".

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

%% Random orthonormal control bases (several draws, averaged)
n_random_draws = 5;
psi_random = cell(1, n_random_draws);
for d = 1:n_random_draws
    rng(1000 + d);
    [Q, ~] = qr(randn(num_parcels, num_parcels), 0);
    psi_random{d} = Q;
end

%% Load per-subject task data
fprintf('Loading task data...\n');
data = load(fullfile(vohryzek_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
fieldNames = fieldnames(data.zstat);
representTask = [3, 11, 19, 30, 41, 44, 47];
n_tasks = numel(representTask);
N_base = 20; % fixed static reconstruction budget
K_values = 1:40; % additional dimensions tested
n_K = numel(K_values);

Y_tasks = cell(1, n_tasks);
for ti = 1:n_tasks
    Y_tasks{ti} = parcellate_average(data.zstat.(fieldNames{representTask(ti)})(cortex_ind, :), parc_cortex);
end
n_subj = size(Y_tasks{1}, 2);

%% Compute residuals (fixed, per subject per task) after static-N=20 reconstruction
fprintf('Computing static-N=%d residuals...\n', N_base);
Psi_base = psi_static(:,1:N_base);
Residual_tasks = cell(1, n_tasks);
for ti = 1:n_tasks
    Y = Y_tasks{ti};
    Y_hat_base = Psi_base * (Psi_base' * Y);
    Residual_tasks{ti} = Y - Y_hat_base; % [180 x 255]
end

%% For each K, compute fraction of residual variance explained by each source
r2_static_ext = nan(n_subj, n_tasks, n_K);
r2_delay = nan(n_subj, n_tasks, n_K);
r2_random = nan(n_subj, n_tasks, n_K);

fprintf('\n=== Residual variance explained, K=1..%d ===\n', max(K_values));
for ti = 1:n_tasks
    R = Residual_tasks{ti}; % [180 x 255]
    var_R = sum(R.^2, 1); % total residual sum-of-squares per subject

    for ki = 1:n_K
        K = K_values(ki);

        % 1. Next static modes (21:20+K)
        Psi_ext = psi_static(:, N_base+1 : N_base+K);
        R_hat = Psi_ext * (Psi_ext' * R);
        r2_static_ext(:,ti,ki) = 1 - sum((R - R_hat).^2, 1)' ./ var_R';

        % 2. Delay-informed modes (1:K)
        Psi_d = psi_delay(:, 1:K);
        R_hat = Psi_d * (Psi_d' * R);
        r2_delay(:,ti,ki) = 1 - sum((R - R_hat).^2, 1)' ./ var_R';

        % 3. Random orthonormal directions (averaged over draws)
        r2_rand_draws = nan(n_subj, n_random_draws);
        for d = 1:n_random_draws
            Psi_r = psi_random{d}(:, 1:K);
            R_hat = Psi_r * (Psi_r' * R);
            r2_rand_draws(:,d) = 1 - sum((R - R_hat).^2, 1)' ./ var_R';
        end
        r2_random(:,ti,ki) = mean(r2_rand_draws, 2);
    end
    fprintf('  task %d/%d done\n', ti, n_tasks);
end

% average across tasks
r2_static_ext_avg = squeeze(mean(r2_static_ext, 2)); % [255 x n_K]
r2_delay_avg = squeeze(mean(r2_delay, 2));
r2_random_avg = squeeze(mean(r2_random, 2));

mean_static_ext = mean(r2_static_ext_avg, 1, 'omitnan');
mean_delay = mean(r2_delay_avg, 1, 'omitnan');
mean_random = mean(r2_random_avg, 1, 'omitnan');

fprintf('\nK=5:  static_ext R2=%.4f, delay R2=%.4f, random R2=%.4f\n', mean_static_ext(5), mean_delay(5), mean_random(5));
fprintf('K=10: static_ext R2=%.4f, delay R2=%.4f, random R2=%.4f\n', mean_static_ext(10), mean_delay(10), mean_random(10));
fprintf('K=20: static_ext R2=%.4f, delay R2=%.4f, random R2=%.4f\n', mean_static_ext(20), mean_delay(20), mean_random(20));
fprintf('K=40: static_ext R2=%.4f, delay R2=%.4f, random R2=%.4f\n', mean_static_ext(40), mean_delay(40), mean_random(40));

%% Paired stats at K=5 and K=10 (delay vs static_ext, delay vs random)
fprintf('\n=== Paired t-tests: delay vs. controls, Bonferroni x4 ===\n');
test_Ks = [5, 10];
for K = test_Ks
    ki = find(K_values == K);
    x = r2_delay_avg(:,ki);
    y_ext = r2_static_ext_avg(:,ki);
    y_rand = r2_random_avg(:,ki);
    valid1 = ~isnan(x) & ~isnan(y_ext);
    valid2 = ~isnan(x) & ~isnan(y_rand);
    [~, p1, ~, st1] = ttest(x(valid1), y_ext(valid1));
    [~, p2, ~, st2] = ttest(x(valid2), y_rand(valid2));
    p1c = min(p1*4, 1); p2c = min(p2*4, 1);
    fprintf('  K=%2d: delay vs next-static-modes: diff=%+.4f, t(%d)=%.3f, p_corr=%.3e\n', K, mean(x(valid1)-y_ext(valid1)), st1.df, st1.tstat, p1c);
    fprintf('  K=%2d: delay vs random:            diff=%+.4f, t(%d)=%.3f, p_corr=%.3e\n', K, mean(x(valid2)-y_rand(valid2)), st2.df, st2.tstat, p2c);
end

%% Figure
fg = figure('Name', 'Residual variance explained', 'Visible', 'off');
plot(K_values, mean_static_ext, 'LineWidth', 1.8); hold on
plot(K_values, mean_delay, 'LineWidth', 1.8);
plot(K_values, mean_random, '--', 'LineWidth', 1.5);
grid on; xlabel('Additional dimensions (K)'); ylabel('Fraction of static-N=20 residual variance explained');
legend({'More static modes (21..20+K)', 'Delay-informed modes (1..K)', 'Random orthonormal directions'}, 'Location', 'southeast');
title('Does delay explain what static EDR+LR (N=20) leaves as residual?');
exportgraphics(fg, fullfile(fig_out_dir, 'residual_variance_explained.png'));

%% Save
save(fullfile(results_out_dir, 'residual_variance_test_results.mat'), ...
    'K_values', 'r2_static_ext_avg', 'r2_delay_avg', 'r2_random_avg', ...
    'mean_static_ext', 'mean_delay', 'mean_random', 'N_base', '-v7.3');
fprintf('\nSaved figure -> %s\n', fullfile(fig_out_dir, 'residual_variance_explained.png'));
fprintf('Saved results -> %s\n', fullfile(results_out_dir, 'residual_variance_test_results.mat'));
fprintf('=== Residual variance test complete ===\n');
