%% Hybrid basis test (Option A): does delay complement geometry+long-range?
%
% Reframing from "does a delay-informed basis replace the static one" (H1,
% now well-triangulated as NOT supported) to "does delay carry
% COMPLEMENTARY information geometry+long-range alone misses" -- i.e.
% not arguing against Vohryzek/Pang, but testing whether delay should be
% incorporated ALONGSIDE their static EDR+LR modes rather than instead of
% them.
%
% At a FIXED total budget of N=20 modes (the proposal's primary
% comparison point), build a hybrid basis: N_static static EDR+LR modes
% (priority order, i.e. kept in full) + N_delay delay-informed
% eigenvector-basis modes (orthogonalized against the static ones via QR,
% so only their UNIQUE/complementary component contributes), with
% N_static + N_delay = 20 always. Sweep N_delay from 0 (pure static,
% reproduces the already-established per-subject baseline exactly -- a
% built-in correctness check) to 20 (pure delay, reproduces the earlier
% delay-only result). If delay is genuinely complementary, accuracy
% should INCREASE somewhere in between; if delay modes are just worse
% substitutes for static ones (as everything so far suggests), accuracy
% should monotonically decrease as N_delay increases.
%
% Per-subject (255 subjects, 7 tasks), matching run_per_subject_statistics.m's
% rigor, not group-averaged.

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
EDR_LRE_km  = EDR_LRE_parc; EDR_LRE_km(1:N3+1:end)  = 0;

%% Build static + delay-informed bases (EDR+LR, the "geometry + long-range" graph)
fprintf('Computing static (Laplacian) EDR+LR eigenmodes...\n');
[psi_static, ~] = normalized_laplacian_eigenmodes(EDR_LRE_km, num_parcels);

fprintf('Computing W_delay eigenvector basis (f_mu=10Hz, speed=5m/s, kappa=6 -- standard reference point)...\n');
dt = 1e-3; f_mu = 10; omega = f_mu*2*pi*ones(N3,1); kappa = 6; speed_basis = 5;
tau = delay_matrix_from_distance(rr_parc, speed_basis, dt);
[psi_delay, ~] = build_delay_eigenvector_basis(EDR_LRE_km, tau, omega, kappa, dt);

%% Load per-subject task data (7 representative tasks)
fprintf('Loading task data...\n');
data = load(fullfile(vohryzek_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
fieldNames = fieldnames(data.zstat);
representTask = [3, 11, 19, 30, 41, 44, 47];
n_tasks = numel(representTask);
N_total = 20; % the proposal's primary comparison point

Y_tasks = cell(1, n_tasks);
for ti = 1:n_tasks
    Y_tasks{ti} = parcellate_average(data.zstat.(fieldNames{representTask(ti)})(cortex_ind, :), parc_cortex); % [180 x 255]
end
n_subj = size(Y_tasks{1}, 2);

%% Sweep the static/delay split at fixed total N=20
N_delay_values = 0:N_total;
n_splits = numel(N_delay_values);
acc_avg = nan(n_subj, n_splits); % per-subject accuracy (7-task average), per split

fprintf('\n=== Hybrid basis sweep: N_static + N_delay = %d, N_delay = 0..%d ===\n', N_total, N_total);
for si = 1:n_splits
    N_delay = N_delay_values(si);
    N_static = N_total - N_delay;

    if N_static > 0 && N_delay > 0
        candidate = [psi_static(:,1:N_static), psi_delay(:,1:N_delay)];
        [Psi, ~] = qr(candidate, 0); % static modes prioritized (first columns), delay modes orthogonalized against them
    elseif N_static > 0
        Psi = psi_static(:,1:N_static);
    else
        Psi = psi_delay(:,1:N_delay);
    end

    acc_this_split = nan(n_subj, n_tasks);
    for ti = 1:n_tasks
        Y = Y_tasks{ti};
        A = Psi' * Y;
        Y_hat = Psi * A;
        acc_this_split(:,ti) = corr_columns(Y, Y_hat)';
    end
    acc_avg(:,si) = mean(acc_this_split, 2);

    fprintf('  N_static=%2d, N_delay=%2d | mean accuracy = %.4f\n', N_static, N_delay, mean(acc_avg(:,si), 'omitnan'));
end

%% Sanity checks against already-established results
fprintf('\nSanity check -- N_delay=0 (pure static) mean accuracy: %.4f (expect ~0.6325, matching run_per_subject_statistics.m)\n', ...
    mean(acc_avg(:,1), 'omitnan'));
fprintf('Sanity check -- N_delay=20 (pure delay) mean accuracy: %.4f (expect ~0.6100, matching run_per_subject_statistics.m)\n', ...
    mean(acc_avg(:,end), 'omitnan'));

%% Statistics: does ANY hybrid split beat pure static (N_delay=0), paired t-test?
baseline = acc_avg(:,1);
fprintf('\n=== Paired t-test: each hybrid split vs. pure static (N_delay=0), Bonferroni x%d ===\n', n_splits-1);
best_p = 1; best_split = 0; best_diff = 0;
for si = 2:n_splits
    x = acc_avg(:,si); y = baseline;
    valid = ~isnan(x) & ~isnan(y);
    [~, p, ~, st] = ttest(x(valid), y(valid));
    p_corr = min(p * (n_splits-1), 1);
    diff = mean(x(valid) - y(valid));
    fprintf('  N_delay=%2d: mean diff = %+.4f, t(%d)=%.3f, p_corr=%.3e %s\n', ...
        N_delay_values(si), diff, st.df, st.tstat, p_corr, ternary_str(diff > 0 && p_corr < 0.05, '<-- SIGNIFICANT IMPROVEMENT', ''));
    if diff > best_diff
        best_diff = diff; best_p = p_corr; best_split = N_delay_values(si);
    end
end
fprintf('\nBest-performing split: N_delay=%d (diff=%+.4f vs pure static, p_corr=%.3e)\n', best_split, best_diff, best_p);
if best_diff > 0 && best_p < 0.05
    fprintf('>>> A hybrid split SIGNIFICANTLY beats pure static -- delay adds complementary information.\n');
else
    fprintf('>>> No hybrid split significantly beats pure static -- no evidence delay adds complementary information at this operating point.\n');
end

%% Figure
fg = figure('Name', 'Hybrid basis sweep', 'Visible', 'off');
mean_acc = mean(acc_avg, 1, 'omitnan');
sem_acc = std(acc_avg, 0, 1, 'omitnan') / sqrt(n_subj);
errorbar(N_delay_values, mean_acc, sem_acc, 'o-', 'LineWidth', 1.8);
grid on; xlabel('Number of delay-informed modes in the 20-mode basis (N_{delay})');
ylabel('Mean per-subject accuracy (7-task average)');
title('Hybrid static-EDR+LR / delay-informed basis: does delay complement geometry+long-range?');
exportgraphics(fg, fullfile(fig_out_dir, 'hybrid_basis_sweep.png'));

%% Save
save(fullfile(results_out_dir, 'hybrid_basis_test_results.mat'), ...
    'N_delay_values', 'acc_avg', 'mean_acc', 'sem_acc', 'best_split', 'best_diff', 'best_p', 'N_total', '-v7.3');
fprintf('\nSaved figure -> %s\n', fullfile(fig_out_dir, 'hybrid_basis_sweep.png'));
fprintf('Saved results -> %s\n', fullfile(results_out_dir, 'hybrid_basis_test_results.mat'));
fprintf('=== Hybrid basis test complete ===\n');

function s = ternary_str(cond, a, b)
    if cond, s = a; else, s = b; end
end
