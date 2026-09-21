%% Does the W_delay eigenvector basis ever exceed static? Sweep f_mu (and verify kappa is irrelevant).
%
% kappa cannot matter for this basis: W_delay = kappa * (phase factor) .* W
% scales the whole matrix by a positive constant, which scales eigenVALUES
% by kappa but leaves eigenVECTORS unchanged (Av=lambda*v => (kappa*A)v =
% (kappa*lambda)v, same v). Verified empirically below, then not swept
% further. f_mu DOES change the eigenvectors nontrivially, since it
% multiplies the edge-varying delay tau(i,j) inside the phase term
% exp(-i*omega*dt*tau(i,j)) -- swept here across several speeds to look
% for any (f_mu, speed) combination where the eigenbasis exceeds static
% at any N, not just N=20.

clear; clc;
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
EDR_conn_km = EDR_conn;     EDR_conn_km(1:N3+1:end) = 0;
EDR_LRE_km  = EDR_LRE_parc; EDR_LRE_km(1:N3+1:end)  = 0;

%% Static baseline
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

%% Verify kappa is irrelevant (cheap sanity check)
dt = 1e-3; f_mu_ref = 10; omega_ref = f_mu_ref*2*pi*ones(N3,1);
tau_ref = delay_matrix_from_distance(rr_parc, 5, dt);
[psi_k1, ~] = build_delay_eigenvector_basis(EDR_conn_km, tau_ref, omega_ref, 1, dt);
[psi_k20, ~] = build_delay_eigenvector_basis(EDR_conn_km, tau_ref, omega_ref, 20, dt);
% compare subspaces (sign/order-robust): principal angles via SVD of psi_k1'*psi_k20
subspace_diff = norm(abs(psi_k1(:,1:20)'*psi_k20(:,1:20)) - eye(20), 'fro');
fprintf('Kappa invariance check: ||psi(kappa=1) vs psi(kappa=20), first 20 modes|| deviation = %.2e (should be ~0)\n', subspace_diff);

%% Sweep f_mu across several speeds -- look for ANY (f_mu, speed, N) exceeding static
f_mu_values = [1, 2, 5, 10, 20, 40, 80, 160];
speeds_test = [1, 4.05, 8.63, 15, 30];
kappa = 6;

best_margin_edr = -Inf; best_edr = struct();
best_margin_edrlr = -Inf; best_edrlr = struct();

fprintf('\n=== f_mu x speed sweep: max(eigenbasis_corr - static_corr) over all N ===\n');
results = [];
for fi = 1:numel(f_mu_values)
    f_mu = f_mu_values(fi);
    omega = f_mu*2*pi*ones(N3,1);
    for si = 1:numel(speeds_test)
        speed = speeds_test(si);
        tau = delay_matrix_from_distance(rr_parc, speed, dt);

        [psi_edr, ~] = build_delay_eigenvector_basis(EDR_conn_km, tau, omega, kappa, dt);
        [psi_edrlr, ~] = build_delay_eigenvector_basis(EDR_LRE_km, tau, omega, kappa, dt);

        [~, corr_edr] = reconstruct_all(psi_edr, y_tasks, max_N);
        [~, corr_edrlr] = reconstruct_all(psi_edrlr, y_tasks, max_N);
        mc_edr = mean(corr_edr, 1);
        mc_edrlr = mean(corr_edrlr, 1);

        margin_edr = max(mc_edr - mean_corr_static_edr);
        margin_edrlr = max(mc_edrlr - mean_corr_static_edrlr);
        [~, n_at_margin_edr] = max(mc_edr - mean_corr_static_edr);
        [~, n_at_margin_edrlr] = max(mc_edrlr - mean_corr_static_edrlr);
        margin20_edr = mc_edr(20) - mean_corr_static_edr(20);
        margin20_edrlr = mc_edrlr(20) - mean_corr_static_edrlr(20);

        results(end+1,:) = [f_mu, speed, margin_edr, margin_edrlr, n_at_margin_edr, n_at_margin_edrlr, margin20_edr, margin20_edrlr]; %#ok<AGROW>
        fprintf('%8g %8.2f | best: %7.4f (N=%3d) %7.4f (N=%3d) | @N=20: %7.4f %7.4f\n', ...
            f_mu, speed, margin_edr, n_at_margin_edr, margin_edrlr, n_at_margin_edrlr, margin20_edr, margin20_edrlr);

        if margin_edr > best_margin_edr
            best_margin_edr = margin_edr;
            best_edr = struct('f_mu', f_mu, 'speed', speed, 'N', n_at_margin_edr, 'margin', margin_edr);
        end
        if margin_edrlr > best_margin_edrlr
            best_margin_edrlr = margin_edrlr;
            best_edrlr = struct('f_mu', f_mu, 'speed', speed, 'N', n_at_margin_edrlr, 'margin', margin_edrlr);
        end
    end
end

fprintf('\nBest EDR margin over static:   %.4f at f_mu=%g Hz, speed=%.2f m/s, N=%d\n', ...
    best_edr.margin, best_edr.f_mu, best_edr.speed, best_edr.N);
fprintf('Best EDR+LR margin over static: %.4f at f_mu=%g Hz, speed=%.2f m/s, N=%d\n', ...
    best_edrlr.margin, best_edrlr.f_mu, best_edrlr.speed, best_edrlr.N);
if best_margin_edr > 0 || best_margin_edrlr > 0
    fprintf('>>> Eigenbasis DOES exceed static at at least one (f_mu, speed, N) combination.\n');
else
    fprintf('>>> Eigenbasis NEVER exceeds static at any tested (f_mu, speed, N) combination.\n');
end

max_margin20_edr = max(results(:,7));
max_margin20_edrlr = max(results(:,8));
fprintf('\nBest margin AT N=20 specifically -- EDR: %.4f, EDR+LR: %.4f\n', max_margin20_edr, max_margin20_edrlr);
if max_margin20_edr > 0 || max_margin20_edrlr > 0
    [~, best_row_edr] = max(results(:,7));
    [~, best_row_edrlr] = max(results(:,8));
    fprintf('  EDR best @N=20 at f_mu=%g, speed=%.2f\n', results(best_row_edr,1), results(best_row_edr,2));
    fprintf('  EDR+LR best @N=20 at f_mu=%g, speed=%.2f\n', results(best_row_edrlr,1), results(best_row_edrlr,2));
else
    fprintf('  At the proposal''s primary comparison point (N=20), eigenbasis never exceeds static anywhere in this grid.\n');
end

save(fullfile(results_out_dir, 'eigbasis_fmu_sweep.mat'), 'results', 'best_edr', 'best_edrlr', 'subspace_diff');
fprintf('Saved -> %s\n', fullfile(results_out_dir, 'eigbasis_fmu_sweep.mat'));

%% Local functions
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
