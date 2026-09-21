%% Sweep natural-frequency heterogeneity to escape the rank-2 shared-carrier limit
%
% Confirmed (sweep_coupling_strength.m, sliding_window_complexity.m):
% with every node sharing one natural frequency, cos(theta(t)) is
% mathematically confined to ~2 effective dimensions regardless of
% coupling strength or time window. This tests whether giving each node
% its own natural frequency (Gaussian spread around the same 10 Hz mean)
% breaks that constraint, at a fixed, moderate coupling (kappa=6, the
% same value used throughout Stage 3).

clear; clc;
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

%% Sweep frequency spread
dt = 1e-3;
T_transient = 5.0; T_steady = 5.0;
t_full = 0:dt:(T_transient + T_steady);
n_transient_steps = round(T_transient/dt);

f_mu = 10;
kappa = 6; speed_basis = 5;
tau = delay_matrix_from_distance(rr_parc, speed_basis, dt);

sigma_f_values = [0, 0.1, 0.5, 1, 2, 5]; % Hz, std of frequency spread around f_mu
n_s = numel(sigma_f_values);
R_edr = nan(1,n_s); R_edrlr = nan(1,n_s);
peff_edr = nan(1,n_s); peff_edrlr = nan(1,n_s);

fprintf('%10s | %10s %10s | %12s %12s\n', 'sigma_f', 'R (EDR)', 'R (EDR+LR)', 'P_eff (EDR)', 'P_eff (EDR+LR)');
for si = 1:n_s
    rng(2); % fixed seed for the frequency draw, so only sigma_f varies
    omega_i = 2*pi*(f_mu + sigma_f_values(si)*randn(N3,1));
    rng(1); theta0 = 2*pi*rand(N3,1); % same initial condition as before

    theta_edr = simulate_dKM_fast(EDR_conn_km, tau, omega_i, kappa, theta0, t_full, dt, 'euler');
    theta_edrlr = simulate_dKM_fast(EDR_LRE_km, tau, omega_i, kappa, theta0, t_full, dt, 'euler');

    R_e = order_parameter(theta_edr(n_transient_steps+1:end,:), N3);
    R_el = order_parameter(theta_edrlr(n_transient_steps+1:end,:), N3);
    R_edr(si) = mean(R_e); R_edrlr(si) = mean(R_el);

    X_edr = cos(theta_edr(n_transient_steps+1:end, :))';
    X_edrlr = cos(theta_edrlr(n_transient_steps+1:end, :))';
    sv_edr = svd(X_edr, 'econ'); sv_edrlr = svd(X_edrlr, 'econ');
    peff_edr(si) = (sum(sv_edr.^2))^2 / sum(sv_edr.^4);
    peff_edrlr(si) = (sum(sv_edrlr.^2))^2 / sum(sv_edrlr.^4);

    fprintf('%10.2f | %10.4f %10.4f | %12.2f %12.2f\n', sigma_f_values(si), R_edr(si), R_edrlr(si), peff_edr(si), peff_edrlr(si));
end

fg = figure('Name', 'Frequency heterogeneity sweep', 'Visible', 'off');
subplot(1,2,1)
plot(sigma_f_values, R_edr, 'o-', 'LineWidth', 1.5); hold on
plot(sigma_f_values, R_edrlr, 's-', 'LineWidth', 1.5); grid on
xlabel('Frequency spread \sigma_f (Hz)'); ylabel('Order parameter R'); ylim([0 1])
legend('EDR', 'EDR+LR', 'Location', 'best'); title('Synchronization vs. frequency spread')
subplot(1,2,2)
plot(sigma_f_values, peff_edr, 'o-', 'LineWidth', 1.5); hold on
plot(sigma_f_values, peff_edrlr, 's-', 'LineWidth', 1.5); grid on
xlabel('Frequency spread \sigma_f (Hz)'); ylabel('Effective # modes (P_{eff})')
legend('EDR', 'EDR+LR', 'Location', 'best'); title('Basis richness vs. frequency spread')
exportgraphics(fg, fullfile(fig_out_dir, 'frequency_sweep.png'));

save(fullfile(project_root, 'Results', 'novel_delay_extension', 'frequency_sweep.mat'), ...
    'sigma_f_values', 'R_edr', 'R_edrlr', 'peff_edr', 'peff_edrlr');
fprintf('\nSaved -> %s\n', fullfile(fig_out_dir, 'frequency_sweep.png'));
