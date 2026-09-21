%% Sliding-window participation ratio: where does spatial complexity live?
%
% The coupling sweep showed P_eff stuck near 2 in the (discarded)
% steady-state tail regardless of kappa. This instead tracks P_eff across
% the WHOLE simulation in a sliding window, to find where (if anywhere)
% the dynamics are genuinely richer than the ~rank-2 free-rotation limit
% -- expected in the transient, based on the wave-like patterns visible
% early in the ring-graph validation raster plot.

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

%% Simulate once, long enough to see the full transient-to-steady-state arc
dt = 1e-3; T = 10.0; t_full = 0:dt:T;
f_mu = 10; omega = f_mu*2*pi*ones(N3,1);
kappa = 6; speed_basis = 5;
rng(1); theta0 = 2*pi*rand(N3,1);
tau = delay_matrix_from_distance(rr_parc, speed_basis, dt);

fprintf('Simulating EDR-only and EDR+LR (kappa=%g, speed=%g m/s, T=%g s)...\n', kappa, speed_basis, T);
theta_edr = simulate_dKM_fast(EDR_conn_km, tau, omega, kappa, theta0, t_full, dt, 'euler');
theta_edrlr = simulate_dKM_fast(EDR_LRE_km, tau, omega, kappa, theta0, t_full, dt, 'euler');

%% Sliding window participation ratio + order parameter
win_size = round(0.5/dt);  % 0.5 s windows
win_step = round(0.1/dt);  % slide by 0.1 s
starts = 1:win_step:(length(t_full) - win_size);
win_times = t_full(starts) + 0.25; % window center

peff_edr_t = nan(1, numel(starts));
peff_edrlr_t = nan(1, numel(starts));
R_edr_t = nan(1, numel(starts));
R_edrlr_t = nan(1, numel(starts));

for wi = 1:numel(starts)
    idx = starts(wi):(starts(wi)+win_size-1);
    Xe = cos(theta_edr(idx,:))';
    Xel = cos(theta_edrlr(idx,:))';
    sv_e = svd(Xe, 'econ'); sv_el = svd(Xel, 'econ');
    peff_edr_t(wi) = (sum(sv_e.^2))^2 / sum(sv_e.^4);
    peff_edrlr_t(wi) = (sum(sv_el.^2))^2 / sum(sv_el.^4);
    R_edr_t(wi) = mean(order_parameter(theta_edr(idx,:), N3));
    R_edrlr_t(wi) = mean(order_parameter(theta_edrlr(idx,:), N3));
end

fprintf('\nPeak P_eff (EDR):   %.2f at t=%.2fs\n', max(peff_edr_t), win_times(find(peff_edr_t==max(peff_edr_t),1)));
fprintf('Peak P_eff (EDR+LR): %.2f at t=%.2fs\n', max(peff_edrlr_t), win_times(find(peff_edrlr_t==max(peff_edrlr_t),1)));

fg = figure('Name', 'Sliding-window complexity', 'Visible', 'off');
subplot(2,1,1)
plot(win_times, R_edr_t, 'LineWidth', 1.5); hold on
plot(win_times, R_edrlr_t, 'LineWidth', 1.5); grid on
ylabel('Order parameter R'); ylim([0 1])
legend('EDR', 'EDR+LR', 'Location', 'best')
title(sprintf('kappa=%g, speed=%g m/s', kappa, speed_basis))
subplot(2,1,2)
plot(win_times, peff_edr_t, 'LineWidth', 1.5); hold on
plot(win_times, peff_edrlr_t, 'LineWidth', 1.5); grid on
xlabel('Time (s)'); ylabel('Effective # modes (P_{eff})')
legend('EDR', 'EDR+LR', 'Location', 'best')
exportgraphics(fg, fullfile(fig_out_dir, 'sliding_window_complexity.png'));
fprintf('Saved -> %s\n', fullfile(fig_out_dir, 'sliding_window_complexity.png'));
