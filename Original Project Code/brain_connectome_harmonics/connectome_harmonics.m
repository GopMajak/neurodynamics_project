% connectome_harmonics.m
%
% Builds the "EDR+LR" connectome from Vohryzek et al. 2024 (PNAS 2025):
% start from a smooth exponential-distance-rule (EDR) baseline fit to the
% real structural connectome, then splice back in the rare long-range (LR)
% connections that are much stronger than that baseline predicts. Compares
% the eigenmodes (eigenmode = a resting/vibration pattern of the network,
% like the ringing modes of a spread-out object) of this connectome to the
% pure-geometry eigenmodes from Pang et al. 2023.
%
% The original derivation script calls a function (calc_parcellate_matrix)
% that isn't included in either paper's code release, so it's reimplemented
% here in functions/calc_parcellate_matrix.m.
%
% Stage 1: quick test at low resolution (360 parcels) to check the pipeline
%          logic is correct. Not expected to match the paper's exact numbers.
% Stage 2: repeat Stage 1 at full resolution (29,696 vertices) - this is
%          the expensive part and is what's actually used in the paper.
% Stage 3: run a Kuramoto simulation to see whether the long-range
%          connections actually change simulated brain dynamics, not just
%          the connectome's structure.

%% Setup

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Code from Vohryzek et al 2024', 'vohryzek2024_EDRLR');

addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

%% Stage 1: load data (Glasser360 parcellation, left hemisphere)

hemisphere = 'lh';
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';
parc_name = 'Glasser360';

% surface geometry, used to get parcel centroid distances
[vertices, faces] = read_vtk(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% cortex mask (which vertices are actually cortex vs. medial wall)
cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
num_vertices = length(cortex);

% parcellation labels (given at full-vertex resolution; keep only cortex vertices)
parc = dlmread(fullfile(edrlr_data_dir, 'Data', 'parcellations', ...
    sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
parcels = unique(parc_cortex(parc_cortex>0));
num_parcels = length(parcels);
fprintf('Loaded %s parcellation: %d parcels\n', parc_name, num_parcels);

% empirical group-average structural connectome (29,696 x 29,696)
% one copy of this file is corrupted/truncated, so load the other copy instead
load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', ...
    'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');

%% Stage 1: parcel centroid distances and parcellated connectome

vertices_cortex = surface_midthickness.vertices(cortex_ind, :);
centroids = zeros(num_parcels, 3);
for p = 1:num_parcels
    centroids(p,:) = mean(vertices_cortex(parc_cortex==parcels(p), :), 1);
end
rr_parc = squareform(pdist(centroids)); % [num_parcels x num_parcels] Euclidean distance

connectome_parc = calc_parcellate_matrix(parc_cortex, avgSC_L);
C_parc = connectome_parc / max(connectome_parc(:));
clear avgSC_L

%% Stage 1: bin connectivity by distance and fit the exponential distance rule

% fewer bins than the paper's 400: with only ~360 parcels there are far
% fewer distance pairs than at vertex resolution, so more bins would leave
% most of them empty
NR = 60;
NSTD = 3;      % how many SDs above the local mean counts as "long-range"
DistRange = 40; % mm; minimum distance to call a strong connection "long-range"

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

% skip near-zero-distance bins (self/adjacent-parcel geometry) and empty bins
fit_start = find(xcoor >= 10, 1);
fit_ind = fit_start:NR;
fit_ind = fit_ind(~isnan(ycoor2(fit_ind)));

% the original script uses lsqcurvefit (needs the Optimization Toolbox);
% fminsearch (built into base MATLAB) minimizes the same squared error instead
expfunc = @(A, x) (A(1)*exp(-A(2)*x));
sse = @(A) sum((expfunc(A, xcoor(fit_ind)) - ycoor2(fit_ind)).^2);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, 'Display', 'off');
A0 = [0.15, 0.18];
Afit = fminsearch(sse, A0, options);
lambda = Afit(2);
yl = Afit(1)*exp(-Afit(2)*xcoor);

fprintf('EDR fit (Glasser360, validation only): A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

%% Stage 1: detect long-range exceptions and build the EDR+LR connectome

Clong = zeros(num_parcels, num_parcels);   % strong AND far - the actual LR exceptions
Clong_all = zeros(num_parcels, num_parcels); % all connections in the exception distance range, for plotting

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

%% Stage 1: plot the validation results

figure('Name', 'Stage 1 validation - EDR fit');
errorbar(xcoor, cellfun(@(x) mean(x,'omitnan'), sc_density), cellfun(@(x) std(x,'omitnan'), sc_density), 'o');
hold on
plot(xcoor, yl, 'r-', 'linewidth', 2)
xlabel('Distance (mm)'); ylabel('Connection strength (normalized)')
legend('binned SC (mean +/- SD)', 'EDR fit')
title(sprintf('Glasser360 validation: lambda = %.4f /mm', lambda))
grid on

figure('Name', 'Stage 1 validation - Connectome comparison');
subplot(2,2,1); imagesc(C_parc); axis square; colorbar; title('Connectome (parcellated)')
subplot(2,2,2); imagesc(EDR_conn); axis square; colorbar; title('EDR fit')
subplot(2,2,3); imagesc(EDR_LRE_parc); axis square; colorbar; title('EDR+LR')
subplot(2,2,4); imagesc(Clong>0); axis square; colorbar; title(sprintf('Long-range exceptions (n=%d)', nnz(triu(Clong,1))))
colormap(flipud(bone))

% keep these parcel-resolution variables around for Stage 3's simulation
% below (named differently from Stage 2's vertex-resolution versions so
% Stage 2 doesn't silently overwrite them)
clear Clong Clong_all index_parc sc_density sc_density_i sc_density_j

%% Stage 2: full vertex-resolution EDR+LR connectome derivation
%
% Same idea as Stage 1, but now at the real 29,696-vertex resolution used
% in the paper. This is the expensive part: several multi-GB matrices and
% a dense eigendecomposition. Variables are cleared as we go to keep
% memory usage manageable, and results are checkpointed to results/ so
% this slow step doesn't need to be redone if something later fails.

if ~exist('results', 'dir'); mkdir('results'); end

num_modes = 200;
NR = 400; NRini = 20; NRfin = 380; NSTD = 3; DistRange = 40; % as in the original paper

fprintf('\n=== STAGE 2: full vertex-resolution EDR+LR derivation ===\n');

% empirical connectome (Pang's intact copy, single precision, 29696x29696)
load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', ...
    'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');
C = avgSC_L / max(avgSC_L(:));
clear avgSC_L

% vertex-level Euclidean distances (single precision to halve memory vs. double)
tic
rr = squareform(pdist(single(vertices_cortex)));
fprintf('Computed %dx%d vertex distance matrix in %.1f s\n', size(rr,1), size(rr,2), toc);

% bin distances and fit the EDR curve, using accumarray instead of a
% per-bin loop - much faster and avoids storing every matrix entry in a cell array
range_dist = max(rr(:));
delta = range_dist / NR;
xcoor = delta/2 + delta*(0:NR-1);

index = uint16(min(floor(double(rr)/delta) + 1, NR));

tic
bin_mean = accumarray(index(:), C(:), [NR,1], @mean, single(NaN));
bin_std  = accumarray(index(:), C(:), [NR,1], @std, single(NaN));
fprintf('Binned connectivity by distance (%d bins) in %.1f s\n', NR, toc);

fit_ind = 25:NR;
fit_ind = fit_ind(~isnan(bin_mean(fit_ind)));
expfunc = @(A, x) (A(1)*exp(-A(2)*x));
sse = @(A) sum((expfunc(A, xcoor(fit_ind)) - double(bin_mean(fit_ind))').^2);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, 'Display', 'off');
Afit = fminsearch(sse, [0.15, 0.18], options);
lambda = Afit(2);
fprintf('EDR fit (full resolution): A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

% long-range exceptions: >3 SD above the local (distance-binned) mean, and >40mm away
tic
local_mean = bin_mean(index);
local_std = bin_std(index);
is_exception = (index >= NRini) & (index <= NRfin) & (rr > DistRange) & (C > local_mean + NSTD*local_std);
clear local_mean local_std
fprintf('Identified %d long-range exceptions in %.1f s\n', nnz(triu(is_exception,1)), toc);

EDR_LRE = single(Afit(1) * exp(-Afit(2) * double(rr)));
EDR_LRE(is_exception) = C(is_exception);
clear is_exception index C rr

save(fullfile('results', 'EDR_LR_connectome_full_resolution.mat'), 'EDR_LRE', 'Afit', 'lambda', '-v7.3');
fprintf('Saved full-resolution EDR+LR connectome to results/EDR_LR_connectome_full_resolution.mat\n');

%% Stage 2: EDR+LR eigenmodes (the expensive eigendecomposition step)

fprintf('\nComputing EDR+LR eigenmodes (dense eig on %dx%d)...\n', size(EDR_LRE,1), size(EDR_LRE,2));
tic
[eig_vec_temp, eig_val] = calc_network_eigenmode_lowmem(EDR_LRE, num_modes);
fprintf('Eigendecomposition done in %.1f s\n', toc);
clear EDR_LRE

% add back the medial wall vertices as zeros, matching Pang/Vohryzek's convention
eig_vec_EDRLR = zeros(num_vertices, num_modes, 'single');
eig_vec_EDRLR(cortex_ind, :) = eig_vec_temp;
clear eig_vec_temp

save(fullfile('results', 'synthetic_EDRLR_eigenmodes_fsLR_32k-lh_200.mat'), ...
    'eig_vec_EDRLR', 'eig_val', '-v7.3');
fprintf('Saved EDR+LR eigenmodes to results/synthetic_EDRLR_eigenmodes_fsLR_32k-lh_200.mat\n');

%% Stage 2: plot EDR+LR eigenmodes on the cortical surface

mode_interest = [2, 3, 4, 5, 7, 16]; % the modes highlighted in PNAS Figure 3
surface_to_plot = surface_midthickness;
data_to_plot = eig_vec_EDRLR(:, mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Stage 2 - EDR+LR eigenmodes';

%% Stage 3: time-delayed Kuramoto simulation, EDR vs EDR+LR connectomes
%
% Checks whether the long-range connections actually change simulated
% brain dynamics, not just the connectome's structure. Uses the
% time-delayed Kuramoto model (dKM) from Budzinski et al. 2023: identical
% oscillators coupled through the weighted network, where the signal from
% node k to node j arrives delayed by tau_jk = distance_jk / conduction_speed.
% Compares the EDR-only connectome (no long-range splice) against EDR+LR
% to isolate what the long-range connections actually do.
%
% Runs at Stage 1's Glasser360 (parcel) resolution, reusing the variables
% Stage 1 deliberately kept around above. Doesn't use Stage 2's full
% vertex-resolution connectome - this simulation is much more expensive
% per timestep, and Stage 2's eigendecomposition on that same matrix
% already came close to running out of memory (see stage2_log.txt).

addpath('simulation');
addpath('analysis');

if ~exist('EDR_conn', 'var') || ~exist('EDR_LRE_parc', 'var') || ~exist('rr_parc', 'var')
    error(['STAGE 3 requires STAGE 1''s parcel-resolution EDR_conn, ' ...
        'EDR_LRE_parc, and rr_parc. Run STAGE 1 first (STAGE 2 need not be run).']);
end

fprintf('\n=== STAGE 3: delayed Kuramoto on EDR vs EDR+LR (Glasser360 lh) ===\n');

N3 = size(EDR_conn, 1);

% zero out the diagonal: it holds Afit(1) (since exp(0)=1 at distance 0),
% which is just a modeling artifact, not a real self-connection
EDR_conn_km = EDR_conn;     EDR_conn_km(1:N3+1:end) = 0;
EDR_LRE_km  = EDR_LRE_parc; EDR_LRE_km(1:N3+1:end)  = 0;

%% Stage 3: single-run comparison (delayed vs non-delayed, EDR vs EDR+LR)

dt3 = 1e-3; T3 = 5.0; t3 = 0:dt3:T3;
f_mu3 = 10;                        % natural frequency (Hz), as in Budzinski's HCP run
omega3 = f_mu3*2*pi*ones(N3,1);
kappa3 = 6;                        % coupling strength; hand-tuned for this connectome, not fit to data
speed_ref3 = 5;                    % m/s, matches Budzinski's HCP reference speed

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

n_tail3 = round(2/dt3); % last 2 s
fprintf('Mean R, last 2s (delayed):    EDR = %.4f   EDR+LR = %.4f\n', ...
    mean(R_edr_delay(end-n_tail3+1:end)), mean(R_edrlr_delay(end-n_tail3+1:end)));
fprintf('Mean R, last 2s (non-delayed): EDR = %.4f   EDR+LR = %.4f\n', ...
    mean(R_edr_nodelay(end-n_tail3+1:end)), mean(R_edrlr_nodelay(end-n_tail3+1:end)));

figure('Name', 'Stage 3 - order parameter time series');
plot(t3, R_edr_delay, 'LineWidth', 1.5); hold on
plot(t3, R_edrlr_delay, 'LineWidth', 1.5);
plot(t3, R_edr_nodelay, '--', 'LineWidth', 1);
plot(t3, R_edrlr_nodelay, '--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Order parameter R(t)'); ylim([0 1]); grid on
legend('EDR (delayed)', 'EDR+LR (delayed)', 'EDR (no delay)', 'EDR+LR (no delay)', 'Location', 'best')
title(sprintf('Glasser360 lh: speed = %g m/s, \\kappa = %g', speed_ref3, kappa3))

save(fullfile('results', 'stage3_dKM_single_run.mat'), ...
    't3', 'R_edr_delay', 'R_edrlr_delay', 'R_edr_nodelay', 'R_edrlr_nodelay', ...
    'speed_ref3', 'kappa3', 'f_mu3', '-v7.3');
fprintf('Saved single-run dKM results to results/stage3_dKM_single_run.mat\n');

%% Stage 3: conduction-speed scan, EDR vs EDR+LR

n_speeds3 = 20;
speeds3 = linspace(1, 30, n_speeds3); % m/s, matches Budzinski's HCP speed-scan range
T_scan3 = 3.0; t_scan3 = 0:dt3:T_scan3;
n_half3 = floor(length(t_scan3)/2); % steady-state window: second half

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

figure('Name', 'Stage 3 - conduction speed scan');
plot(speeds3, R_edr_scan, 'o-', 'LineWidth', 1.5); hold on
plot(speeds3, R_edrlr_scan, 's-', 'LineWidth', 1.5);
xlabel('Conduction speed (m/s)'); ylabel('Steady-state order parameter R'); ylim([0 1]); grid on
legend('EDR', 'EDR+LR', 'Location', 'best')
title('Glasser360 lh: synchrony vs. conduction speed, EDR vs EDR+LR')

save(fullfile('results', 'stage3_dKM_speed_scan.mat'), ...
    'speeds3', 'R_edr_scan', 'R_edrlr_scan', 'kappa3', 'f_mu3', 'T_scan3', '-v7.3');
fprintf('Saved speed-scan dKM results to results/stage3_dKM_speed_scan.mat\n');
