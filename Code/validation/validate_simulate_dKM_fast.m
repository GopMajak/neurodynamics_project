%% Validate simulate_dKM_fast.m against the reference simulate_dKM.m
%
% "Original Project Code/brain_connectome_harmonics/simulation/simulate_dKM_fast.m"
% claims to be a vectorized, faster version of Budzinski et al. 2023's
% reference simulate_dKM.m that "gives identical results". This checks
% that claim directly on a small synthetic network before trusting it for
% real work, using the same small ring-graph setup already validated
% against Budzinski's own analytical prediction
% (Code/validation/validate_km_waves_ring.m).

clear; clc;
this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir));
budzinski_dir = fullfile(project_root, 'Reference Code Material', 'budzinski2023_dKM');
recovered_dir = fullfile(project_root, 'Original Project Code', 'brain_connectome_harmonics');

addpath(fullfile(budzinski_dir, 'graphs'));      % ring_graph, delay_matrix (reference)
addpath(fullfile(budzinski_dir, 'simulation'));  % simulate_dKM (reference, slow)
addpath(fullfile(recovered_dir, 'simulation'));  % simulate_dKM_fast (to validate)
addpath(fullfile(recovered_dir, 'functions'));   % delay_matrix_from_distance

rng(1);
N = 60; k = 20; % small ring, keeps the slow reference tractable
dt = 0.001; T = 1.5; t = 0:dt:T;
f_mu = 10;
omega = f_mu*2*pi*ones(N,1);
kappa = 0.5;
speed = 400;

a = ring_graph(N, k);
tau_ref = delay_matrix(a, speed, dt); % reference: ring-specific delay matrix (integer timesteps already)
theta0 = 2*pi*(rand(N,1) - 0.5);

fprintf('Simulating with reference simulate_dKM.m (slow, nested loops)...\n');
tic; theta_ref = simulate_dKM(a, tau_ref, omega, kappa, theta0, t, dt, 'euler'); t_ref = toc;
fprintf('  %.2f s\n', t_ref);

fprintf('Simulating with simulate_dKM_fast.m (vectorized)...\n');
tic; theta_fast = simulate_dKM_fast(a, tau_ref, omega, kappa, theta0, t, dt, 'euler'); t_fast = toc;
fprintf('  %.2f s\n', t_fast);

max_abs_diff = max(abs(theta_ref(:) - theta_fast(:)));
fprintf('\nMax absolute difference between reference and fast: %.3e\n', max_abs_diff);
fprintf('Speedup: %.1fx\n', t_ref / t_fast);

if max_abs_diff < 1e-9
    fprintf('PASS: simulate_dKM_fast.m is numerically identical to the reference.\n');
else
    fprintf('FAIL: outputs diverge -- do not trust simulate_dKM_fast.m without further checking.\n');
end
