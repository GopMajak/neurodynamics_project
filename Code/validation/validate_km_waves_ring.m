%% Validation: reproduce Budzinski et al. 2023's ring-graph dKM test case
%
% This is a lightly patched copy of Reference Code Material/budzinski2023_dKM/
% km_waves_ring.m (small, self-contained, N=100 nodes -- safe to run
% locally). Patches vs. the original:
%   - path setup relative to this project instead of assuming cwd is the
%     budzinski2023_dKM folder
%   - figures saved to PNG (headless-safe) instead of just displayed
%   - the interactive movie/animation section (pause()-based) removed --
%     it produces no artifact useful for batch validation
%   - added a quantitative check of the paper's core claim: that
%     conduction delay shifts the dominant eigenmode (largest real part
%     of the delayed coupling operator's eigenvalues) away from the
%     trivial synchronous mode (mode index 1, k=0)
%
% All simulation/analysis logic (ring_graph, delay_matrix,
% circulant_eigensystem, simulate_KM, simulate_dKM) is untouched,
% original Budzinski et al. 2023 code.

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir));
budzinski_dir = fullfile(project_root, 'Reference Code Material', 'budzinski2023_dKM');
addpath(fullfile(budzinski_dir, 'helper_functions'));
addpath(fullfile(budzinski_dir, 'analysis'));
addpath(fullfile(budzinski_dir, 'simulation'));
addpath(fullfile(budzinski_dir, 'graphs'));

fig_out_dir = fullfile(project_root, 'Figures', 'validation');
if ~exist(fig_out_dir, 'dir'); mkdir(fig_out_dir); end

%% Parameters (identical to original km_waves_ring.m)
dt = 0.001; T = 3.0; t = 0:dt:T; f_mu = 10;
N = 100; k = 50; method = 'euler';
omega = f_mu*2*pi*ones(N,1);
epsilon = 0.5;   % coupling strength
speed = 400;     % conduction velocity

fprintf('Building ring graph (N=%d, k=%d)...\n', N, k);
a = ring_graph(N, k);

fprintf('Computing heterogeneous delay matrix (speed=%d)...\n', speed);
tau = delay_matrix(a, speed, dt);

fprintf('Solving analytical (circulant) eigensystems...\n');
A = epsilon .* a;
[v_nd, d_nd] = circulant_eigensystem(A);           % non-delayed
W = epsilon .* (exp((-1i * omega * dt) .* tau) .* a);
[v, d] = circulant_eigensystem(W);                  % delayed

theta0 = 2*pi*(rand(N,1) - 0.5);

fprintf('Simulating KM (no delay)...\n');
tic; thetaKM = simulate_KM(a, omega, epsilon, theta0, t, dt, method, 0); fprintf('  %.1f s\n', toc);

fprintf('Simulating dKM (with delay)...\n');
tic; thetaDKM = simulate_dKM(a, tau, omega, epsilon, theta0, t, dt, method); fprintf('  %.1f s\n', toc);

%% Quantitative check: does delay shift the dominant mode away from k=0?
eig_real_nd = real(diag(d_nd));
eig_real_d = real(diag(d));
[~, dominant_mode_nd] = max(eig_real_nd);
[~, dominant_mode_d] = max(eig_real_d);

fprintf('\n=== Core claim check ===\n');
fprintf('Non-delayed: dominant mode index = %d (eigenvalue real part = %.4f)\n', ...
    dominant_mode_nd, eig_real_nd(dominant_mode_nd));
fprintf('Delayed:     dominant mode index = %d (eigenvalue real part = %.4f)\n', ...
    dominant_mode_d, eig_real_d(dominant_mode_d));
if dominant_mode_nd == 1 && dominant_mode_d ~= 1
    fprintf(['PASS: delay shifts spectral dominance away from the trivial\n' ...
             '      synchronous mode (index 1), consistent with Budzinski et al. 2023.\n']);
else
    fprintf(['CHECK NEEDED: expected non-delayed dominance at mode index 1 and\n' ...
             '      delayed dominance elsewhere -- got nd=%d, d=%d. Inspect eigenvalue\n' ...
             '      plots before trusting this configuration.\n'], dominant_mode_nd, dominant_mode_d);
end

%% Figures (saved instead of displayed -- headless run)
fg1 = figure('Visible', 'off'); hold on;
imagesc(1:N, t, thetaKM); colormap bone;
xlabel('nodes'); ylabel('time (s)'); title('no delay');
set(gca, 'fontname', 'arial', 'fontsize', 14, 'linewidth', 2);
xlim([1 N]); ylim([0 3]);
ax = gca; ax.YDir = 'reverse';
exportgraphics(fg1, fullfile(fig_out_dir, 'ring_no_delay_raster.png'));

fg2 = figure('Visible', 'off'); hold on;
imagesc(1:N, t, thetaDKM); colormap bone;
xlabel('nodes'); ylabel('time (s)'); title('with delays');
set(gca, 'fontname', 'arial', 'fontsize', 14, 'linewidth', 2);
xlim([1 N]); ylim([0 3]);
ax = gca; ax.YDir = 'reverse';
exportgraphics(fg2, fullfile(fig_out_dir, 'ring_with_delay_raster.png'));

fg3 = figure('Visible', 'off'); hold on;
plot(real(diag(d_nd)), imag(diag(d_nd)), 'o-', 'MarkerEdgeColor', '#0072BD', 'MarkerSize', 10, 'linewidth', 2);
plot(real(diag(d)), imag(diag(d)), 's-', 'MarkerEdgeColor', '#D95319', 'MarkerSize', 10, 'linewidth', 2);
legend({'no delay', 'with delay'});
xlabel('real part'); ylabel('imaginary part'); title('eigenvalues in the complex plane');
set(gca, 'fontname', 'arial', 'fontsize', 14, 'linewidth', 2);
exportgraphics(fg3, fullfile(fig_out_dir, 'eigenvalues_complex_plane.png'));

fg4 = figure('Visible', 'off'); hold on;
plot(real(diag(d_nd)), 'o-', 'MarkerEdgeColor', '#0072BD', 'MarkerSize', 8, 'linewidth', 2);
plot(real(diag(d)), 's-', 'MarkerEdgeColor', '#D95319', 'MarkerSize', 8, 'linewidth', 2);
legend({'no delay', 'with delay'});
ylabel('real part'); xlabel('modes'); title('eigenvalue real part by mode');
set(gca, 'fontname', 'arial', 'fontsize', 14, 'linewidth', 2);
exportgraphics(fg4, fullfile(fig_out_dir, 'eigenvalues_real_part_by_mode.png'));

fprintf('\nFigures saved -> %s\n', fig_out_dir);
fprintf('=== Validation run complete ===\n');
