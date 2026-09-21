function theta_final = simulate_dKM_batch( w, tau, omega, kappa, theta0_batch, n_steps, dt )
%
% Same delayed Kuramoto model as simulate_dKM.m, but runs R separate
% simulations at once (each with its own initial conditions, sharing the
% same network/delays/frequencies) and only returns the final phase
% state -- useful when you don't need the full time series, just the
% end result for many random starting conditions.
%
% Keeps only a short rolling window of phase history (a "circular
% buffer") instead of the full time series, since only the most recent
% history is ever needed. Requires tau >= 1 everywhere w > 0 (a
% connection can't have zero delay).
%
% INPUT
% w - weight matrix (NxN)
% tau - delay matrix (timesteps, integer, NxN)
% omega - frequencies (Nx1) (rad/s)
% kappa - coupling strength
% theta0_batch - initial conditions (N x R)
% n_steps - total number of Euler steps to integrate (global step count)
% dt - timestep (s)
%
% OUTPUT
% theta_final - wrapped phase state at step n_steps (N x R)
%

N = size(w,1); R = size(theta0_batch,2);
max_tau = max(tau(:));
start_timestep = max_tau + 1;
buflen = max_tau + 2; % buffer only needs to hold the longest delay, plus one

B = zeros(buflen, N, R);
tvec = (0:start_timestep-1)' * dt;
for jj = 1:N
    B(1:start_timestep, jj, :) = theta0_batch(jj,:) + tvec .* omega(jj);
end

rowidx = repmat( (1:N)', 1, N );
omega_col = omega(:);

pos = @(g) mod(g-1, buflen) + 1; % maps a timestep to its slot in the circular buffer

for ii = start_timestep:n_steps

    idx2d = pos( ii - tau );                 % N x N, buffer slot for each delayed timestep
    lin2d = idx2d + (rowidx - 1) * buflen;    % N x N, linear index into the buffer
    prev_pos = pos(ii-1);

    theta_prev = reshape( B(prev_pos, :, :), [N, R] );   % phases at the previous step, N x R

    planeR = reshape( 0:(R-1), [1,1,R] ) * (buflen*N);
    lin3d = lin2d + planeR;                                % extend the index across all R realizations
    M = B(lin3d);                                           % delayed phases, N x N x R (kk,jj,r)

    theta_prev_b = reshape( theta_prev, [1, N, R] );        % broadcast so each kk sees node jj's previous phase
    S = w .* sin( M - theta_prev_b );                        % coupling term for every pair, N x N x R
    dth = omega_col + kappa .* reshape( sum(S,1), [N,R] );   % total phase update, N x R

    new_theta = theta_prev + dth * dt;
    B( pos(ii), :, : ) = reshape( new_theta, [1,N,R] );

end

theta_final = angle( exp( 1i * new_theta ) ); % wrap phases into [-pi, pi]
