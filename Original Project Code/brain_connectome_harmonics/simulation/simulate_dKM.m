function theta = simulate_dKM( w, tau, omega, kappa, theta0, time, dt, method )
%
% Simulates the time-delayed Kuramoto model (dKM): a network of coupled
% oscillators where the signal from node k to node j arrives delayed by
% tau(k,j) timesteps, instead of instantly. This is the straightforward
% (unoptimized, double for-loop) version -- see simulate_dKM_fast.m for
% a faster vectorized version with identical results.
%
% INPUT
% w - weight matrix (NxN)
% tau - delay matrix (empty if no delays) (s)
% omega - frequencies (Nx1) (rad/s)
% theta0 - initial condition (rad)
% time - time axis (s)
% dt - timestep (s)
% method - integration method
%

N = size(w,1);
theta = zeros( length(time), N ); theta(1,:) = theta0;

if strcmp( method, 'euler' )

    % before the largest delay has elapsed, there's no history to look
    % back on yet, so just run the oscillators forward at their own
    % natural frequency
    start_timestep = max( tau(:) ) + 1;
    for jj = 1:N
        theta( 1:start_timestep, jj ) = theta0(jj) + time(1:start_timestep).*omega(jj);
    end

    for ii = start_timestep:length(time)

        for jj = 1:N

            % each neighbor kk's influence on node jj uses kk's phase
            % from tau(kk,jj) steps ago, not its current phase
            dth = omega(jj);
            for kk = 1:N
                dth = dth + kappa .* w(kk,jj) .* sin( theta(ii-tau(kk,jj),kk) - theta(ii-1,jj) );
            end
            theta(ii,jj) = theta(ii-1,jj) + (dth * dt);

        end

    end

end

theta = angle( exp( 1i*theta ) ); % wrap phases into [-pi, pi]
