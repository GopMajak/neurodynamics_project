function theta = simulate_dKM_fast( w, tau, omega, kappa, theta0, time, dt, method )
%
% Same delayed Kuramoto model as simulate_dKM.m, but with the inner
% double for-loop over node pairs replaced by matrix operations. Gives
% identical results, just much faster. Requires tau >= 1 everywhere
% w > 0 (a connection can't have zero delay).
%
% INPUT
% w - weight matrix (NxN)
% tau - delay matrix (timesteps, integer)
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

    rowidx = repmat( (1:N)', 1, N );
    omega_row = omega(:)';

    for ii = start_timestep:length(time)

        lin = sub2ind( [length(time), N], ii - tau, rowidx ); % index into theta's history, per delay
        M = theta( lin );                                      % delayed phases seen by each node
        S = w .* sin( M - theta(ii-1,:) );                     % coupling term for every connection
        dth = omega_row + kappa .* sum( S, 1 );                % total phase update for each node
        theta(ii,:) = theta(ii-1,:) + (dth * dt);

    end

end

theta = angle( exp( 1i*theta ) ); % wrap phases into [-pi, pi]
