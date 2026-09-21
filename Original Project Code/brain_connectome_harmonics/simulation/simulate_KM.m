function theta = simulate_KM( w, omega, epsilon, theta0, time, dt, method, phi )
%
% Simulates the standard (non-delayed) Kuramoto model: a network of
% coupled oscillators, each pulling its neighbors' phases toward its own.
% Supports three integration methods: a simple hand-written Euler step,
% or MATLAB's ode45/ode113 solvers (both call KM.m for the ODE itself).
%
% INPUT
% w - adjacency matrix (NxN)
% omega - frequencies (Nx1) (rad/s)
% epsilon - coupling strength
% theta0 - initial condition (rad)
% time - time axis (s)
% dt - timestep (s)
% method - integration method: 'euler', 'ode45', or 'ode113'
% phi - phase-lag
%

N = size(w,1);
theta = zeros( length(time), N ); theta(1,:) = theta0;

if strcmp( method, 'euler' )
    
    for ii = 2:length(time)
        
        previous_state = theta(ii-1,:);
        
        for jj = 1:N
            dth = omega(jj) + epsilon * nansum( w(jj,:) .* sin( previous_state - previous_state(jj) - phi ), 2 );
            theta(ii,jj) = previous_state(jj) + (dth * dt);
        end
        
    end
    
elseif strcmp( method, 'ode45' )
    
    opts = odeset( 'reltol', 1e-10, 'abstol', 1e-10 );
    [ ~ , theta ] = ode45( @(t,y) KM( time, y, N, omega, epsilon, w, phi ), time, theta0, opts );
    
elseif strcmp( method, 'ode113' )
    
    opts = odeset( 'reltol', 1e-10, 'abstol', 1e-10 );    
    [ ~ , theta ] = ode113( @(t,y) KM( time, y, N, omega, epsilon, w, phi ), time, theta0, opts );
    
end

theta = angle( exp( 1i*theta ) ); % wrap phases into [-pi, pi]
