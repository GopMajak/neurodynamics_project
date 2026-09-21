function [o] = KM( t, y, N, omega, k, w, alpha )
%
% Right-hand side of the standard Kuramoto model's ODE, for use with
% MATLAB's ode45/ode113 solvers (called from simulate_KM.m). Computes the
% rate of change of each oscillator's phase given the current phases y.
%
% INPUT
% N         number of oscillators
% omega     vector containing instantaneous angular frequencies
% k         coupling strength
% w         N x N weight matrix (double)
%
% OUTPUT
% o         rate of change of each oscillator's phase (dtheta/dt)
%

o = zeros( N, 1 );
theta = y; % current phases

for ii = 1:N
	idx = ~isnan(w(ii,:)); % skip any missing (NaN) connections
	o(ii) = omega(ii) + k * sum( w(ii,idx) * ( sin( theta(idx) - theta(ii) - alpha ) ) );
end
