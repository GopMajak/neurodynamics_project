function r = order_parameter( ang, NN )
% Computes the Kuramoto order parameter: a measure of how in-sync a group
% of oscillators is at each point in time (0 = totally desynchronized,
% 1 = perfectly in phase).
%
% Inputs:
%   ang - phase angles over time (rows = time, columns = nodes)
%   NN  - number of nodes
%
% Output:
%   r - order parameter over time, values in [0, 1]

r = (1/NN) * sum( exp(1i*ang), 2 );
r = abs(r);
