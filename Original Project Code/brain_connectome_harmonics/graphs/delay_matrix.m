function tau = delay_matrix( a, v, dt )
% Builds a delay matrix for a 1D ring network: how many timesteps it
% takes a signal to travel between two connected nodes, based on their
% distance around the ring and the conduction speed.
%
% Inputs:
%   a  - adjacency matrix
%   v  - conduction speed
%   dt - simulation timestep
%
% Output:
%   tau - delay matrix, in number of timesteps

N = size( a, 1 ); tau = zeros( N );
for ii = 1:N
    for jj = 1:N
        if ( a(ii,jj) > 0)
            d = abs( ii - jj ); d = min( N - d, d ); tau(ii,jj) = floor( (d/v) ./ dt );
        end
    end
end
