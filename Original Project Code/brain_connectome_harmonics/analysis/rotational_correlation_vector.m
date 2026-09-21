function [rho,ang] = rotational_correlation_vector( x, p )
% Checks how much a value at each node correlates with that node's
% angular position around the sagittal plane (i.e. does x change
% systematically as you rotate around the brain).
%
% Inputs:
%   x - column vector of values, one per node
%   p - node positions in 3D space (nodes x 3)
%
% Outputs:
%   rho - correlation between x and angular position
%   ang - angular position computed for each node

assert( iscolumn(x) == 1, 'row vector input required, x' )

N = length(x); ang = nan( N, 1 ); cp = mean( p );
for ii = 1:N, ang(ii) = cart2pol( p(ii,1) - cp(1), -( p(ii,3) - cp(3) ) ); end

rho = circ_corrcl( ang, x );
