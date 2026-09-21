function [rotp,rotm] = rotation_operator( th, ang )
% Measures how much a set of phases rotates together in the positive vs.
% negative direction, given each node's angular position in space.
%
% Inputs:
%   th  - phase values over time (rows = time, columns = nodes)
%   ang - angular position of each node
%
% Outputs:
%   rotp - amount of rotation in the positive direction
%   rotm - amount of rotation in the negative direction

assert( ismatrix(th) ); assert( isvector(ang) );

rotp = nan( 1, size(th,1) ); rotm = nan( 1, size(th,1) );
for ii = 1:size(th,1)
    rotp(ii) = (1/size(th,2)) * sum( exp(1i*th(ii,:)) .* exp(1i*ang) );
    rotm(ii) = (1/size(th,2)) * sum( exp(-1i*th(ii,:)) .* exp(1i*ang) );
end
rotp = abs(rotp); rotm = abs(rotm);
