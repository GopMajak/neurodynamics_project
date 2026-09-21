function eigenvector_circle( x, p, cax, cm, num_ev )
% Plots a value at each node's 3D position, colored by that value, with
% a title showing which simulation/eigenvector number this is.
% Same idea as km_plot.m, plus a title showing the mode/simulation number.
%
% Inputs:
%   x      - data to plot, one value per node
%   p      - node positions in 3D space (nodes x 3)
%   cax    - color axis range, [min max]
%   cm     - (optional) colormap, 256 x 3. Defaults to bone(256).
%   num_ev - eigenvector/simulation number, shown in the plot title

assert( isvector(x) == 1, 'vector input required, x' )

if ( nargin > 3 ), assert( all( size(cm) == [256 3] ) )
else, cm = bone(256); end

view( [0 90] ); N = length( x );

fg1 = figure; hold on; axis image; axis off;
title(['simulation #',num2str(num_ev)]);
ax = gca;
ax.TitleHorizontalAlignment = 'left';

h = nan( N, 1 );
for ii = 1:N
	h(ii) = plot3( p(ii,1), p(ii,2), p(ii,3), 'k.', 'markersize', 50 );
end

cdata = zeros( N, 3 ); data = x;
data( data > cax(2) ) = cax(2); data( data < cax(1) ) = cax(1);

for jj = 1:N
    colorID = max( 1, sum( data(jj) > linspace(cax(1),cax(2),size(cm,1)) ) );
    cdata(jj,:) = cm( colorID, : );
    set( h(jj), 'color', cdata(jj,:) )
end
