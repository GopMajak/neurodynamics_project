function stats = compute_distance_bin_stats(vertices, W, n_bins, dist_range, block_size)
%COMPUTE_DISTANCE_BIN_STATS Bin an edge-weight matrix W by pairwise
%   Euclidean distance between rows of VERTICES.
%
%   Processes VERTICES in row-blocks so the full NxN distance matrix is
%   never materialized at once (N = 29,696 cortical vertices would need
%   ~7 GB as a single dense double matrix).
%
%   stats.edges   : n_bins+1 bin edges (mm)
%   stats.centers : n_bins bin centers (mm)
%   stats.mean    : mean W within each bin
%   stats.std     : std of W within each bin
%   stats.count   : number of (i,j) pairs, i<j, within each bin

if nargin < 5
    block_size = 2000;
end

n = size(vertices, 1);
edges = linspace(dist_range(1), dist_range(2), n_bins + 1);
centers = (edges(1:end-1) + edges(2:end)) / 2;

bin_sum = zeros(n_bins, 1);
bin_sumsq = zeros(n_bins, 1);
bin_count = zeros(n_bins, 1);

sq_norms = sum(vertices.^2, 2);

for start_row = 1:block_size:n
    end_row = min(start_row + block_size - 1, n);
    rows = (start_row:end_row)';

    block_vertices = vertices(rows, :);
    d2 = sq_norms(rows) - 2*(block_vertices * vertices') + sq_norms';
    d2(d2 < 0) = 0;
    D_block = sqrt(d2); % numel(rows) x n

    % keep only j > global row index (upper triangle, no self-pairs)
    for r = 1:numel(rows)
        D_block(r, 1:rows(r)) = -1;
    end

    W_block = W(rows, :);

    in_range = D_block >= dist_range(1) & D_block <= dist_range(2);
    bidx = discretize(D_block(in_range), edges);
    w = W_block(in_range);

    valid = ~isnan(bidx);
    bidx = bidx(valid);
    w = w(valid);

    bin_sum = bin_sum + accumarray(bidx, w, [n_bins, 1]);
    bin_sumsq = bin_sumsq + accumarray(bidx, w.^2, [n_bins, 1]);
    bin_count = bin_count + accumarray(bidx, 1, [n_bins, 1]);
end

stats.edges = edges;
stats.centers = centers;
stats.mean = bin_sum ./ max(bin_count, 1);
stats.std = sqrt(max(bin_sumsq ./ max(bin_count, 1) - stats.mean.^2, 0));
stats.count = bin_count;
end
