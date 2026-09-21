function stats = compute_distance_bin_median(vertices, W, n_bins, dist_range, block_size)
%COMPUTE_DISTANCE_BIN_MEDIAN Bin nonzero edge weights by pairwise distance
%   and compute robust per-bin statistics (median, trimmed mean) in
%   addition to the plain mean -- to check whether a small fraction of
%   anomalously strong long-range "exception" connections (the very thing
%   Vohryzek et al. 2025 define as LR exceptions) are contaminating a
%   plain per-bin MEAN and distorting the EDR exponential fit.
%
%   NEW FILE -- does not modify COMPUTE_DISTANCE_BIN_STATS.M or
%   COMPUTE_DISTANCE_BIN_STATS_V2.M; all three coexist for comparison.
%
%   Collects raw nonzero weights per bin (memory cost ~ a few GB at real
%   connectome scale -- fine on the server, not intended to run locally).

if nargin < 5
    block_size = 2000;
end

n = size(vertices, 1);
edges = linspace(dist_range(1), dist_range(2), n_bins + 1);
centers = (edges(1:end-1) + edges(2:end)) / 2;

bin_values = cell(n_bins, 1);
for b = 1:n_bins
    bin_values{b} = zeros(0, 1, 'single');
end

sq_norms = sum(vertices.^2, 2);

for start_row = 1:block_size:n
    end_row = min(start_row + block_size - 1, n);
    rows = (start_row:end_row)';

    block_vertices = vertices(rows, :);
    d2 = sq_norms(rows) - 2*(block_vertices * vertices') + sq_norms';
    d2(d2 < 0) = 0;
    D_block = sqrt(d2);

    for r = 1:numel(rows)
        D_block(r, 1:rows(r)) = -1;
    end

    W_block = W(rows, :);

    in_range = D_block >= dist_range(1) & D_block <= dist_range(2) & W_block > 0;
    bidx = discretize(D_block(in_range), edges);
    w = single(W_block(in_range));

    valid = ~isnan(bidx);
    bidx = bidx(valid);
    w = w(valid);

    for b = 1:n_bins
        bin_values{b} = [bin_values{b}; w(bidx == b)]; %#ok<AGROW>
    end
end

stats.edges = edges;
stats.centers = centers;
stats.mean_nonzero = zeros(n_bins, 1);
stats.median_nonzero = zeros(n_bins, 1);
stats.trimmed_mean_nonzero = zeros(n_bins, 1); % excludes top 5% within bin
stats.count_nonzero = zeros(n_bins, 1);

for b = 1:n_bins
    v = double(bin_values{b});
    stats.count_nonzero(b) = numel(v);
    if isempty(v)
        continue
    end
    stats.mean_nonzero(b) = mean(v);
    stats.median_nonzero(b) = median(v);
    thresh = prctile(v, 95);
    stats.trimmed_mean_nonzero(b) = mean(v(v <= thresh));
end
end
