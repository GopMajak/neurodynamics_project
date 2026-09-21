function stats = compute_distance_bin_stats_v2(vertices, W, n_bins, dist_range, block_size)
%COMPUTE_DISTANCE_BIN_STATS_V2 Like COMPUTE_DISTANCE_BIN_STATS, but also
%   tracks nonzero-only statistics per bin, to diagnose whether the
%   zero-inflation typical of tractography-derived structural connectomes
%   (most vertex pairs have zero detected streamlines) distorts the EDR
%   exponential fit when zeros are included in the per-bin mean.
%
%   NEW FILE -- does not modify or replace COMPUTE_DISTANCE_BIN_STATS.M
%   (Stage 1 original), so both are preserved for side-by-side comparison.
%
%   stats.mean_all     : mean W within each bin, including zero entries
%   stats.mean_nonzero  : mean W within each bin, zero entries excluded
%   stats.frac_nonzero  : fraction of pairs in each bin with W > 0
%   stats.std_all / std_nonzero : corresponding standard deviations

if nargin < 5
    block_size = 2000;
end

n = size(vertices, 1);
edges = linspace(dist_range(1), dist_range(2), n_bins + 1);
centers = (edges(1:end-1) + edges(2:end)) / 2;

sum_all = zeros(n_bins, 1);
sumsq_all = zeros(n_bins, 1);
count_all = zeros(n_bins, 1);

sum_nz = zeros(n_bins, 1);
sumsq_nz = zeros(n_bins, 1);
count_nz = zeros(n_bins, 1);

sq_norms = sum(vertices.^2, 2);

for start_row = 1:block_size:n
    end_row = min(start_row + block_size - 1, n);
    rows = (start_row:end_row)';

    block_vertices = vertices(rows, :);
    d2 = sq_norms(rows) - 2*(block_vertices * vertices') + sq_norms';
    d2(d2 < 0) = 0;
    D_block = sqrt(d2);

    for r = 1:numel(rows)
        D_block(r, 1:rows(r)) = -1; % exclude self + lower triangle
    end

    W_block = W(rows, :);

    in_range = D_block >= dist_range(1) & D_block <= dist_range(2);
    bidx = discretize(D_block(in_range), edges);
    w = W_block(in_range);

    valid = ~isnan(bidx);
    bidx = bidx(valid);
    w = w(valid);

    % all-pairs stats (zeros included)
    sum_all = sum_all + accumarray(bidx, w, [n_bins, 1]);
    sumsq_all = sumsq_all + accumarray(bidx, w.^2, [n_bins, 1]);
    count_all = count_all + accumarray(bidx, 1, [n_bins, 1]);

    % nonzero-only stats
    is_nz = w > 0;
    sum_nz = sum_nz + accumarray(bidx(is_nz), w(is_nz), [n_bins, 1]);
    sumsq_nz = sumsq_nz + accumarray(bidx(is_nz), w(is_nz).^2, [n_bins, 1]);
    count_nz = count_nz + accumarray(bidx(is_nz), 1, [n_bins, 1]);
end

stats.edges = edges;
stats.centers = centers;

stats.mean_all = sum_all ./ max(count_all, 1);
stats.std_all = sqrt(max(sumsq_all ./ max(count_all,1) - stats.mean_all.^2, 0));
stats.count_all = count_all;

stats.mean_nonzero = sum_nz ./ max(count_nz, 1);
stats.std_nonzero = sqrt(max(sumsq_nz ./ max(count_nz,1) - stats.mean_nonzero.^2, 0));
stats.count_nonzero = count_nz;

stats.frac_nonzero = count_nz ./ max(count_all, 1);
end
