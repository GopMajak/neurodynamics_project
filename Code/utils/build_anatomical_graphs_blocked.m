function [A_edr_binary, A_edr_continuous, A_edr_lr, report] = build_anatomical_graphs_blocked( ...
    vertices, SC, alpha_short, alpha_fit, lr_stats, sd_thresh, min_lr_distance, weight_floor, block_size, seed)
%BUILD_ANATOMICAL_GRAPHS_BLOCKED Construct the EDR binary, EDR continuous,
%   and EDR+LR vertex-level adjacency matrices in a single memory-safe
%   blocked pass over pairwise distances.
%
%   vertices        : cortical vertex coordinates (N x 3, mm)
%   SC              : N x N structural connectome (avgSC_L)
%   alpha_short     : EDR decay for standalone EDR binary/continuous
%                     (0.12 mm^-1, matches Pang et al. 2023)
%   alpha_fit       : EDR decay for the base of EDR+LR, fitted to this
%                     connectome (~0.162 mm^-1, Vohryzek et al. 2025)
%   lr_stats        : distance-bin mean/std of SC (COMPUTE_DISTANCE_BIN_STATS),
%                     used to flag long-range exceptions
%   sd_thresh       : SD threshold above the per-bin mean (3)
%   min_lr_distance : minimum distance (mm) to qualify as "long-range" (40)
%   weight_floor    : minimum EDR weight to keep an edge (sparsification)
%   block_size      : rows processed per block
%   seed            : RNG seed for the stochastic EDR binary draw
%
%   EDR+LR construction (documented modeling choice -- the papers describe
%   the ingredients but not an explicit combination formula): edges are
%   the union of the fitted EDR-continuous base (weight = exp(-alpha_fit*d))
%   and the detected LR exceptions, which are inserted at weight 1 (i.e.
%   treated as being as strong as a direct/local connection), consistent
%   with Fig. 1A-C of Vohryzek et al. 2025 depicting LR exceptions as
%   additional discrete edges layered on top of the local EDR mesh.

if nargin < 10
    seed = 42;
end
rng(seed);

n = size(vertices, 1);
sq_norms = sum(vertices.^2, 2);

cont_i = cell(0, 1); cont_j = cell(0, 1); cont_v = cell(0, 1);
bin_i = cell(0, 1); bin_j = cell(0, 1); bin_v = cell(0, 1);
lr_i = cell(0, 1); lr_j = cell(0, 1); lr_v = cell(0, 1);

n_lr_total = 0;
block_idx = 0;

for start_row = 1:block_size:n
    block_idx = block_idx + 1;
    end_row = min(start_row + block_size - 1, n);
    rows = (start_row:end_row)';

    block_vertices = vertices(rows, :);
    d2 = sq_norms(rows) - 2*(block_vertices * vertices') + sq_norms';
    d2(d2 < 0) = 0;
    D_block = sqrt(d2); % numel(rows) x n

    for r = 1:numel(rows)
        D_block(r, 1:rows(r)) = Inf; % exclude self + lower triangle
    end

    %% EDR continuous + binary (short-range alpha, Pang-consistent)
    p_short = exp(-alpha_short * D_block);
    keep_cont = p_short >= weight_floor;
    [rr, cc] = find(keep_cont);
    cont_i{end+1} = rows(rr); cont_j{end+1} = cc; cont_v{end+1} = p_short(keep_cont); %#ok<AGROW>

    draws = rand(size(D_block));
    keep_bin = keep_cont & (draws < p_short);
    [rr2, cc2] = find(keep_bin);
    bin_i{end+1} = rows(rr2); bin_j{end+1} = cc2; bin_v{end+1} = ones(numel(rr2), 1); %#ok<AGROW>

    %% EDR+LR base (fitted alpha) + long-range exceptions from SC
    p_fit = exp(-alpha_fit * D_block);
    keep_fit = p_fit >= weight_floor;

    W_block = SC(rows, :);
    bidx = discretize(D_block, lr_stats.edges);
    valid_bin = ~isnan(bidx);
    thresh = inf(size(D_block));
    thresh(valid_bin) = lr_stats.mean(bidx(valid_bin)) + sd_thresh * lr_stats.std(bidx(valid_bin));
    is_lr = valid_bin & (D_block > min_lr_distance) & (W_block > thresh);

    keep_edrlr = keep_fit | is_lr;
    edrlr_weight = p_fit;
    edrlr_weight(is_lr) = 1;

    [rr3, cc3] = find(keep_edrlr);
    lin_idx = sub2ind(size(D_block), rr3, cc3);
    lr_i{end+1} = rows(rr3); lr_j{end+1} = cc3; lr_v{end+1} = edrlr_weight(lin_idx); %#ok<AGROW>

    n_lr_total = n_lr_total + nnz(is_lr);

    fprintf('  block %d: rows %d-%d done (cumulative LR exceptions: %d)\n', ...
        block_idx, start_row, end_row, n_lr_total);
end

A_edr_continuous = make_symmetric_sparse(cont_i, cont_j, cont_v, n);
A_edr_binary     = make_symmetric_sparse(bin_i, bin_j, bin_v, n);
A_edr_lr         = make_symmetric_sparse(lr_i, lr_j, lr_v, n);

report.n_lr_exceptions = n_lr_total;
report.density_binary = nnz(A_edr_binary) / (n*(n-1));
report.density_continuous = nnz(A_edr_continuous) / (n*(n-1));
report.density_edrlr = nnz(A_edr_lr) / (n*(n-1));
end

function A = make_symmetric_sparse(ic, jc, vc, n)
i = cat(1, ic{:});
j = cat(1, jc{:});
v = cat(1, vc{:});
A = sparse([i; j], [j; i], [v; v], n, n);
end
