function [A_edr_binary, A_edr_continuous, A_edr_lr, report] = build_anatomical_graphs_blocked_v2( ...
    vertices, C_norm, alpha_short, A_fit, lambda_fit, bin_edges, bin_mean, bin_std, ...
    nr_ini, nr_fin, nstd, min_lr_distance, weight_floor, block_size, seed)
%BUILD_ANATOMICAL_GRAPHS_BLOCKED_V2 Corrected EDR binary/continuous/EDR+LR
%   construction, following the faithful-replication fit validated in
%   diagnose_edr_fit_faithful.m (A=0.0658, lambda=0.1616 vs. paper's
%   0.066/0.162) and the exact LR-exception procedure read from
%   Vohryzek's own longrange_derivation_project_laplacian_v2.m:
%     - LR exceptions are detected on the GLOBAL-MAX-NORMALIZED
%       connectome (C_norm = C / max(C(:))), binned into 400 equal bins
%       over [0, max_distance], using per-bin mean+NSTD*std as the
%       threshold, restricted to bins NR_INI:NR_FIN and distance >
%       MIN_LR_DISTANCE.
%     - EDR+LR = EDR-continuous base (A_fit*exp(-lambda_fit*d))
%       everywhere EXCEPT at exception pairs, which are replaced by their
%       ACTUAL normalized connectome weight (not a fixed placeholder
%       value, correcting the earlier documented assumption in
%       build_anatomical_graphs_blocked.m).
%   EDR binary/continuous (short-range, alpha_short=0.12) are unchanged
%   from the original Stage 1 procedure -- that part was already correct
%   (Pang et al. 2023's alpha is a fixed constant, not refit from data).
%
%   NEW FILE -- does not modify build_anatomical_graphs_blocked.m; both
%   are preserved so the original (bugged-fit) and corrected results can
%   be compared side by side.

if nargin < 15
    seed = 42;
end
rng(seed);

n = size(vertices, 1);
sq_norms = sum(vertices.^2, 2);

cont_i = cell(0,1); cont_j = cell(0,1); cont_v = cell(0,1);
bin_i = cell(0,1); bin_j = cell(0,1); bin_v = cell(0,1);
lr_i = cell(0,1); lr_j = cell(0,1); lr_v = cell(0,1);

n_lr_total = 0;
block_idx = 0;

for start_row = 1:block_size:n
    block_idx = block_idx + 1;
    end_row = min(start_row + block_size - 1, n);
    rows = (start_row:end_row)';

    block_vertices = vertices(rows, :);
    d2 = sq_norms(rows) - 2*(block_vertices * vertices') + sq_norms';
    d2(d2 < 0) = 0;
    D_block = sqrt(d2);

    for r = 1:numel(rows)
        D_block(r, 1:rows(r)) = Inf; % exclude self + lower triangle
    end

    %% EDR continuous + binary (short-range, alpha=0.12, unchanged)
    p_short = exp(-alpha_short * D_block);
    keep_cont = p_short >= weight_floor;
    [rr, cc] = find(keep_cont);
    cont_i{end+1} = rows(rr); cont_j{end+1} = cc; cont_v{end+1} = p_short(keep_cont); %#ok<AGROW>

    draws = rand(size(D_block));
    keep_bin = keep_cont & (draws < p_short);
    [rr2, cc2] = find(keep_bin);
    bin_i{end+1} = rows(rr2); bin_j{end+1} = cc2; bin_v{end+1} = ones(numel(rr2), 1); %#ok<AGROW>

    %% EDR+LR: corrected fit + faithful exception detection
    edr_conn = A_fit * exp(-lambda_fit * D_block);

    C_block = C_norm(rows, :);
    bidx = discretize(D_block, bin_edges);
    valid_exc_bin = ~isnan(bidx) & bidx >= nr_ini & bidx <= nr_fin;

    threshold = inf(size(D_block));
    threshold(valid_exc_bin) = bin_mean(bidx(valid_exc_bin)) + nstd * bin_std(bidx(valid_exc_bin));
    is_exception = valid_exc_bin & (C_block > threshold) & (D_block > min_lr_distance);

    final_weight = edr_conn;
    final_weight(is_exception) = C_block(is_exception);

    keep_edrlr = is_exception | (final_weight >= weight_floor);
    [rr3, cc3] = find(keep_edrlr);
    lin_idx = sub2ind(size(D_block), rr3, cc3);
    lr_i{end+1} = rows(rr3); lr_j{end+1} = cc3; lr_v{end+1} = final_weight(lin_idx); %#ok<AGROW>

    n_lr_total = n_lr_total + nnz(is_exception);

    fprintf('  block %d: rows %d-%d done (cumulative LR exceptions: %d)\n', ...
        block_idx, start_row, end_row, n_lr_total);
end

A_edr_continuous = make_symmetric_sparse(cont_i, cont_j, cont_v, n);
A_edr_binary     = make_symmetric_sparse(bin_i, bin_j, bin_v, n);
A_edr_lr         = make_symmetric_sparse(lr_i, lr_j, lr_v, n);

n_total_pairs = n*(n-1)/2;
report.n_lr_exceptions = n_lr_total;
report.pct_lr_exceptions = 100 * n_lr_total / n_total_pairs;
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
