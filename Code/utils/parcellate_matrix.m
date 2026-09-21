function parc_matrix = parcellate_matrix(matrix, parc_labels)
%PARCELLATE_MATRIX Average an NxN vertex-level matrix into parcel x parcel
%   blocks, by taking the mean of every (vertex-in-p, vertex-in-q) entry.
%   PARC_LABELS (vertices x 1; label 0 excluded) must index MATRIX's rows
%   and columns.

parc_labels = parc_labels(:);
labels = unique(parc_labels(parc_labels > 0));
num_parcels = numel(labels);
parc_matrix = zeros(num_parcels, num_parcels);

idx_cell = cell(num_parcels, 1);
for p = 1:num_parcels
    idx_cell{p} = find(parc_labels == labels(p));
end

for p = 1:num_parcels
    for q = p:num_parcels
        block = matrix(idx_cell{p}, idx_cell{q});
        val = mean(block(:));
        parc_matrix(p, q) = val;
        parc_matrix(q, p) = val;
    end
end
end
