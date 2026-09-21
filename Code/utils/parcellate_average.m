function parc_data = parcellate_average(data, parc_labels)
%PARCELLATE_AVERAGE Average vertex-wise data within parcels.
%   PARC_DATA = PARCELLATE_AVERAGE(DATA, PARC_LABELS) averages rows of
%   DATA (vertices x features) within each parcel defined by PARC_LABELS
%   (vertices x 1; label 0 = unlabeled/medial wall, excluded). Returns
%   PARC_DATA as (num_parcels x features), ordered by increasing label.
%   DATA and PARC_LABELS must be indexed over the same vertex set.

parc_labels = parc_labels(:);
labels = unique(parc_labels(parc_labels > 0));
num_parcels = numel(labels);
parc_data = zeros(num_parcels, size(data, 2));

for p = 1:num_parcels
    idx = parc_labels == labels(p);
    parc_data(p, :) = mean(data(idx, :), 1);
end
end
