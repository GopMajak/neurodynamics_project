function matrix_parcellated = calc_parcellate_matrix(parc, matrix_input)
% calc_parcellate_matrix.m
%
% Shrinks a vertex-by-vertex matrix (e.g. a connectome or distance matrix)
% down to a parcel-by-parcel matrix, by averaging all the vertex entries
% that fall inside each pair of parcels. Vertices labeled parc==0 (not
% assigned to any parcel) are left out.
%
% Written using sparse matrix multiplication instead of nested loops so
% it stays fast even at full-brain vertex resolution.
%
% Inputs: parc         : parcellation labels aligned to matrix rows/cols [Nx1]
%         matrix_input : vertex-level matrix [NxN]
%
% Output: matrix_parcellated : parcellated matrix [num_parcels x num_parcels]

parcels = unique(parc(parc>0));
num_parcels = length(parcels);
num_vertices = size(matrix_input, 1);

parcel_of_vertex = zeros(num_vertices, 1);
for p = 1:num_parcels
    parcel_of_vertex(parc==parcels(p)) = p;
end
valid = parcel_of_vertex > 0;

% P(p,v) = 1 if vertex v belongs to parcel p, else 0
P = sparse(parcel_of_vertex(valid), find(valid), 1, num_parcels, num_vertices);

counts = full(sum(P, 2)); % number of vertices in each parcel
sum_matrix = P * double(matrix_input) * P';
count_matrix = counts * counts';

matrix_parcellated = full(sum_matrix) ./ count_matrix;
