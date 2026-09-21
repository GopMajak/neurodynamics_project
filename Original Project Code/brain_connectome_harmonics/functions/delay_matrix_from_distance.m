function tau = delay_matrix_from_distance(rr, speed, dt)
% delay_matrix_from_distance.m
%
% Converts a pairwise distance matrix (e.g. distances between parcel
% centroids or vertices) into a matrix of conduction delays, measured in
% integer simulation timesteps. Works for any distance matrix, unlike
% graphs/delay_matrix.m which only handles a 1D ring layout.
%
% Delays are floored to a minimum of 1 timestep wherever rr > 0, because
% simulate_dKM/simulate_dKM_fast need every delay to be at least 1 (a
% node can't be coupled to itself in the same timestep).
%
% Inputs: rr    : pairwise distance matrix (mm) [N x N]
%         speed : conduction velocity (m/s)
%         dt    : timestep (s)
%
% Output: tau   : delay matrix in integer timesteps [N x N]

tau = floor((rr/1000) / speed / dt); % rr (mm) -> m, then time/dt -> timesteps
has_dist = rr > 0;
tau(has_dist) = max(tau(has_dist), 1);
end
