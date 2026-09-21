function [eig_vec, eig_val] = calc_network_eigenmode_lowmem(network, num_modes)
% calc_network_eigenmode_lowmem.m
%
% Computes the same eigenmodes as calc_network_eigenmode_dense.m (the
% normalized graph Laplacian of `network`), but without ever building the
% full N x N Laplacian matrix. At full connectome resolution (~30,000
% nodes) that matrix would take several GB of memory, and a dense eig()
% on top of it is enough to crash MATLAB -- even though we only actually
% need a couple hundred of the possible eigenpairs.
%
% Instead this uses eigs(), MATLAB's iterative solver for "give me just
% the top/bottom K eigenpairs", with a matrix-free function handle
% (Bmult) so it only ever multiplies by `network` itself, never forms the
% full Laplacian. It solves for the largest eigenvalues of a related
% matrix M instead of the smallest eigenvalues of the Laplacian, because
% the two are mathematically related (Laplacian eigenvalue = 1 - M's
% eigenvalue) and asking eigs() for the largest ones is much cheaper.
%
% Use calc_network_eigenmode_dense.m instead for small networks (e.g. a
% parcellated connectome) -- it's simpler and just as fast there.
%
% Inputs: network   : symmetric connectivity matrix [N x N]
%         num_modes : number of modes to return (int)
%
% Outputs: eig_vec  : eigenvectors (eigenmodes) [N x num_modes]
%          eig_val  : eigenvalues [num_modes x 1]

if nargin < 2
    num_modes = size(network,1);
end

N = size(network,1);

% Remove self-connections
network(1:(N+1):end) = 0;

d = double(sum(network, 2));
dhalf = 1 ./ sqrt(d);

% isolated nodes (no connections, d=0) need special handling so they
% don't produce 0/0 = NaN -- treat them as contributing their own
% isolated eigenvalue instead
iso = (d == 0);
dhalf(iso) = 0;
iso_mask = double(iso);

% network is stored in single precision to save memory; math is done in
% double precision for accuracy
Bmult = @(x) iso_mask.*x + dhalf .* double(network * single(dhalf .* x));

opts.issym = true;
opts.isreal = true;
opts.p = min(N, 4*num_modes);   % wider search space = more reliable convergence
opts.maxit = 1000;

[eig_vec, eig_val_M] = eigs(Bmult, N, num_modes, 'largestreal', opts);

eig_val_all = 1 - diag(eig_val_M);
[eig_val, sort_idx] = sort(eig_val_all, 'ascend');
eig_vec = single(eig_vec(:, sort_idx));
