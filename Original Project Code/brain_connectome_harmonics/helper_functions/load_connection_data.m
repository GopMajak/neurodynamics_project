function [N,w,p,d,tau,frontal_idx] = load_connection_data( filepath, v, dt )
% Loads an HCP structural connectivity dataset and builds the matrices
% needed for simulation: connection weights, node positions, distances,
% and conduction delays. Also finds which nodes are frontal-lobe regions.
%
% Inputs:
%   filepath - path to the .mat connectivity data file
%   v        - conduction speed (m/s)
%   dt       - simulation timestep (s)
%
% Outputs:
%   N           - number of nodes
%   w           - weight (connectivity) matrix
%   p           - node positions in 3D space
%   d           - distance matrix between nodes
%   tau         - delay matrix, in timesteps
%   frontal_idx - indices of nodes in frontal-lobe regions

C = load( filepath );

N = size( C.CIJ_fbden_average, 1 );

% weight matrix: zero entries mean "not connected", not "zero-weight",
% so they're set to NaN instead of 0
w = C.CIJ_fbden_average; w(w==0) = NaN;

p = C.roi_xyz_avg';

d = C.CIJ_edgelength_average;

% delay = distance (converted from mm to m) / conduction speed,
% then converted from seconds to number of timesteps
tau = round( (d./1000) ./ (v*dt) );

% Frontal ROIs = areas anterior to the precentral gyrus (primary motor
% area), following Hagmann et al. 2008. Included here because Andrillon
% et al. 2011 found a sharp frequency transition across the supplemental
% motor area, so it's useful to be able to pull out just the frontal nodes.
frontal_lbls = { 'rSF'; 'rFP'; 'rMOF'; 'rRAC'; 'rCAC'; 'rCMF'; 'rRMF'; 'rPOPE'; 'rPTRI'; 'rPORB'...
    ; 'lSF'; 'lFP'; 'lMOF'; 'lRAC'; 'lCAC'; 'lCMF'; 'lRMF'; 'lPOPE'; 'lPTRI'; 'lPORB'; ...
    'lLOF'; 'rLOF' };

anat = strtrim( cellstr(C.anat_lbls) );

match = cellfun( @(x) ismember(x, frontal_lbls), anat, 'UniformOutput', 0 );
r = find( cell2mat(match) );

frontal_idx = find( ismember(C.roi_lbls, r) );
