%% Build the missing Euclidean_distance_LRE.mat
%
% PNAS_Figure_2_main.m (original Vohryzek et al. 2025 repo) unconditionally
% loads 'Connectome_derivation/Euclidean_distance_LRE.mat' before its
% parLoop branch, but this file is not included in the repo/OSF bundle
% we have. Its only used variable in the parLoop=0 path (the path we run,
% since we're using the precomputed Results/*.mat) is dead weight -- it's
% never referenced again in that branch -- but the script will still
% crash on the load if the file doesn't exist.
%
% This script computes the one thing that IS a well-defined, deterministic
% quantity derivable from data we already have: FC_LRE_euc_dist_matrix,
% the parcellated (Glasser360, 180x180) pairwise Euclidean distance matrix
% between left-cortex vertices, exactly as done inline in the authors' own
% FC_longrange_derivation.m (surface_dist = squareform(pdist(...)); then
% parcellated). No fMRI-derived quantity (e.g. FC_emp_HCP_255_subject_mean)
% is fabricated here since it is unused in the parLoop=0 reproduction path.

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(this_dir)); % Code/reproduction -> Code -> project root
addpath(fullfile(project_root, 'Code', 'utils'));

vohryzek_dir = fullfile(project_root, 'Reference Code Material', 'Vohryzeketal2024', 'vohryzek2024_EDRLR');

surface_file = fullfile(vohryzek_dir, 'Data', 'template_surfaces', 'fsLR_32k_midthickness-lh.vtk');
cortex_mask_file = fullfile(vohryzek_dir, 'Data', 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt');
parc_file = fullfile(vohryzek_dir, 'Data', 'parcellations', 'fsLR_32k_Glasser360-lh.txt');

out_dir = fullfile(vohryzek_dir, 'Connectome_derivation');
out_file = fullfile(out_dir, 'Euclidean_distance_LRE.mat');

fprintf('Loading surface + cortex mask + parcellation...\n');
[vertices, ~] = read_vtk_surface(surface_file);
cortex = dlmread(cortex_mask_file); %#ok<DLMRD>
cortex_ind = find(cortex);
parc = dlmread(parc_file); %#ok<DLMRD>

fprintf('Computing pairwise Euclidean distance for %d cortical vertices...\n', numel(cortex_ind));
surface_dist = squareform(pdist(vertices(cortex_ind, :)));

fprintf('Parcellating to Glasser360 (%d parcels)...\n', numel(unique(parc(parc>0))));
FC_LRE_euc_dist_matrix = parcellate_matrix(surface_dist, parc(cortex_ind)); %#ok<NASGU>

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save(out_file, 'FC_LRE_euc_dist_matrix');
fprintf('Saved -> %s\n', out_file);
