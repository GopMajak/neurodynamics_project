# Vohryzek et al. 2024 (PNAS) Replication — Status & Plan

Reference: J. Vohryzek, Sanz-Perl Y, Kringelbach ML, Deco G. "Human brain dynamics
are shaped by rare long-range connections over and above cortical geometry."
PNAS (2025). Code: `Reference Code Material/Vohryzeketal2024/vohryzek2024_EDRLR/`

Two distinct goals — pick the one you mean when you prompt about this later:

- **Goal A — Reproduce the published figures.** Just run the existing scripts
  against the precomputed results already bundled in the repo. No raw data
  needed. Already possible today.
- **Goal B — Full from-scratch replication.** Regenerate the intermediate
  `.mat` results yourself from raw HCP data (including last night's
  rfMRI download). Several pieces still missing/broken. See below.

## Goal A: Reproduce figures from precomputed results (ready now)

`PNAS_Figure_2_main.m` and `PNAS_Figure_3_main_part_ABC.m` both default to
`parLoop = 0` / `parfor_run = 0`, which just loads:
- `Results/long_3T/*.mat` (Figure 2)
- `Results/long_3T_task/*.mat` (Figure 3)

All of these files are present in the repo already. To actually run the
scripts you still need to:
1. Fix `repo_dir` at the top of each script (currently hardcoded to the
   original author's machine, e.g. `/project_laplacian/BrainEigenmodes_EDRLR-main`).
2. Have MATLAB (tested on R2023b) with the helper functions on path
   (`calc_network_eigenmode`, `calc_parcellate`, `calc_eigendecomposition`,
   `calc_triu_ind`, `read_vtk`, `daviolinplot`, `shadedErrorBar`, `redblue`,
   `draw_surface_*` — check these all resolve; some come from the Pang et al.
   BrainEigenmodes repo, already present under `Reference Code Material/pang2023_BrainEigenmodes/`).

## Goal B: Full from-scratch replication — inventory

### Already have
- 255/255 subjects' `rfMRI_REST1_LR` raw CIFTI (`*.dtseries.nii`) — downloaded,
  verified complete (`data/empirical/rfMRI_raw/`).
- Precomputed synthetic eigenmode connectomes: EDR binary, EDR continuous,
  Geometry, EDR+LR (`Connectome_derivation/Connectomes/synthetic_*_eigenmodes_fsLR_32k-lh_200.mat`).
- Structural connectome: `Data/empirical/S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat`.
- Task fMRI z-stats: `Data/empirical/S255_tfMRI_ALLTASKS_raw_lh.mat`.
- Surfaces (`fsLR_32k_midthickness-{lh,rh}.vtk`), cortex masks, Glasser360
  parcellation, `subject_list_HCP.txt` (255 subject IDs).

### Missing / blocking
1. **Per-subject hemisphere-extracted timeseries.** Reconstruction scripts
   (`PNAS_Figure_2_main.m` parLoop=1 branch, `FC_longrange_derivation.m`)
   expect `subject_<ID>_rfMRI_timeseries-lh.mat` (a `timeseries` variable,
   left-cortex vertices x time), not the raw whole-brain CIFTI. **No
   conversion script exists yet in this repo** — needs to be written
   (load CIFTI → apply `fsLR_32k_cortex-lh_mask.txt` → save per subject).
2. **`synthetic_distance_connectome.mat`** (variable `surface_dist`, pairwise
   Euclidean distance between cortical vertices) — referenced by
   `longrange_derivation_project_laplacian_v2.m` but not found anywhere in
   the project. Can likely be recomputed directly from the surface
   (`pdist` on `surface_midthickness.vertices(cortex_ind,:)`, as done inline
   in `EDR_connectomes_derivation.m`).
3. **`Euclidean_distance_LRE.mat` / `EDR_LR_derivation_v2.mat`** (variable
   `EDR_LRE`) — referenced by `PNAS_Figure_2_main.m` and
   `EDRLR_connectome_derivation.m`, but not present. `FC_longrange_derivation.m`
   *computes* `EDR_LRE` in memory but never saves it — that save step is
   missing from the script as given and needs to be added.
4. **Filename mismatches between derivation output and figure-script input.**
   `EDR_connectomes_derivation.m` saves e.g.
   `synthetic_EDRconnectome_binary_eigenmodes-lh_200.mat`; figure scripts look
   for `synthetic_EDRbinary_eigenmodes-fsLR_32k-lh_200.mat`. Filenames need to
   be reconciled (rename on save, or on load) before the pipeline chains
   correctly.
5. **Hardcoded absolute paths** throughout every script (`/project_laplacian/...`,
   `/Users/jakub/Datasets/HCP/HCP255/...`) need to be repointed at this
   project's directories.
6. **OSF repository** linked in the README
   (`Reference Code Material/Vohryzeketal2024/vohryzek2024_EDRLR/README.md`) —
   the link itself is malformed there (display text and href point to two
   different OSF IDs: `asntf` vs `3qjp5`). Worth checking that OSF page
   directly — it may already contain the missing distance/derivation `.mat`
   files, making steps 2–3 unnecessary.

## Suggested order of operations for Goal B

1. Resolve the OSF link and check what it actually contains (may shortcut
   steps below).
2. Write the CIFTI → per-subject left-hemisphere timeseries extraction
   script (blocker for everything downstream that touches raw fMRI).
3. Compute/recover `surface_dist` (`synthetic_distance_connectome.mat`).
4. Run `FC_longrange_derivation.m`, adding a `save(...)` for `EDR_LRE`.
5. Run `EDRLR_connectome_derivation.m` and `EDR_connectomes_derivation.m`,
   fixing output filenames to match what the figure scripts expect (or fix
   the figure scripts' load paths — pick one convention and apply
   consistently).
6. Fix hardcoded `repo_dir` / dataset paths in every script to point at this
   project.
7. Run `PNAS_Figure_2_main.m` and `PNAS_Figure_3_main_part_ABC.m` with
   `parLoop`/`parfor_run = 1` to regenerate `Results/*.mat` from scratch, then
   compare against the bundled precomputed versions as a sanity check.

## Notes
- Only `REST1_LR` (single run, left hemisphere only) is used — consistent
  with the Pang et al. 2023 BrainEigenmodes convention this codebase builds on.
- 255 is the fixed cohort size throughout (matches `subject_list_HCP.txt`).
