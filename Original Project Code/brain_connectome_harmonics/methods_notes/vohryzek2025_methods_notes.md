# Methods notes: Vohryzek et al. 2025 (PNAS), "Human brain dynamics are shaped
# by rare long-range connections over and above cortical geometry"

This paper's Methods section maps almost one-to-one onto my code
(`connectome_harmonics.m`, `step4_5_full_resolution_reconstruction.m`,
`step5_static_reconstruction.m`, the `figure_*` scripts), so every subsection
here is directly relevant — unlike the Pang notes, nothing gets skipped.

---

## 1. Data sources (all borrowed from Pang et al. 2023)

> "We used a subset of 255 participants... consistent with work of Pang et
> al." / "...provided by work of Pang et al. (5) from https://osf.io/xczmp/
> in 'S255_tfMRI_ALLTASKS_raw_lh.mat'..."

**In plain terms:** This paper doesn't collect new brain-scan data — it
reuses the exact same 255-person HCP dataset, connectome, and cortical
surface mesh that Pang et al. published, so results between the two papers
are directly comparable.

**My thoughts:** This explains why my `eigenmode_toolbox_dir` (the
`pang2023_BrainEigenmodes` folder) shows up constantly even in scripts that
are "about" Vohryzek's method — the task-fMRI data, the empirical
connectome, and the surface mesh all physically live in Pang's data release,
not Vohryzek's. Vohryzek's own folder (`edrlr_data_dir`) mainly contributes
the *precomputed results files* (`Results/long_3T/*.mat`,
`Results/long_3T_task/*.mat`) and templates/parcellations — worth remembering
this split when tracing where a `load(...)` call's file actually comes from.

---

## 2. The EDR (exponential distance rule) fit

> "C_EDR(i,j) = A·e^(-λ·r(i,j))... we have generated 400 bins of equal
> Euclidean distance taking the bins spanning 10 to 170 mm (thus excluding
> the first 25 bins)... The estimation yielded A = 0.066 and λ = 0.162 mm⁻¹."

**In plain terms:** Same idea as in Pang et al.: connection strength decays
roughly exponentially with distance. Here they bin all vertex-pairs by
distance (400 bins), throw out the bins for very short distances (too noisy/
dominated by direct neighbors), and fit a curve to what's left.

**My thoughts:** This is almost exactly `connectome_harmonics.m` Stage 2:
`NR = 400`, `fit_ind = 25:NR`. The paper's prose ("excluding the first 25
bins") reads like the fit should start at bin 26, which looked like a
possible off-by-one against my `fit_ind = 25:NR` (starts at bin 25) — but I
checked the original author's own script
(`Code from Vohryzek et al 2024/vohryzek2024_EDRLR/Connectome_derivation/
longrange_derivation_project_laplacian_v2.m`, line 81) and it fits with
`xcoor(25:end)` — i.e. bin 25 onward, exactly matching my code. So the paper
text is just a slightly loose description of the same thing; no mismatch,
no fix needed here. My own fitted A/λ still won't exactly equal their
reported 0.066/0.162 (different data sample, `fminsearch` vs. `lsqcurvefit`),
but the *procedure* is confirmed identical.

---

## 3. Long-range (LR) exceptions

> "We defined connectivity exceptions as three SD above the mean for a given
> distance bin that are longer than 40 mm."

**In plain terms:** Most connections follow the smooth exponential falloff
above. A small number don't — they're much stronger than the falloff
predicts for their distance. If a connection is both (a) unusually strong for
its distance bin (>3 standard deviations above that bin's mean) and (b)
actually far away (>40mm), it counts as a "long-range exception."

**My thoughts:** This is `NSTD = 3` and `DistRange = 40` in every one of my
scripts that builds an EDR+LR connectome — direct 1:1 match, no ambiguity
here. This is also the exact quantity my `diag_lr_magnitude.m` diagnostic
script investigates: it doesn't just count these exceptions, it asks how
*big* they actually are relative to the smooth baseline, which matters
because "statistically rare" (>3 SD) doesn't automatically mean "practically
large."

---

## 4. Four connectome/mode types: Geometry, EDR binary, EDR continuous, EDR+LR

> "EDR binary... a binary adjacency matrix... EDR continuous... all the
> connections and their weights are kept... EDR+LR... we combine the EDR
> continuous with LR exceptions to the EDR."

**In plain terms:** The paper compares four different ways of building the
network whose eigenmodes get used for reconstruction:
- **Geometry** — ignores the connectome entirely, uses pure cortex shape (Pang's method).
- **EDR binary** — a yes/no synthetic network following the smooth exponential rule (thresholded to be 0/1).
- **EDR continuous** — same smooth exponential rule, but keeping the actual weighted strengths.
- **EDR+LR** — EDR continuous, but with the rare long-range exceptions spliced back in from the real data.

Comparing all four isolates exactly what the long-range connections add on
top of pure geometry and on top of a smooth wiring rule.

**My thoughts:** This is precisely the `graph_names = {'Geometry', 'EDR
binary', 'EDR continuous', 'EDR+LR'}` list built in both `step4_5` and
`step5`. Useful to keep in mind while reading my own code: "EDR binary" and
"EDR continuous" always use the fixed λ=0.12 (Pang's number, see the Pang
notes), while "EDR+LR" uses *my own freshly-fit* λ for its smooth baseline
before splicing in the LR exceptions — these are deliberately *not* the same
λ, and conflating them would silently erase the comparison the whole analysis
depends on. My code comments already flag this, but it's worth remembering
*why* the paper needs both: EDR binary/continuous exist to replicate Pang's
own comparison conditions, while EDR+LR's baseline needs to be a fair,
freshly-fit baseline specific to the long-range-splice question.

---

## 5. Laplacian decomposition (getting the harmonic modes)

> "L_norm = D^(-1/2) L D^(-1/2), with L = D - A... the harmonic modes were
> computed as eigenvectors of Δ_A ψ_k(x_i) = λ_k ψ_k(x_i)."

**In plain terms:** Once you have any of the four adjacency matrices above,
you get its "modes" the same way as the connectome eigenmodes in Pang et
al. — normalize the graph Laplacian and take its eigenvectors.

**My thoughts:** Same underlying math as Pang's connectome-eigenmode section
(see Pang notes, section 5) — my `calc_network_eigenmode_dense.m` /
`calc_network_eigenmode_lowmem.m` serve both papers' methods with no changes
needed between them, since the normalized-Laplacian eigenproblem is
identical regardless of which of the four adjacency matrices goes in.

---

## 6. Reconstructing brain activity & measuring error

> "For the spontaneous fMRI, we... focused on... reconstructing the LR FC
> derived as a subset of connections with high-correlation values (>0.5
> correlation) and a long Euclidean distance (>40mm)... For the task-based
> fMRI, we calculated the reconstruction error as the mse distance between
> the empirical and reconstructed task-activation maps."

**In plain terms:** Same reconstruct-with-N-modes idea as Pang et al., but
two different ways of scoring how good the reconstruction is depending on
data type: for resting-state functional connectivity, they zoom in
specifically on the long-range functional connections (the ones that are
both strongly correlated *and* far apart) rather than the whole FC matrix;
for task maps, they just use plain mean-squared-error (MSE) between the real
and reconstructed activation map.

**My thoughts:** This explains a detail I hadn't fully connected before:
`figure_2_fMRI_reconstruction.m` loads `recon_mse_parc_SC_exceptions_version`
— the "SC_exceptions" naming refers to exactly this >0.5-correlation-and->40mm
long-range functional-connection restriction, applied by the original authors
when they generated these precomputed files. My script never applies that
filter itself; it's already baked into the `.mat` file I'm loading. Good to
note in my thesis methods so it's clear that filter isn't something I
implemented — I'm just consuming their precomputed result. By contrast, my
`step4_5`/`step5`/`figure_3_tfMRI_reconstruction.m` task-fMRI reconstructions
use plain `corr`/`immse` with no such distance/correlation restriction,
matching the task-data half of this same Methods paragraph.

---

## Summary: things worth double-checking or writing up explicitly

1. ~~The `fit_ind = 25:NR` vs. "excluding the first 25 bins" off-by-one~~ —
   checked against the original author's script, confirmed `fit_ind = 25:NR`
   matches exactly (section 2). No action needed.
2. Why my code deliberately keeps λ=0.12 (EDR binary/continuous) and my own
   fitted λ (EDR+LR baseline) as two *different* numbers, never conflated
   (section 4).
3. The long-range functional-connection restriction (>0.5 corr, >40mm) is
   baked into Vohryzek's precomputed resting-state MSE files, not something
   my own code applies (section 6) — only the task-fMRI comparisons in my
   code compute reconstruction error from scratch.
