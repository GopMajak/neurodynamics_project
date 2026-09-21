# Methods notes: Pang et al. 2023 (Nature), "Geometric constraints on human brain function"

Reading these notes as a grad student trying to actually replicate this paper's
pipeline (not just cite it) — for each Methods subsection I pulled out what
matters, restated it simply, and added my own notes on how it connects to my
code in this project.

Only the subsections that this project's code actually implements or depends
on are covered here. Pang et al. also cover an NFT wave model, a biophysical
neural-mass (BEI) model, resting-state dynamics fitting, time-lagged
(lag-thread) analysis, and eigenmodes of subcortical/hippocampal structures —
none of that is used anywhere in this project, so I skipped it.

---

## 1. Deriving the geometric eigenmodes

> "We employed the LaPy python library... to derive the geometric eigenmodes
> of the human cortex... a triangular surface mesh... comprising 32,492
> vertices in each hemisphere... We used the first 200 modes."

**In plain terms:** They take the cortex surface (like a bumpy 3D sheet) and
solve a math problem called the Laplace-Beltrami eigenvalue problem on it —
basically the same math that tells you the resonant vibration patterns of a
drum, but for the shape of the brain instead of a flat drum skin. Each
solution ("eigenmode") is a wave-like pattern spread across the whole cortex.
Mode 1 is just a constant (no wiggles), and later modes have more and more
fine wiggles. They keep the first 200 modes because adding more barely
improves how well you can reconstruct real brain activity.

**My thoughts:** I never actually run this step myself — this is exactly what
Pang et al.'s precomputed geometric eigenmode files give me for free
(`data/template_eigenmodes/..._emode_200.txt` and
`data/results/basis_geometric_..._evec_200.txt`, loaded directly in
`step4_5_full_resolution_reconstruction.m` and `step5_static_reconstruction.m`).
Good to know I don't need LaPy or the finite-element solve myself — that's a
big chunk of implementation complexity I get to skip because they released
the outputs, not just the code.

---

## 2. Reconstructing brain activity from eigenmodes

> "yₜ(r,t) = Σⱼ aⱼ(t) ψⱼ(r)... the amplitudes can be obtained by integrating
> over the cortical surface... aⱼ(t) = ∫ y(r,t) ψⱼ(r) dr."

**In plain terms:** Once you have a set of eigenmodes (a "basis," like a set
of building blocks), you can approximate any real brain activity map as a
weighted sum of those modes. The "weight" for each mode is found the same way
you'd project a vector onto an axis — an inner product / dot product between
the real data and that mode's pattern. Using more modes in the sum gives you
a better (but never perfect) approximation of the real data.

**My thoughts:** This is the core operation my `step4_5` and `step5` scripts
run over and over: project real task-fMRI data onto the first *N* eigenmodes
(`calc_eigendecomposition(..., 'matrix')` — one of Pang's own helper
functions, not something I wrote), then reconstruct with `basis * beta` and
compare to the real map. Worth remembering: this is a **least-squares
projection**, not a Fourier transform, even though the mental picture ("sum
of modes with different wavelengths") is very Fourier-like — that distinction
matters if I ever try to explain why this isn't literally spatial-frequency
filtering.

---

## 3. HCP data used

> "We analysed data from 255 unrelated healthy individuals... task-evoked
> fMRI measured in seven task domains... and task-free resting-state fMRI...
> mapped onto the fsLR-32k CIFTI space with 32,492 vertices in each
> hemisphere."

**In plain terms:** All the brain-activity data in this paper comes from the
Human Connectome Project (HCP): 255 people scanned doing 7 different tasks
plus resting-state scans, all resampled onto the same standard cortical
surface mesh so that "vertex #12345" means the same physical location in
everyone's brain.

**My thoughts:** This is exactly the `S255_tfMRI_ALLTASKS_raw_lh.mat` file
(47 task contrasts, not 7 — the 7 *domains* expand into 47 individual
contrasts, e.g. "motor" splits into several specific movements) that
`step4_5`, `step5`, and `figure_3D_surface_reconstruction.m` all load. Good
to keep straight: "7 tasks" (paper's headline number) and "47 contrasts"
(what's actually in the file and what my scripts loop over) are both correct,
just different levels of granularity.

---

## 4. Parcellation used

> "We presented results using the HCP-MMP1 parcellation with 180 regions per
> hemisphere (we term this the Glasser360 parcellation)."

**In plain terms:** Instead of working vertex-by-vertex (32,492 of them),
results are often summarized by averaging within 180 named anatomical
regions per hemisphere — a coarser, more interpretable resolution. "Glasser360"
just means 180 regions x 2 hemispheres = 360 total, even though this project
only ever uses the left hemisphere's 180.

**My thoughts:** This is the `parc_name = 'Glasser360'` used everywhere in my
scripts for evaluating reconstruction error, and separately (at parcel
resolution) for the Stage 1 pilot EDR fit in `connectome_harmonics.m` and
`step5_static_reconstruction.m`. Worth remembering *why* I parcellate before
measuring error even when eigenmodes are computed at full vertex resolution
(`step4_5`): it's not just for speed, it matches how Pang/Vohryzek actually
report their numbers, so my results are comparable to theirs.

---

## 5. Connectome eigenmodes (graph Laplacian)

> "Lψ = -λψ, where L' is the normalized graph Laplacian... L' = D^(-1/2) L
> D^(-1/2), with L = ½[(D-A) + (D-A)ᵀ]."

**In plain terms:** Instead of using cortex *shape*, you can instead build a
network out of the brain's actual wiring (the connectome) and ask the same
"what are the resonant vibration modes of this network" question. A node's
"degree" *D* is how much total connection weight it has; the Laplacian *L*
roughly measures how different a node's activity is from its neighbors, and
the eigenmodes of the *normalized* version of that matrix are the network's
analogue of the geometric eigenmodes above.

**My thoughts:** This is the exact formula my `functions/calc_network_eigenmode_dense.m`
and `functions/calc_network_eigenmode_lowmem.m` implement — I use `Lnorm = I -
(dhalf .* network .* dhalf')` which is algebraically the same D^(-1/2) L
D^(-1/2) construction, just written to avoid ever building the dense D
matrix. One thing to double check: the paper explicitly symmetrizes with
½[(D-A)+(D-A)ᵀ] before normalizing — in my case the input networks (`C`,
`EDR_LRE`) are already symmetric by construction (built from a symmetric
distance matrix), so I skip that step. Worth a sentence in my thesis methods
explaining that this symmetrization is a no-op for my case, not something I
overlooked.

---

## 6. EDR (exponential distance rule) eigenmodes

> "...we fitted the variation of weights as a function of Euclidean distance,
> d... by an exponential function of the form e^(-αd)... resulting in an
> optimal empirical parameter value of α = 0.12."

**In plain terms:** On average, two points in the brain that are farther
apart tend to have weaker structural connections, and that falloff looks
close to a clean exponential curve. Instead of using the real (noisy)
connectome, you can build a synthetic one that follows this idealized
exponential rule perfectly, then ask whether *that* synthetic connectivity
is enough to explain brain function.

**My thoughts:** This is exactly the `lambda_pang = 0.12` constant hard-coded
in `step4_5_full_resolution_reconstruction.m` and `step5_static_reconstruction.m`
— now I know precisely where that number comes from (Pang's own fit,
reported in this Methods section) and *why* I use it instead of my own fresh
EDR fit for the "EDR binary"/"EDR continuous" comparison graphs: it's what
makes those two graphs comparable to Pang's own results, not just an
arbitrary choice. Also worth noting: their exponential fit used MATLAB's
`lsqcurvefit` (Optimization Toolbox) — I don't have that toolbox, so my own
EDR fits (Stage 1/2 of `connectome_harmonics.m`) use `fminsearch` minimizing
the same sum-of-squared-errors instead. Same objective, different solver —
worth a one-line caveat in my methods write-up.
