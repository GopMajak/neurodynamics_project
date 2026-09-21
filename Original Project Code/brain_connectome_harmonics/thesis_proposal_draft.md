# Thesis Proposal — Draft 2

**Thesis title (for registrar office):** [DRAFT — discuss with supervisor] Conduction Delays and the Dimensionality of Brain Dynamics: A Delayed-Kuramoto Extension of the EDR+LR Harmonic-Mode Framework

**Course description (for registrar office):** [DRAFT — discuss with supervisor] This project tests whether distance-dependent conduction delay, introduced via a time-delayed Kuramoto model, reduces the number of structural eigenmodes required to reconstruct human brain activity, and whether this delay-related mechanism accounts for part of the reconstruction advantage previously attributed to rare long-range cortical connections.

---

## Introduction

In systems neuroscience, a fundamental question is understanding how the brain's structure and connectivity shape spatiotemporal dynamics. The underlying mechanisms by which structural features translate into the dynamic patterns of activity still remain unclear. In recent work this problem has been explored through the harmonic decomposition of anatomical connectivity, essentially viewing the brain's structural graph as a static substrate whose eigenmodes are used to reconstruct observed functional activity. Pang et al. (2023) showed that eigenmodes derived from cortical *geometry* alone reconstruct functional activity better than eigenmodes of the empirical structural connectome — a result that, on its face, suggests that the topologically complex long-range connections present in the true connectome contribute little beyond what geometry already explains.

Vohryzek et al. (2025) revisit this question with a more surgical construction: rather than using the full empirical connectome, they build a dense exponential-distance-rule (EDR) baseline and splice back in only the rare long-range (LR) connections that are statistical outliers (>3 SD above the local mean, >40mm). This EDR+LR basis reconstructs both functional connectivity and task-evoked activity with *fewer* modes than geometry, EDR-only, or the full connectome — specifically within the first 1–20 modes, where most of the reconstructible structure lives. The apparent tension with Pang is reconciled by construction, not contradiction: it is not long-range connectivity in general that helps, but the rare, spatially specific LR exceptions layered on top of an EDR baseline. Their results suggest a low-dimensional set of structurally-derived modes is enough to capture most of the brain's spatiotemporal organization — and they explicitly flag temporally-evolving descriptions of brain dynamics as an unaddressed direction for future work.

This is the gap this project addresses: time. The EDR+LR framework treats connectivity as instantaneous, taking no account of the finite conduction delays that arise from the physical distances signals travel along axonal pathways. Budzinski et al. (2023) show that in a time-delayed Kuramoto model, these delays are not a negligible detail but an active driver of spatiotemporal structure: heterogeneous, distance-dependent delay rotates the eigenspectrum of the coupling operator, shifting dominance away from the trivial synchronous mode and toward specific higher-order modes that produce traveling-wave patterns consistent with those observed empirically in human cortical recordings. Koller, Schirner, and Ritter (2024) extend this logic directly onto empirical structural connectomes: they show that conduction delay is *necessary* — not merely distance-dependent coupling — for a connectome's long-range/instrength structure to produce correctly-directed traveling waves, and that removing delay abolishes the effect entirely. Notably, they explicitly identify the relationship between geometry, connectome structure, and this delay-dependent dynamics as an open question, directly pointing back to Pang's result.

This proposed project asks: does introducing distance-dependent conduction delay into the harmonic-mode framework improve reconstruction of brain dynamics — specifically, does it reduce the number of modes required — and is the previously reported importance of rare long-range connections partly attributable to delay-related temporal structure that the static framework does not represent?

## Significance

Pang (2023) and Vohryzek (2025) currently offer two partially conflicting stories about what anatomical feature of the connectome matters for brain function: geometry (Pang) versus a specific sparse long-range structure layered on geometry (Vohryzek). Both are static analyses. If conduction delay — a physical, measurable property of axonal pathways — can independently account for some or all of the efficiency gain Vohryzek attributes to long-range connections, this reframes the question from *which anatomical feature* matters to *which anatomical feature matters once time is taken into account*, and offers a mechanistic bridge between the static harmonic-mode literature (Pang, Vohryzek) and the dynamical traveling-wave literature (Budzinski, Koller), which have so far developed largely independently. Practically, a demonstration that delay shapes the dimensionality of reconstructible brain activity would motivate conduction-velocity-sensitive measures (e.g., myelination, axon caliber) as a target in conditions where white-matter conduction is disrupted (e.g., demyelinating disease, aging), beyond the connectivity-strength measures such studies typically focus on.

## Method

To test whether conduction delay alters which anatomical features best account for brain spatiotemporal dynamics, this project extends the harmonic-mode framework of Vohryzek et al. (2025) by incorporating the time-delayed Kuramoto model (dKM) of Budzinski et al. (2023), and applying both to the same HCP dataset and reconstruction targets.

### Dataset

This project uses the same HCP data as Vohryzek et al. (2025): 255 participants with complete resting-state and task-based fMRI (47 contrasts across 7 domains), plus a 100-participant validation subset, parcellated to Glasser360 (left hemisphere, 180 parcels). Structural connectivity comes from the same high-resolution group-average dMRI connectome (29,696 cortical vertices, left hemisphere, medial wall excluded) provided by Pang et al. (2023).

### Structural graph construction

As in Vohryzek et al. (2025), four anatomical graph representations are built following their published parameters: geometry (Laplace-Beltrami eigenmodes of the cortical surface), EDR binary, EDR continuous, and EDR+LR (EDR continuous plus rare long-range exceptions: >3 SD above the expected weight at a given distance, and >40mm). This project's primary comparisons focus on EDR continuous (hereafter "EDR-only") and EDR+LR, the two conditions that isolate the effect of the rare long-range exceptions; geometry and EDR binary are retained as replication checks against Vohryzek's published static results.

### Delay-informed reconstruction pipeline

Vohryzek's "modes" are eigenvectors of a static, real, symmetric graph Laplacian, used to project empirical fMRI data and measure reconstruction accuracy as a function of the number of modes retained (N). Budzinski's dKM, by contrast, does not define modes this way: its dynamics are governed by a complex, delay-rotated coupling operator whose eigenvectors are not orthogonal in the usual sense and cannot be substituted directly into Vohryzek's projection formula. This project resolves that mismatch with a two-step pipeline, run separately for the EDR-only and EDR+LR graphs:

1. **Static baseline** (replication): compute the graph Laplacian eigenmodes of each structural graph, project the 47 HCP task maps and resting-state FC onto the first N modes, and compute reconstruction accuracy (correlation and mean squared error, following Vohryzek) as a function of N. This reproduces Vohryzek's published accuracy-vs-N curves and establishes the baseline this project's delayed condition is compared against.
2. **Delay-informed basis** (novel): simulate the dKM on the same structural graph across a range of conduction speeds, using distance-dependent delay τ_jk = d_jk / v. From the resulting steady-state phase time series, extract a low-dimensional, orthonormal spatial basis via PCA/SVD (complex-valued PCA where appropriate, following established methodology in the cortical traveling-wave literature). Project the same empirical fMRI targets onto this delay-informed basis and compute the identical accuracy-vs-N curve as in the static case.

Because the delay-informed basis is orthonormal by construction (a property of PCA/SVD), it supports the same projection and accuracy-vs-N machinery as Vohryzek's static eigenmodes, without requiring an analytic solution to the non-Hermitian eigenvector problem posed by Budzinski's operator directly.

### Implementation pipeline

| Stage | What | Tool | File |
|---|---|---|---|
| 0 | Get processed HCP data: structural connectome, cortical surface mesh, Glasser360 parcellation, 47 task contrasts + resting-state FC | HCP + Pang et al. (2023) / Vohryzek et al. (2025) public data release | already provided as processed `.mat`/`.vtk` files (no raw FreeSurfer reconstruction or DTI tractography needed — Pang and Vohryzek already publish the processed connectome and surfaces) |
| 1 | Load geometry, cortex mask, parcellation | MATLAB | `connectome_harmonics.m` (Setup) — existing |
| 2 | Validate EDR fit + long-range exception detection, Glasser360 parcel resolution | MATLAB | `connectome_harmonics.m` (Stage 1) — existing |
| 3 | Build full vertex-resolution EDR-only and EDR+LR connectome graphs | MATLAB | `connectome_harmonics.m` (Stage 2, connectome construction) — existing |
| 4 | Compute connectome Laplacian + static eigenmodes (harmonics) for both graphs | MATLAB | `functions/calc_network_eigenmode_lowmem.m` — existing |
| 5 | Reconstruct empirical task/FC data from static eigenmodes, accuracy vs. N (replicate Vohryzek Fig. 1) | MATLAB | `step5_static_reconstruction.m` — to be written |
| 6 | **Extension:** time-delayed Kuramoto (dKM) simulation on both graphs across conduction speeds | MATLAB | `simulation/simulate_dKM_fast.m`, `connectome_harmonics.m` (Stage 3) — existing |
| 7 | **Extension:** extract delay-informed basis (PCA/SVD of simulated steady-state dynamics) | MATLAB | `step7_delay_informed_basis.m` — to be written |
| 8 | **Extension:** reconstruct empirical task/FC data from delay-informed basis, accuracy vs. N | MATLAB | `step8_delay_reconstruction.m` — to be written |
| 9 | **Extension:** statistical comparison, static vs. delay-informed × EDR-only vs. EDR+LR (H1/H2) | MATLAB | `step9_compare_conditions.m` — to be written |
| — | **Validation:** reproduce Budzinski's own ring-graph and 998-region test cases before running on the EDR/EDR+LR graphs | MATLAB | `budzinski2023_dKM/` — existing (partial) |

### Analysis plan

The design is a 2×2 comparison — {static, delay-informed} × {EDR-only, EDR+LR} — evaluated via paired comparisons (following Vohryzek's own paired t-tests across the 47 task contrasts) at matched N, primarily N=20 (Vohryzek's reported range for peak EDR+LR advantage), with the full accuracy-vs-N curve reported for each condition. Four comparisons address the two hypotheses below:

- Delayed EDR+LR vs. static EDR+LR — does delay reduce N *within* the long-range-containing graph?
- Delayed EDR-only vs. static EDR-only — does delay alone, without long-range exceptions, already improve efficiency?
- **Delayed EDR-only vs. static EDR+LR** — does delay alone recover the efficiency Vohryzek attributed to long-range connections?
- Delayed EDR+LR vs. delayed EDR-only — does long-range structure still confer an advantage once delay is present?

### Hypotheses

**Hypothesis 1 (primary):** The delay-informed basis of the EDR+LR connectome will reconstruct empirical task-evoked activity and resting-state FC with fewer modes (equivalent accuracy at lower N) than the static EDR+LR eigenmodes.
*Null:* Static and delay-informed EDR+LR bases require statistically indistinguishable N for matched reconstruction accuracy.

**Hypothesis 2 (secondary):** The delay-informed basis of the EDR-only connectome will reconstruct empirical activity as efficiently (statistically indistinguishable N) as the static EDR+LR basis — indicating that the previously reported advantage of rare long-range connections is at least partly attributable to delay-related temporal structure that the static framework omits.
*Null:* The delayed EDR-only basis remains less efficient than the static EDR+LR basis, indicating long-range connections retain an effect independent of delay.

## Timeline

Project runs September 9, 2026 – end of April 2027. Fall/Winter reading-break weeks fall inside the ranges below as buffer but aren't dated precisely — adjust once the 2026–27 academic calendar is confirmed.

| Stage | Dates | Task | How? |
|---|---|---|---|
| 1 | Sep 9 – Sep 25, 2026 | Data collection: acquire and organize the HCP dataset (255-participant resting-state + task fMRI, 47 contrasts across 7 domains; 100-participant validation subset; Glasser360 parcellation; high-resolution group-average dMRI connectome) | HCP database (already processed — no raw FreeSurfer/tractography needed) |
| 2 | Sep 28 – Nov 6, 2026 | Reproducing prior results: replicate Vohryzek et al. (2025)'s static Laplacian-eigenmode reconstruction — build the connectome graphs (geometry, EDR binary, EDR continuous, EDR+LR: local + rare long-range edges), compute the connectome Laplacian + harmonics (eigendecomposition), and reproduce their published accuracy-vs-N curves | MATLAB |
| 3 | Nov 9 – Dec 18, 2026 | Introducing new code: develop the novel delay-informed reconstruction pipeline (dKM simulation across conduction speeds + PCA/SVD basis extraction from the simulated dynamics) | MATLAB |
| 4 | Jan 4 – Jan 15, 2027 | Validate new pipeline: reproduce Budzinski's own ring-graph and 998-region test cases; confirm the delay-informed basis converges to the static basis as delay → 0 | MATLAB |
| 5 | Jan 18 – Jan 29, 2027 | Pilot run: full static/delay-informed × EDR-only/EDR+LR comparison on a subset of task contrasts, to confirm the analysis and stats pipeline before scaling up | MATLAB |
| 6 | Feb 1 – Feb 26, 2027 | Run the complete reconstruction analysis: all 47 task contrasts + resting-state FC, both graphs, static and delay-informed; paired statistical comparisons at matched N; sensitivity analysis across conduction speeds | MATLAB |
| 7 | Ongoing (Sep 9, 2026 – Apr 16, 2027; concentrated Mar 1 – Apr 16) | Writing/Revisions: thesis draft assembly and supervisor/committee feedback rounds | Word / LaTeX |
| 8 | Apr 12 – Apr 19, 2027 | Presentation Preparation: defense slides | Slides (PowerPoint/Keynote) |
| 9 | Apr 19 – Apr 30, 2027 | Defense & final submission | — |

## References

Budzinski, R. C., Nguyen, T. T., Benigno, G. B., Đoàn, J., Mináč, J., Sejnowski, T. J., & Muller, L. E. (2023). Analytical prediction of specific spatiotemporal patterns in nonlinear oscillator networks with distance-dependent time delays. *Physical Review Research*, 5(1), Article 013159. https://doi.org/10.1103/PhysRevResearch.5.013159

Koller, D. P., Schirner, M., & Ritter, P. (2024). Human connectome topology directs cortical traveling waves and shapes frequency gradients. *Nature Communications*, 15, 3570. https://doi.org/10.1038/s41467-024-47860-x

Pang, J. C., Aquino, K. M., Oldehinkel, M., Robinson, P. A., Fulcher, B. D., Breakspear, M., & Fornito, A. (2023). Geometric constraints on human brain function. *Nature (London)*, 618(7965), 566–574. https://doi.org/10.1038/s41586-023-06098-1

Vohryzek, J., Sanz-Perl, Y., Kringelbach, M. L., & Deco, G. (2025). Human brain dynamics are shaped by rare long-range connections over and above cortical geometry. *Proceedings of the National Academy of Sciences - PNAS*, 122(1), Article e2415102122. https://doi.org/10.1073/pnas.2415102122
