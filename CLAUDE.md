# Neurodynamics Project — Handoff / Context

Thesis: "Conduction Delays and the Dimensionality of Brain Dynamics: A Delayed-Kuramoto
Extension of the EDR+LR Harmonic-Mode Framework." Replicates Vohryzek et al. 2025 (PNAS,
EDR+LR harmonic-mode framework, builds on Pang et al. 2023 Nature) and extends it with
Budzinski et al. 2023's delayed-Kuramoto (dKM) model for conduction delays.

## Working setup (read before doing anything)

- **Two machines, two roles.** This machine (desktop or laptop) is for editing code and
  reviewing results. Heavy compute and all real HCP data live ONLY on the server `parana`
  (SSH: `gop@parana`), which has 64GB RAM / 24 cores. Local machines do NOT have the real
  data (~130GB of HCP `.nii`/`.mat` files, deliberately excluded from git — see `.gitignore`).
- **Always be explicit about which terminal a command belongs in**: "Terminal A — SSH
  session on parana" vs "Terminal B — local PowerShell." The user has repeatedly pasted
  commands into the wrong terminal; always state which one before giving a command.
- **New MATLAB files must sync to the server via `scp` before running there.** After
  editing/creating a `.m` file locally, it must be copied to
  `gop@parana:~/Neurodynamics_Project/...` (matching relative path) before it can run.
  Long-running server jobs should use:
  ```bash
  nohup matlab -nodisplay -nosplash -nodesktop -r "run('...'); exit" < /dev/null > logs_X.txt 2>&1 &
  disown
  ```
  The `< /dev/null` is required — without it, MATLAB spins forever on
  "Error reading character from command line" once detached from the terminal.
- **Never overwrite prior code versions.** Every fix/variant is a NEW file
  (e.g. `run_stage1.m` → `run_stage1_v2.m`), so original results stay next to novel results
  for side-by-side comparison. This is a standing, explicit user instruction.
- **Log every tested parameter/sweep/run to `Experiments/Tracking.xlsx`**, not just in chat
  — for supervisor review. This is a standing instruction (see project memory).
- **Git repo**: https://github.com/GopMajak/neurodynamics_project (private). Tracks only
  code, figures, notes, and `Tracking.xlsx` — never data, never PDFs, never `.mat` files
  (see `.gitignore`). Workflow: edit → `git add`/`commit`/`push` → `git pull` on the other
  machine → `scp` the relevant files to `parana` to actually run anything.
  Note: this session's auto-mode classifier blocks `git push`/`git remote add` as
  "out-of-place publication" — those commands need to be run by the user directly, not
  by Claude.

## Current state of the science (as of 2026-09-21)

**Baseline reproduction:** Vohryzek's own code reproduces the paper's published figures
successfully (`Code/reproduction/`). A corrected independent EDR fit
(`Code/stage1_anatomical_graphs_v2/run_stage1_v2.m`) matches the paper's values
(A=0.0658, lambda=0.1616, matching the recovered prior-session log almost exactly).

**Delay as a linear reconstruction basis — ROBUSTLY NULL.** Extensive testing (hybrid
basis test, residual variance test, generalization check across graph type and parameter
point, per-subject statistics on 255 subjects × 7 tasks) all converged on the same
conclusion: delay-informed eigenvector bases do NOT provide additional linear
reconstruction capacity for empirical task-activation maps beyond static geometry+long-range
harmonic modes. Static EDR+LR at N=20 modes remains the best reconstruction basis found.
See Tracking.xlsx sheets: "Hybrid Basis Test (Option A)", "Residual Variance Test (Option B)",
"Residual Variance Generalization", "Per-Subject Statistics".

**Pivot: delay as a generative dynamical feature — FIRST ROBUST POSITIVE FINDING.**
Rather than testing delay as a static linear basis, `Code/novel_delay_extension/
simulate_FC_delay_vs_nodelay.m` and its follow-up `simulate_FC_delay_robustness.m`
simulate the delayed-Kuramoto model directly on the EDR-only/EDR+LR structural graphs
(identical oscillator frequency across nodes, f_mu=10Hz, kappa=6) and compare the
resulting simulated FC to empirical group-average resting-state FC (255 HCP subjects,
Fisher-z averaged; extraction pipeline in `extract_parcellate_rfMRI_all_subjects.m`).

Result, confirmed robust across 5 random seeds x 2 graphs (10/10 combinations, sign test
p~0.001):
- Delay beats a no-delay Kuramoto control on the SAME graph, peaking sharply at
  **2.5-3.0 m/s** conduction speed (replicated almost exactly across seeds).
- EDR+LR peaks at 3.0 m/s: mean r=0.345 +/-0.044 vs no-delay 0.231 +/-0.036.
- EDR-only peaks at 2.5 m/s: mean r=0.325 +/-0.041 vs no-delay 0.137 +/-0.120.
- EDR+LR beats EDR-only throughout (consistent with geometry+long-range+delay narrative).
- The no-delay control is the SEED-UNSTABLE one (for EDR-only, SD is nearly as large as
  the mean, and it goes negative for one seed) -- delay doesn't just improve the mean, it
  stabilizes the model against a real failure mode.

See Tracking.xlsx sheets: "Simulated FC Delay vs NoDelay", "Simulated FC Robustness".

## Open threads / possible next steps

- Formalize the multi-seed robustness result with a proper statistical test (currently
  just sign test + descriptive mean/SD across 5 seeds).
- Check whether 2.5-3 m/s is physiologically plausible for unmyelinated cortico-cortical
  fibers (cross-reference Budzinski et al. 2023's cited conduction-velocity ranges).
- Reconcile this parcel-level (180-node, Euclidean centroid distance) peak of ~2.5-3 m/s
  with the earlier 998-node HCP DSI validation's peak of ~9 m/s (see "Speed Scan - HCP 998"
  sheet) -- likely a graph-resolution / distance-metric difference, not yet investigated.
- Consider whether to write this up as the thesis's key novel contribution: delay as a
  generative/dynamical complement to geometry+long-range structure (not a replacement,
  and not useful as a linear reconstruction basis) -- this reframes the earlier "delay
  needs to be taken into account" framing on solid empirical footing.
