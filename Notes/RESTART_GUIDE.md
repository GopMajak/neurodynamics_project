# Restart Guide — Neurodynamics Project (post system-restart recovery)

## 1. What's confirmed safe locally (Desktop\Neurodynamics_Project)

- `Code\Ring_km.m`, `Code\km_dynamics_hcp.m`, `Code\km_cv_approach.m`
- `Code\simulation\`, `Code\analysis\`, `Code\graphs\`, `Code\helper_functions\` (incl. `expokit\`), `Code\colorcet\`
- `Code\DSI_release2_2011.mat`
- `Figures\` (Figures 1–7, 3 PNGs each), `Papers\` (~24 PDFs), `Experiments\Tracking.xlsx`
- **Empty:** `Notes\`, `Results\` — anything that lived only there before the restart is genuinely gone.

## 2. What's missing — and NOT safe to fabricate

A separate, more advanced working copy of this project existed in another Claude Code session
(Remote Control, session ref `528ba4`, "Replicate Vohryzek et al 2024 connectome harmonics
results") running on a different machine. Per that session's chat transcript, it contained:

- **`connectome_harmonics.m`** — the actual thesis pipeline replicating Vohryzek et al. 2024
  (EDR+LR method built on the Pang et al. 2023 toolbox). Three stages:
  - Stage 1: load connectome + surface
  - Stage 2: compute eigenmodes (`num_modes = 200` full-scale)
  - Stage 3: simulate dynamics (single run `5.0s`, speed scan `20 speeds × 3.0s`)
  - A `TEST_MODE` flag (line ~15) was added: subsamples 29,696 → 3,000 vertices, drops to
    10 modes, shortens Stage 3 runs (`1.0s` / `3 speeds × 0.5s`), and prefixes all outputs with
    `TEST_` so smoke-test results can't collide with real ones.
  - Loads `S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat` → `avgSC_L`
    (the HCP 255-subject group-average left-hemisphere connectome from Pang et al. 2023's
    toolbox — the correct dataset for the Vohryzek replication).
- **`budzinski2023_dKM/km_dynamics_hcp.m`** — a separate standalone script reproducing
  Budzinski et al. 2023's own figures. Correctly uses `DSI_release2_2011.mat` (Hagmann's DSI
  connectome) — this is *not* a bug, it's the dataset that paper actually used. Note: your local
  copy has `km_dynamics_hcp.m` at the top level of `Code\`, not inside a `budzinski2023_dKM\`
  subfolder — the other session had reorganized it.

**I only have a narrative description of these files from a pasted chat transcript — not their
actual source code.** I'm not going to invent a from-scratch rewrite of a scientific pipeline
and present it as if it were the recovered original; that risks silently introducing wrong math
(binning, normalization, exact parameter values) into your thesis results.

## 3. Recovery paths, in order of preference

1. **Wait for the other machine to reconnect.** I already sent a message to that session
   (queued, delivered once it's back online) asking it to share its current file contents.
   When it replies, I can pull the real files down and place them correctly.
2. **Check cloud sync.** If that other machine backs up to OneDrive, Google Drive, iCloud, or a
   git remote, the real files may already be sitting there. Worth checking before rebuilding
   anything by hand.
3. **Rebuild from source materials**, if neither of the above pans out:
   - Get the Pang et al. 2023 toolbox (public GitHub — search "connectome harmonics Pang 2023")
     and the `S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat` file it
     ships with.
   - Re-derive `connectome_harmonics.m` stage-by-stage against the actual Vohryzek et al. 2024
     methods section, using the toolbox's example scripts as a base — I can help write this once
     you have the toolbox and paper in hand, so the math is checked against a real source instead
     of a transcript summary.
   - Move `km_dynamics_hcp.m` into a `Code\budzinski2023_dKM\` subfolder if you want the same
     layout as the other session had.

## 4. Status

- Message sent to session `528ba4` — **queued, awaiting that machine reconnecting.**
- No files fabricated. Nothing overwritten locally.
