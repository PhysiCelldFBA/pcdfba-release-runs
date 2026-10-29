# Colony: results synthesis (pre-grown study)

Status 2026-10-29. Everything below comes from the 13 runs in `runs_def.csv`, all on the pinned code (`review-runs-v1`, commit 7f2b1d0a, MN5-built binary 226118d62b231202; see `provenance.csv`). Numbers are in `results/colony_summary.csv` (per run at 2 min, 6 min and the last snapshot), `results/pg_profiles.csv` (radial profiles) and `results/cap_check_*.csv`.
The full-length growth runs (COL-REF, COL-SEED, COL-REP) are **not run yet**; nothing here is a statement about the 600 min colony.

## Design

A disc of 87,685 cells (R0 200 um, hexagonal packing) is placed at t = 0 with 4 mM glucose and 5 mM O2 on the box faces (dx 2, dz 2, 2D), necrosis rate x10, 10 min, ledger on. This isolates the zonation (necrotic core, acetate-using ring) from the slow growth phase and costs 1 h per run instead of 10 h.
Observables: dead fraction, necrotic-core radius (outer edge of the contiguous region that is at least 50% dead, 10 um annuli), acetate ring (annuli where at least 20% of living cells take up acetate), and the fraction of living cells whose recorded uptake bound is below the Michaelis-Menten bound ("capped", ratio < 0.99). Cap number R = q X dt / c is defined in `scripts/cap_check.py`.

## Findings

1. **Seed and repeat spread are negligible.** Seeds 50, 51, 52 and the seed-50 repeat give 41.8% dead at 10 min in all four, core radius 130 um, ring 120-130 to 190 um. (The only variation is one ring edge moving by one annulus.) The cells are placed identically in this design, so this is the spread of the dynamics only.

2. **The structure depends strongly on dt at the shipped setting and converges only for dt <= 0.0005 min.**

   | dt (min) | dead % at 2 min | dead % at 6 min | core / ring at the last snapshot (um) | capped glucose / acetate (last snapshot, %) |
   |---|---|---|---|---|
   | 0.01 | 0.1 | 2.0 | none / 10-200 (10 min) | 96 / 95 |
   | 0.004 | 1.0 | 12.2 | 100 / 80-190 (10 min) | 79 / 92 |
   | 0.002 | 2.9 | 24.6 | 120 / 110-190 (10 min) | 55 / 89 |
   | 0.001 (shipped) | 6.5 | 37.3 | 130 / 120-190 (10 min) | 29 / 60 |
   | 0.0005 | 13.1 | 45.6 | 140 / 140-190 (10 min) | 2 / 2 |
   | 0.00025 | 23.0 | 48.5 | 140 / 140-190 (6 min) | 0 / 0 |
   | 0.0001 | 38.8 | - | 130 / 140-180 (2 min) | 0 / 0 |

   Between dt 0.0005 and 0.00025 the late state agrees (core 140 um in both, ring 140-190 um, dead fraction 45.6% vs 48.5% at 6 min, 6% relative). dt 0.001 is not converged (core 130 um, 41.8% dead at 10 min against 47.4% at dt 0.0005). dt 0.01 gives no necrotic core.

3. **What changes with dt is the uptake cap, not the kinetics.** The fraction of capped cells falls from about 95% (dt 0.01) to 2% (dt 0.0005) and 0% (dt <= 0.00025), and the core radius tracks it. With the cap inactive the median uptake equals the Michaelis-Menten bound, and the late structure no longer changes with dt. The cap number at dt 0.0005 is below 1 for most cells (median R 0.3 glucose, 0.7 acetate); at dt 0.001 it is 0.7 and 1.6-1.7.

4. **The time at which the core forms keeps shortening as dt falls, even with the cap inactive.** At 2 min the dead fraction is 13.1, 23.0 and 38.8% for dt 0.0005, 0.00025 and 0.0001, with no sign of levelling. dt 0.00025 is at its plateau by about 5 min and dt 0.0001 reaches 38.8% by 2 min, which suggests the same end state is reached earlier, but the dt 0.0001 run stops at 2 min and cannot confirm it. Timing of necrosis (when) is therefore not a converged quantity at any dt tested; the end-state structure (where) is.

5. **Voxel size matters less than dt at dt 0.001.** dx 1 and dx 4 give 40.2% and 37.3% dead at 10 min (core 130 and 120 um) against 41.8% (130 um) at dx 2: changes of 4% and 11% in dead fraction, 0 and 8% in core radius. The cap is active in 35-73% of cells in all three, so this comparison mixes discretisation with cap effects (see the open items).

6. **Necrosis rate x1 instead of x10 changes the state at 10 min** (22.2% dead, core 100 um at 10 min, against 41.8% and 130 um), as expected for a death rate ten times lower. The runs do not show whether x1 reaches the same core given more time; the x10 speed-up is therefore an acceleration whose equivalence is not demonstrated.

7. **Cost.** Wall time on one 112-core node for 87,685 cells: dt 0.01 6 min, 0.004 17 min, 0.002 32 min, 0.001 57-67 min, 0.0005 80 min (all for 10 simulated min); dt 0.00025 60 min for 6 simulated min; dt 0.0001 45 min for 2 simulated min. Cost scales as 1/dt.

## What can be said to the reviewers (R2.M1.5, R2.M1.6, R2.M1.8, R2.M3)

- The colony zonation (necrotic core surrounded by an acetate-using ring) is present at every dt up to 0.004 and disappears only at dt 0.01. It is a qualitative result of the model.
- Core radius, ring inner edge and dead fraction are quantitative results that depend on dt at the shipped setting through the uptake cap, and converge once the cap is inactive (dt <= 0.0005 in this colony). The shipped setting (dt 0.001) lies 7% below that limit in core radius and about 12% in dead fraction at 10 min.
- Stochastic spread in this design is below the resolution of the observables (annulus width 10 um).
- The cap number R (unit-test experiment UT-CAPRULE) predicts where the cap acts; the colony data are the in-situ check.

## Open items and replan

| Item | Status | Proposal |
|---|---|---|
| COL-PG (seeds + repeat) | done | none |
| COL-DT | done for dt 0.0001-0.01; dt 0.00025 stops at 6 min and dt 0.0001 at 2 min (2 h debug limit) | optional: dt 0.0001 to 6 min in the standard queue (about 135 min), to confirm the plateau |
| COL-DX | done at dt 0.001 only, where the cap is active in 35-73% of cells | add dx {1, 4} at dt 0.0005, 6 min (dx 4 about 1 h; dx 1 about 85 min), so dx is read with the cap inactive |
| COL-NEC | done at dt 0.001, 10 min only | x1 at dt 0.001 to 40 min in the standard queue (about 4 h), to see whether x1 reaches the x10 core; otherwise state x10 as an acceleration |
| COL-REF / SEED / REP | not run | decision needed on dt (0.0005 recommended by the evidence above, about 20 h per run; 0.001 about 10 h). The full-growth colony starts small, so its cap fraction may differ from the pre-grown one: run the reference at both dt for one seed and check |
| COL-BND, COL-KIN, COL-PFBA | not run | run with the converged dt; to keep the cost down use R0 100 um after one check run (new experiment COL-R0: R0 100 vs 200 at the same settings, 6 min) |
| COL-FVA | not run | the pre-grown snapshots on the pin already carry the `_lb` columns; the replay can start on COL-PG_s50 and COL-DT_dt0.0005_s50 now |
| COL-SCALE | wall times above | threads scaling not run |
| Ledger (R2.M1.7) | recorded in all runs | residuals not yet analysed |
