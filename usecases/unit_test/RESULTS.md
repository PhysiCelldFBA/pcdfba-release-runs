# Unit test: results synthesis

Status 2026-10-29. All 70 runs of `runs_def.csv` were run on the pinned code (`review-runs-v1`, commit 7f2b1d0a; binary f397aa14e87ace45 built from it) and all finished (rc 0); `provenance.csv` has one row per run. They were run **locally** (8-core machine, single-thread runs 6 at a time), not on MN5; the 16-thread run of UT-THREADS was oversubscribed on 8 cores, so its timing means nothing (its result does). The full printout of the checks is `results/checks_local_pin.txt`, the CONV-BOX comparison is `results/ut_convbox.csv`, BASE and THREADS are in `results/ut_base_threads.csv`.
Not done: UT-DOMAIN (needs the addon hotfix), UT-TRIM, UT-ELEM for the four other models.

## Results

- **dx/dt invariance with the cap inactive (check A, 45 comparisons against the finest run printed in `results/checks_local_pin.txt` (dx 0.5, dt 0.001), limit 0.5%; the cap-number analysis and UT-CONV-BOX use the plateau dx 4, dt 0.001, where the cap cannot act, and give the same pass/fail pattern: dx >= 1 within 0.072% in volume, 0.046% in glucose taken):** 40 pass. The five failures are all dx 0.5 um at dt 0.01 and 0.02 (volume 12% and 49% from the plateau, glucose taken 16% and 63%; at dt 0.02 two of four cells remain). dx 1, 2 and 4 pass at every dt from 0.001 to 0.02 (volume within 0.054%, glucose taken within 0.041%). dx 0.5 is kept as the failure case.
- **Cap active (glucose D = 0, 1 mM; check B):** 31 of 32 pass; the voxel of the cell is emptied to within 0.035%. One failure: dx 1, dt 0.001, lowest glucose -3.5e-4 mM against a limit of -2.7e-5 mM (2 mu_max dt c0).
- **Ledger residual (check C):** 48 inactive and 47 capped checks pass; one failure: dx 0.5, dt 0.001, oxygen residual 2.3e-4 of gross against a limit of 2.7e-5 (the next largest is 4.6e-5, within its limit).
- **Carbon balance (check D):** 16 of 16 pass (largest difference 0.05%, limit 0.5%). **Elemental balance (E):** e_coli_core, 0 of 74 internal reactions unbalanced.
- **UT-CONV-BOX (8 cells growing to 32 in a 16 um box, 50 mM glucose, 26 mM left at 120 min):** compared with the plateau at 30, 60, 90 and 120 min, volume agrees within 0.03% for every dx >= 1 and every dt, and for dx 0.5 at dt 0.001 and 0.005 (0.046%); cell counts are identical. The six comparisons outside the 5% criterion (of 48 at 30, 60 and 90 min) are all dx 0.5 with dt 0.01 and 0.02: volume -9% and -20% at 30 min, -23% and -37% at 60 min, -24% and -51% at 90 min, with cells lost at 60 and 90 min.
- **Voxel sharing (1, 4, 8 cells, mixed):** all pass; glucose stays above -8e-6 mM; the 4-cell and mixed cases draw identical glucose (relative difference 0), so cells that cannot take a substrate do not dilute the others.
- **UT-ORDER (12 cells in one voxel, 5 creation orders, 1 and 4 threads):** all agree with the identity order at 1 thread to at most 1.1e-15 (criterion 1e-12); the lowest concentration is -3.4e-6 mM (solver tolerance).
- **UT-THREADS (1, 4, 16 threads):** identical in every printed digit (volume 45.0398, 32 cells, mean growth rate 0.732744).
- **UT-BASE (5 identical single-thread runs, default test):** identical results. The default test uses up its glucose (about 3e-61 at the end), so its end-state growth rate is 0.

## Cap number (UT-CAPRULE)

R = n q m_dw dt / (c V_vox) = q X dt / c is the fraction of a voxel's substrate that the cells in it ask for in one step (n cells in the voxel at the start; q at the run's own glucose; `scripts/cap_rule.py`, `results/ut_vs_plateau.csv`, `results/fig_cap_rule.png`). For these cells R = 4.5 n dt / dx^3 at 10 mM glucose.

| Experiment | R passing | R failing |
|---|---|---|
| UT-DXDT-A (1 cell, 10 mM; limit 0.5% at 120 min) | up to 0.18 (deviation <= 0.12%) | 0.36 (-12.3%), 0.72 (-48.7%) |
| UT-CONV-BOX (8 cells, 50 mM; limit 5% at 30-90 min) | up to 0.29 (<= 0.044%) | 0.58 (24.4%), 1.15 (51.3%) |

Both experiments change between 0.29 (largest R that passes) and 0.36 (smallest that fails), so the initial threshold is **R* of about 0.3**, with a margin of roughly 0.2-0.5 because the grid has few points around it. Below R* the deviation from the plateau is under 0.15%, above it at least 12%. R grows with the biomass in the voxel, so an initial R of 0.3 reaches R = 1, where the cap starts to act on a cell, once the biomass has grown a few times. In the colony the cap acts on a cell at its current R (`../colony/RESULTS.md`).

## Safe zone

Two conditions, both necessary:
1. **Voxel at least as large as the largest cell.** A cell is a point source in the voxel that holds its centre, which is only consistent if it fits in a voxel. Voxel volume (dx dy dz) >= the largest volume a living cell reaches before dividing (about 2.3-2.9 um3 for the cells here; take it from the output, not from the configured starting volume). With several cell types the largest sets the voxel.
2. **Initial R below about 0.3** for the substrates the cells take up (above that, the uptake cap and not the kinetics sets the dynamics; see the colony).

Outside the first condition: dx 0.5 um (voxel 0.06 times the largest cell) fails check A at dt >= 0.01 and loses cells in CONV-BOX; check C exceeds its limit at dt 0.001. dx 1 um (0.5 times the cell) passes A and CONV-BOX but check B overdraws one voxel by 0.035% of its content at dt 0.001. These two exceptions are taken as outside the validity of the framework, not as a reason to change the limits, which are kept as written.

| Use case | Voxel (um) | Voxel volume (um3) | Configured cell volume (um3) | Largest living cell measured (um3) | Voxel / largest cell |
|---|---|---|---|---|---|
| Unit test, dx 2 (default) | 2 x 2 x 2 | 8 | 1.3 | 2.29 (UT-BASE) | 3.5 |
| Unit test, dx 1 | 1 x 1 x 1 | 1 | 1.3 | 1.95 (CONV-BOX) | 0.5 |
| Unit test, dx 0.5 | 0.5^3 | 0.125 | 1.3 | 1.95 | 0.06 |
| Colony (shipped and pre-grown) | 2 x 2 x 2 | 8 | 1.3 | 1.49 (pre-grown, 10 min) | 5.4 |
| Acetate switch, shipped | 2 x 2 x 2 | 8 | 1.3 | 2.31 (earlier batch) | 3.5 |
| Acetate switch, ACE-REF design | 4 x 4 x 4 | 64 | 1.3 | 2.31 (earlier batch) | 28 |
| Cross-feeding | 1.25^3 | 1.95 | 1.3 / 1.8 | 2.2 / 2.9 (earlier run) | **0.67 (outside)** |
| Tumour | 15^3 | 3375 | 2494 | not measured | 1.35 (configured) |
| Motility | 5 x 5 x 10 | 250 | 8 | 8.0 (earlier run) | 31 |

Sources: voxel and configured volumes from the shipped configurations. Measured volumes come from the pinned unit-test and colony runs where they exist; the acetate, cross-feeding and motility values come from earlier runs (read only for cell size) and must be re-read from the definitive runs. The colony value is from the 10 min pre-grown design; a full-length run is likely to show larger cells (2.3 um3 in the unit test), still inside the zone. The tumour has no measurement yet. Condition 2 is not yet computed for the use cases other than the colony.

## Reading

The time step and the voxel size do not change the single-cell and shared-voxel results as long as the voxel is large enough for the cells and the initial cap number is below 0.3. The colony's own dependence on dt comes from the cap at its density (`../colony/RESULTS.md`).

## Open

1. UT-DOMAIN (needs the addon hotfix), UT-TRIM, UT-ELEM for the four other models.
2. Condition 2 (R) for the other use cases; the cross-feeding setting at dx 1.25 is outside condition 1 and its main-figure setting needs a decision.
3. Voxel/cell volumes of the acetate, cross-feeding, motility and tumour use cases to be re-read from definitive runs.
