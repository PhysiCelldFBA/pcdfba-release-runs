# Tumour (vessel-to-core column): results synthesis

Status 2026-11-05. Pinned code `review-runs-v1` (commit 7f2b1d0a, MN5 build, binary hash 4ff530093806a919), model MCF7 Zielinski 2017 core, 56 threads. Tables are in `results/` (`tum_all_*` dt 0.01 reference and series, `tum5_*` second reference, `sens_*` sensitivity block); scripts in `scripts/` (`design_def.py`, `design_sens.py`, `analyse_tum.py`). Run folders are not in git.

## Design

Box 360 x 120 x 30 um, two layers of voxels in z, dx = dy = dz = 15 um (shipped). Vessel at x = 0 (Dirichlet): O2 0.5 mM, glucose 0.5 mM, the amino acids and other nutrients with a fixed vessel value 1.0 mM; Dirichlet 0 on the x = 360 face for O2, glucose and lactate (as shipped; the y and z faces are closed). 240 cells placed from `cells_column.csv`, 72 h. Pre-equilibration of the fields by `diffusion_steps` BioFVM steps before the cells start (300 at dt 0.01, i.e. 3 min of diffusion; scaled as 300 x 0.01 / dt). Necrosis: R_ATPM lower bound threshold 1.07; `death_rate_increase` is added at each dFBA update.

Observables (fixed 2026-11-04, except the two marked): living edge = 95th percentile of x of the living cells (**chosen after seeing that the largest x is set by one cell**), proliferative edge = 95th percentile of x of living cells with growth rate above half of the maximum, 0.5 x 0.0137 1/h (the threshold is the pre-set one; a first analysis used 1e-3 1/h and was replaced), dead fraction = dead / all cells (dead cells are never removed), cell number. Criterion for agreement with the reference: both edges within 15 um (one voxel), dead fraction and cell number within 10% relative.

## Findings

1. **dt: the necrosis ramp must be scaled, and the proliferative edge then converges for dt <= 0.005 min.** `death_rate_increase` is added once per dFBA update, so the ramp is increase/dt per min^2. With the increase scaled as 1.67e-7 x dt/0.01 (and `diffusion_steps` scaled) the dead fraction is the same at all dt, and the remaining dependence is in the proliferative edge:

   | dt (min) | prolif. edge 72 h (um) | living edge (um) | dead % | cells | max cap number |
   |---|---|---|---|---|---|
   | 0.02 | 70 | 177 | 48.3 | 319 | 0.35 |
   | 0.01 (shipped) | 100 | 183 | 48.3 | 354 | 0.17 |
   | 0.005 | 119-127 (122 in TUM5-REF) | 191 | 47.0-47.3 | 368 | 0.08 |
   | 0.0025 | 127 | 185 | 48.0 | 375 | 0.04 |

   dt 0.005 against 0.0025: edges 8 um, dead 1.5% and cells 1.9% relative apart (pass). dt 0.01 against 0.005 is 19-27 um (fail). The cap acts on 0% of the cells in all dt runs, so this is **not** the voxel cap. The mechanism is not tested. The cap number (R = q X dt / c) halves with each halving of dt and the edge stops moving once it is below about 0.1, which is consistent with a splitting error that is controlled by R, but that was not isolated. dt 0.005 is therefore the second reference (TUM5-*); the shipped dt 0.01 is the first (TUM-*).

2. **Voxel size (dt 0.005, TUM5-DX).** dx = 10, 15 (ref), 20, 30 um:

   | dx | prolif. edge (um, 72 h) | change | living edge change (um) | dead % change (rel) | cells change (rel) | cells capped for O2 (max over snapshots) | verdict |
   |---|---|---|---|---|---|---|---|
   | 10 | 109.5 | -12.4 | -5.9 | +0.4% | -2.2% | 2.2% | pass |
   | 20 | 133.4 | +11.4 | +1.0 | -1.4% | +0.3% | 1.0% | pass |
   | 30 | 141.5 | +19.5 | +12.1 | -6.3% | +2.4% | 6.6% | fail (proliferative edge > 15 um) |

   The edge keeps moving with dx by about one voxel per 5 um step (110, 122, 133, 141 um): within 15 um for dx 10-20, **not converged as dx -> 0**. At dt 0.01 the dependence on dx was non-monotone (77, 100, 116, 114 um for dx 10, 15, 20, 30) and the cap acted on 1-14% of the cells; at dt 0.005 it is monotone and smaller. Voxel volumes are 1500, 3375, 6000 and 13500 um^3 against a reference cell volume of 2494 um^3.

3. **Seeds, repeat, pFBA (dt 0.005).** Seeds 51 and 52 and the repeat of seed 50: proliferative edge 122, 126, 127, 125 um (spread about 5 um), living edge 188-192 um, dead 46.5-47.8%, cells 366-370 (all pass; the repeat is not bit-identical, as expected with 56 threads). pFBA as secondary objective: proliferative edge +4.6 um, living edge +0.2 um, dead -0.2%, cells +0.8% (pass). Exchange fluxes of by-products (lactate, amino acids) were **not** compared between FBA and pFBA in this analysis. The same checks at dt 0.01 give 98-101 um (seeds, repeat, pFBA), also within one voxel.

4. **Conservation.** Ledger residual of the cells about 1e-14 of the gross flux for all substrates in the reference runs. Cap number at the dt 0.005 reference at most 0.08 and no cell capped (the cap acts only in the dx runs above).

5. **Sensitivity block (dt 0.005, seed 50, 24 h; 64 runs, 15 of which end with no living cell).** The proliferative edge at 24 h is within 2.1 um of its 72 h value in all four 72 h corner runs and the living edge within 14 um, so 24 h gives the zones; the dead fraction keeps rising (dead cells stay), so it is not a 72 h value. In the first 12-24 h window the edge of individual runs wanders by 5-15 um, so differences below that in the tables are not interpreted.

   Vessel O2 x vessel glucose (25 points; the centre is TUM5-REF at 24 h: living edge 192, proliferative edge 129 um, dead 34%), proliferative edge (um), 24 h; "-" = no proliferating cell:

   | O2 \ glucose (mM) | 0.1 | 0.25 | 0.5 | 1.0 | 5.0 |
   |---|---|---|---|---|---|
   | 0.05 | - (all dead) | - (all dead) | - (all dead) | - (all dead) | 13 |
   | 0.10 | - (all dead) | - (all dead) | - (all dead) | - (13 alive) | 43 |
   | 0.25 | 25 | 28 | 41 | 56 | 111 |
   | 0.5 | 114 | 120 | 129 (ref) | 131 | 182 |
   | 1.0 | 114 | 186 | 202 | 214 | 213 |

   Dead fraction at 24 h: 100% (0.05-0.1 mM O2, glucose <= 0.5), 55-65% (0.25), 20-35% (0.5), 3-14% (1.0). Oxygen sets the zones; glucose matters when oxygen is not limiting (0.5 mM O2: 129 -> 182 um from 0.5 to 5 mM glucose; 1.0 mM O2: saturates at about 214 um). **The column dies out or has no proliferating zone for a vessel O2 <= 0.1 mM unless glucose is 5 mM.** The reference value, 0.5 mM, is above the oxygen of tissue (0.05-0.13 mM was noted in `experiments.csv`; **the literature value has not been checked**), so the conclusion "zonation of a column of this depth" holds for a vessel that supplies at least 0.25 mM.

   Sobol (16 points, six transport factors x0.5-2, log-uniform, at vessel O2 0.5 and 0.1 mM): at 0.5 mM O2 15 of 16 points have a proliferative zone (26-194 um, median 126 um), none dies out; at 0.1 mM O2 7 of 16 die out and 11 of 16 have no proliferative zone. At 0.5 mM O2 the Pearson correlation of log2 factor with the living edge: O2 Vmax +0.68, O2 Km -0.64, glucose diffusion +0.63, glucose Vmax +0.14, glucose Km +0.06, O2 diffusion -0.02. 16 points are a screen, not variance indices; the correlations are indicative.

   Necrosis threshold of R_ATPM 0.80 / 0.96 / 1.07 (ref) / 1.18 / 1.34: proliferative edge 158 / 142 / 129 / 113 / 102 um, living edge 234 / 207 / 192 / 181 / 163 um, dead 19 / 29 / 34 / 38 / 46%. Monotone; the necrotic fraction is sensitive to this threshold (+-25%: -15 and +12 percentage points).

   Cap number at most 0.24 and no capped cell in any of the 64 runs. The ledger column of `sens_runs.csv` is an artefact (relative residual of ~1e293) in the 15 runs with no living cell because the gross flux is 0; `sens_summary.csv` leaves it blank there; one run (TUM-SOB_o0.5_p01) has 1e-3, the rest 0.

## Open items

| Item | Status |
|---|---|
| Mechanism of the dt dependence of the proliferative edge (not the cap) | not tested |
| dx not converged as dx -> 0 (about one voxel per 5 um) | stated as a limit; dx 10 and 20 within the criterion, dx 30 fails |
| Proliferative edge dips at 72 h in the reference, repeat and dt 0.005 runs (122, 125, 119 um against 126-131 um from 6 to 48 h) | not explained |
| FBA vs pFBA by-product fluxes | not compared |
| Vessel O2 of 0.5 mM versus tissue values | literature value not checked; the grid shows the column dies at <= 0.1 mM |
| Tumour figure (Fig. 5) regenerated from TUM5-REF | not done |
| Elemental balance (UT-ELEM) of the MCF7 model, FVA replay | not done |
| Table S5.4 lists glucose at x_max as 0.2 mM; the run configuration (config/shipped) has 0 | to check with the author |
