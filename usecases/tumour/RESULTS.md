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

   dt 0.005 against 0.0025: edges 8 um, dead 1.5% and cells 1.9% relative apart (pass). dt 0.01 against 0.005 is 19-27 um (fail). The cap acts on 0% of the cells in all dt runs, so this is **not** the voxel cap. **Where it comes from (`scripts/compare_dt_profiles.py`, `results/tum_dt_profiles.csv`):** the dependence is already in the oxygen field at 6 h. O2 at 100 um is 0.168, 0.202, 0.234, 0.257 mM for dt 0.02, 0.01, 0.005, 0.0025 (glucose 0.237-0.241 mM, no dt dependence), and the proliferative edge is 72, 112, 130, 142 um at 6 h; the profiles change little afterwards (the 0.0025 edge falls from 142 to 127 um by 72 h). The dt 0.00125 run (TUM-DT125, 24 h, 2026-11-05) completes the series: O2 at 100 um at 6 h is 0.168, 0.202, 0.234, 0.257, 0.271 mM for dt 0.02 ... 0.00125, successive differences 0.034, 0.033, 0.023, 0.014 mM (ratios 0.96, 0.70, 0.60), so the oxygen field is converging at about the rate of a first-order error but is still 0.014 mM (5%) from the 0.0025 value at the finest step. **The proliferative edge is converged between dt 0.0025 and 0.00125** (6 / 12 / 24 h: 141.8 / 141.3 / 142.6 against 144.5 / 141.2 / 144.3 um, at most 2.7 um). **At 24 h dt 0.005 is 13.6 um (against 0.0025) and 15.3 um (against 0.00125) below the finer runs, at the one-voxel criterion, so the 24 h sensitivity block (dt 0.005) has a common offset of about -14 um in the absolute proliferative edge; comparisons within the block are not affected.** At 72 h dt 0.005 is 5 um from dt 0.0025 (the dip at 72 h of item 6). So the cause is the O2 uptake-diffusion coupling (splitting), not the cell response and not the cap. The cap number (R = q X dt / c) halves with each dt halving. dt 0.005 is therefore the second reference (TUM5-*); the shipped dt 0.01 is the first (TUM-*).

2. **Voxel size (dt 0.005, TUM5-DX).** dx = 10, 15 (ref), 20, 30 um:

   | dx | prolif. edge (um, 72 h) | change | living edge change (um) | dead % change (rel) | cells change (rel) | cells capped for O2 (max over snapshots) | verdict |
   |---|---|---|---|---|---|---|---|
   | 10 | 109.5 | -12.4 | -5.9 | +0.4% | -2.2% | 2.2% | pass |
   | 20 | 133.4 | +11.4 | +1.0 | -1.4% | +0.3% | 1.0% | pass |
   | 30 | 141.5 | +19.5 | +12.1 | -6.3% | +2.4% | 6.6% | fail (proliferative edge > 15 um) |

   The edge keeps moving with dx by about one voxel per 5 um step (110, 122, 133, 141 um): within 15 um for dx 10-20, **not converged as dx -> 0**. At dt 0.01 the dependence on dx was non-monotone (77, 100, 116, 114 um for dx 10, 15, 20, 30) and the cap acted on 1-14% of the cells; at dt 0.005 it is monotone and smaller. Voxel volumes are 1500, 3375, 6000 and 13500 um^3 against a reference cell volume of 2494 um^3.

3. **Seeds, repeat, pFBA (dt 0.005).** Seeds 51 and 52 and the repeat of seed 50: proliferative edge 122, 126, 127, 125 um (spread about 5 um), living edge 188-192 um, dead 46.5-47.8%, cells 366-370 (all pass; the repeat is not bit-identical, as expected with 56 threads). pFBA as secondary objective: proliferative edge +4.6 um, living edge +0.2 um, dead -0.2%, cells +0.8% (pass). **By-product fluxes do not agree as well** (`scripts/compare_pfba_zones.py`, `results/tum5_pfba_*.csv`; zone means of the living cells at 72 h, x < 100 and 100-190 um, against the largest difference of REF with REP and the two seeds): O2 +9.3% / -0.2% and glucose +5.8% / +1.5% (proliferative / hypoxic), CO2 -21% / -14%, lactate +41% (proliferative), ammonium +24% / -26%, glutamine -45% / +17%; 16 of 22 substrates differ by more than 20% in at least one zone and 33 of 43 zone means exceed the seed/repeat spread (0.1-41%, mostly < 5%). Signs are kept. The bin profiles keep their shape for O2, glucose, CO2, lactate (correlation 0.97, 0.98, 0.96, 0.85) and not for glutamate, glycine, tryptophan (-0.36, -0.17, 0.32; fluxes of about 0.01 mmol/gDW/h and few cells per bin). So the zonation is robust to the objective, the individual by-product fluxes of Fig. 5B are one optimum (FVA replay still to do). The same checks at dt 0.01 give 98-101 um (seeds, repeat, pFBA), also within one voxel.

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
| Mechanism of the dt dependence of the proliferative edge | located in the O2 field already at 6 h (finding 1); with dt 0.00125 the O2 field converges (differences 0.034, 0.033, 0.023, 0.014 mM at 100 um) and the edge is converged between 0.0025 and 0.00125 (<= 2.7 um); the 24 h block at dt 0.005 is offset by about -14 um |
| dx not converged as dx -> 0 (about one voxel per 5 um) | stated as a limit; dx 10 and 20 within the criterion, dx 30 fails |
| Proliferative edge dips at 72 h in the reference, repeat and dt 0.005 runs (122, 125, 119 um against 126-131 um from 6 to 48 h) | not explained |
| FBA vs pFBA by-product fluxes | compared (finding 3): magnitudes differ beyond the noise, FVA replay to quantify the alternative optima |
| Vessel O2 of 0.5 mM versus tissue values | literature value not checked; the grid shows the column dies at <= 0.1 mM |
| Tumour figure (Fig. 5) regenerated from TUM5-REF | done (`scripts/make_figure.py`, `results/fig_tumour_ref.{png,pdf}`); per-cell rates use dry mass = 30% of the volume (fluid fraction 0.70); the submitted figure and Supplementary Fig. S9 used 75% and are 2.5 times higher; S8 and S9 regenerated from TUM5-REF (`scripts/make_supp_figs.py`, `results/supp/`); the submitted O2 row of S9 also lacked the cell volume |
| Elemental balance (UT-ELEM) of the MCF7 model | not done; FVA replay done 2026-11-06 (section at the end) |
| Table S5.4 lists glucose at x_max as 0.2 mM; the run has 0 (sink, all substrates) | corrected in `Tumour_edits.tex`; the other rows of the table agree with the run |

## FVA replay (2026-11-06; `common/fva_replay.py`)

Method: for a sample of living cells of one snapshot the linear problem of the cell is rebuilt (SBML model of its cell type, lower bound of every recorded exchange set to the `_lb` column, objective bound `max_growth_rate` not binding), solved with Gurobi (gurobipy 13.0.3 through cobra 0.32.1), and flux variability analysis is run on the recorded exchanges at fraction_of_optimum 0.99999 and 0.9999 (GLPK gives the same range limits to 2e-5, checked on all 29 tables before they were redone with Gurobi). Checks: the replay objective equals the recorded `growth_rate` (median 3.6e-11, max 3.1e-10 at 72 h, relative); infeasible replays are the cells with zero recorded growth; the recorded flux lies inside its FVA range in 100% of the cells. A range is called undetermined when it exceeds 1e-2 mmol/gDW/h; `width_ratio` is the median of range(0.9999) / range(0.99999): 10 for a unique optimum (the range is the tolerance), about 1 for a true degeneracy. The replay is independent of the Clp solver of the run. Tables in `results/fva/` (`*_cells.csv` one row per sampled cell, fraction and exchange; `*_summary.csv`).

**TUM5-REF_s50, 24 h (189 feasible of 196 living cells) and 72 h (181 of 195).** Share of cells with an undetermined exchange by distance from the vessel (`common/fva_zones.py`, `results/fva/zones_snap24.csv`, `zones_snaplast.csv`), 72 h: CO2, NH4+ and glutamine 100% in every zone (90-100% for lactate), width_ratio about 1 (a degeneracy that does not shrink with the tolerance), median range 0.17 (CO2), 0.09 (NH4+), 0.05 (lactate), 0.03 (glutamine) mmol/gDW/h, widest 0.58 (CO2); O2 and glucose are undetermined within 60 um of the vessel (100%; growth is capped by `max_growth_rate`, the uptake is not needed) and 48% at 60-120 um (69% / 55% at 24 h), and determined at 120-200 um (0%) where uptake is limiting. All 17 amino-acid, choline, Pi and ornithine exchanges are determined at the 1e-2 threshold (ranges up to 1e-2 for asn, glu, gly, orn, ser, tyr, which are not binding). Infeasible replays (7 at 24 h, 14 at 72 h) are cells with zero recorded growth. Replay objective: median 4.8e-11, max 4.5e-3 relative at 24 h (a few cells; not examined). This matches the earlier tumour replay (by-products degenerate, O2/glucose determined where limiting); the earlier statement that O2 and glucose are determined from 60 um is true only at 120 um and beyond in these runs. The pFBA contrast for the by-products (finding 3) is the same degeneracy.

## Unmodelled uptake (2026-11-06; `common/free_uptake_budget.py`, `results/free_uptake_budget_*.csv`)

Exchange and demand reactions with a negative lower bound that are not BioFVM densities (not in the YAML / `<transport_model>`) take up substrates from outside the simulation. Their uptake is kept open (author's decision, 2026-11-06; reported as non-limiting medium components, draft in `Review_round_01/Free_uptake_edits.tex`) except for those that no cell uses. Values are from the replay of the reference run, in `results/free_uptake_budget_<run>.csv` (per exchange: share of cells using it, uptake per gDW of biomass formed, amount drawn by the population, mean concentration in the domain). *E. coli* core: NH4 5.5 and Pi 3.7 mmol per gDW formed (100% of the cells); H+ and H2O used in the acetate regime (48% and 36% of the cell-times at the fermentation / reuse phases), solvent and buffer. Closing NH4 and Pi gives zero growth; H+ and H2O cannot be closed in the acetate switch (growth 0.89 of the open value at 100 min, infeasible at 300 min). Growth with every unused exchange closed equals the open growth (min and max ratio 1.000000000).
Tumour model (MCF7 core, TUM5-REF_s50, 24 h and 72 h): DM_Tyr_ggn_c_ (lower bound -1000, a demand reaction for tyrosyl-glycogenin, not a medium component) is used by 100% of the cells (0.034 mmol per gDW formed) and closing it gives zero growth; DM_gudac_c_ (guanidinoacetate, lower bound -0.0096) is used by 29% of the cell-times and closing it leaves growth unchanged; H+ and H2O (lower bound -0.02) can both be closed without change in growth. Decision pending on closing gudac, H+ and H2O in the tumour model (changes the by-product fluxes, not growth; would need a tumour rerun to quantify).
