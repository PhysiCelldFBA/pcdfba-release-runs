# Metabolism-driven motility: results on the pinned code

Status 2026-11-07. Reference and seed runs done (8 runs); convergence (dt, dx), pFBA and the sensitivity block done (34 runs, MN5 job 48710127, 2 min 5 s, all rc 0); FVA replay and the figure are not done.

## Design (`scripts/design_def.py`, `runs_def.csv`, `overrides_def.csv`, `experiments.csv`)

Two variants of the example, each from its shipped settings file (they differ in Km(glucose), `glucose_threshold_low` and `motility_cost_at_vmax` only): **MOT-REF**, motility-only (Km 0.5 mM, threshold 0.5 mM, cost 48.9 mmol/gDW/h) and **MOT-GR**, growth recovery (Km 0.02 mM, threshold 0.35 mM, cost 66.4). One cell at x = -140 um in a 300 x 100 um 2D domain (dx 5 um), glucose boundary 0.65 mM, 600 min, dt_diffusion = intracellular_dt = 0.001 min (harmonised configuration, PR #16), 6 threads, ledger on, snapshots every 5 min. The saved settings of the runs equal the rebuilt configs (`common/check_run.py`: OK). Seeds 50, 51, 52 and a repeat of seed 50 for each variant. Binary 847df0a2939f47f1 (MN5-built from the pinned commit 7f2b1d0a, job 48645327, 8 runs at once in 2 min 43 s; `binaries.tsv`, `provenance.csv`). Analysis: `scripts/analyse_mot.py` -> `results/mot_series.csv`, `results/mot_summary.csv`.

## Results

| run | x of the founder cell at 10 h (um) | cells at 10 h | first division (min) | mean speed (um/min) | share of time in ATP mode |
|---|---|---|---|---|---|
| MOT-REF s50 (and repeat) | 126.1 | 1 | none | 0.443 | 0.99 |
| MOT-REF s51 | 127.1 | 1 | none | 0.445 | 0.99 |
| MOT-REF s52 | 124.8 | 1 | none | 0.441 | 0.99 |
| MOT-GR s50 (and repeat) | 135.8 | 4 | 425 | 0.461 | 0.46 |
| MOT-GR s51 | 135.6 | 4 | 415 | 0.464 | 0.86 |
| MOT-GR s52 | 137.7 | 4 | 425 | 0.464 | 0.98 |

1. **Motility-only variant (Km 0.5 mM): the regime holds.** The cell stays in ATP mode for the whole run, never divides, moves slowly during the first 5 h (-137.6 um at 60 min, -105 at 300 min) and then fast (-42 at 420, +9 at 480, +126 um at 600 min). Position at 10 h 124.8-127.1 um (spread 1.8% over three seeds); the earlier harmonised runs gave +122 um (122-125 um over three seeds, 2.4%). The same-seed repeat is identical to the seed-50 run here (126.1 um), so on 6 threads a single-cell run is reproducible.
2. **Growth-recovery variant (Km 0.02 mM): the sequence holds, the timing is later.** The founder reaches the right-hand part of the domain (x about 135 um, near the wall at 150) at about 300 min in ATP mode, switches to biomass mode (seed 50: mode share 0 from 300 min; seed 51: part of the cells stay in ATP mode, share 0.5-0.75), divides at 415-425 min and has 4 cells at 10 h in all four runs. **The earlier harmonised runs gave the first division at 355-365 min and 5-12 cells at 10 h (5-16 over four seeds); on the pinned code it is 60 min later and 4 cells.** The cause was not tested (earlier runs used the hot-fixes build before the pin: the cap, the race conditions and the equal-step rule are in the pin); the difference is larger than the spread of the three seeds here (415-425 min) and is a shift in timing, not in the sequence (motility, arrival, switch, division). The mode share over time differs between seeds (0.46, 0.86, 0.98) because the cells switch between modes near the glucose threshold.
3. The manuscript text of Fig. 7 (position +122 um at 10 h, first division at 355-365 min, 5-12 cells) must be updated to the numbers above, and the figure regenerated from MOT-REF_s50 and MOT-GR_s50 (`Fig. 7` script not yet ported to this folder).

## Not done
Convergence in dt (0.002, 0.0005 min) and dx; pFBA and FVA replay (objective per cell mode: ATP in mode 1, biomass in mode 0; the recorded column `metabolic_mode`); sensitivity (Km x glucose boundary, ATP cost, glucose D); figure; the share of cells in a voxel and the cap number of the runs (single cell: the cap is expected to be inactive); `cap_check.py`-style check of the cap fraction for this use case.

## Block 2 (2026-11-07; MN5 job 48710127; `scripts/design_def.py`, experiments MOT-DT/DX/PFBA/SENS-REF/GR, seed 50, 600 min; `results/mot_summary.csv`)

Position of the founder cell at 10 h (um; wall at +150), cells at 10 h, first division (min):

| run | MOT-REF (motility-only) | MOT-GR (growth recovery) |
|---|---|---|
| reference (dt 0.001, dx 5) | +126.1, 1 cell | +135.8, 4 cells, 425 |
| dt 0.002 | +140.0, 1 | +147.9, 4, 430 |
| dt 0.0005 | **-140 (never moves, speed 0)**, 1 | +133.3, 4, 430 |
| dx 2.5 | +142.9, 1 | +147.9, 4, 440 |
| dx 10 | **-140 (never moves)**, 1 | +130.4, 4, 415 |
| pFBA | +126.1, 1 (identical) | +135.8, 4, 425 (identical) |

1. **pFBA changes nothing** (position, speed, cells and division identical to the digit in both variants): the regimes do not depend on the choice among alternative optima.
2. **MOT-GR is robust to dt and dx in regime** (4 cells, first division 415-440 min in every run; position 130-148 um, near the wall, so position at 10 h is a weak observable: the wall at 150 caps it). Against the criterion fixed in experiments.csv (position within max(5%, 2 x seed spread)): dt 0.001 is within 1.9% of the finest step (dt 0.0005), dt 0.002 is not (+11%); dx 5 and 10 are 8% and 12% below dx 2.5; the first division is within 30 min for dt and dx (415-440 against 430 / 440).
3. **MOT-REF is not converged: the motile regime is lost at dt 0.0005 and at dx 10** (the founder stays at x = -140 with speed 0 for the whole run; ATP-mode share 0.99), while dt 0.002 and dx 2.5 reach the wall. The reference (+126 um) sits on one side of a threshold: of the 16 REF runs of the block that are not the reference (dt 2, dx 2, D 2, cost 2, Km x boundary 8), 6 never move (dt 0.0005, dx 10, D x0.5, and the grid points Km x1 / boundary x0.5, Km x2 / x0.5, Km x2 / x1), 1 barely moves (cost x2: -105 um), 6 reach +123 to +143 um and 3 (boundary x2) stay near the middle of the domain and divide. The cause is not diagnosed (the single-cell run has no cap issue expected; the ATP surplus above the basal flux at the first positions has not been read from the series).
4. **Sensitivity (seed 50):** boundary glucose x2 (1.3 mM) removes the migration in both variants (founder at -7 to -13 um, 8-16 cells, first division 235-325 min: the cell grows where it is); boundary x0.5 with Km x0.5-2 gives the wall without division in GR (1 cell) and no migration or a short one in REF; the motility cost x0.5 gives 2 cells (division at 525 min) in REF and no change in GR, x2 slows REF (-105 um); D(glucose) x0.5 stops REF and shifts GR by -4%, x2 reaches the wall in both. The runs that keep the regime of the variant (reach the wall and, for GR, divide): GR 5 of the 12 sensitivity runs reach the wall and divide (9 reach the wall, 8 divide); REF 4 of 12 reach the wall.
5. Open: the diagnosis of item 3 (series of the speed, the ATP flux and the glucose at the founder cell at dt 0.0005 against 0.001); seeds 51 and 52 at dt 0.0005 and dx 10; intermediate dt (0.00075, 0.0014); FVA replay (objective by mode); the figure from MOT-REF / MOT-GR.

## Block 3: why MOT-REF loses its motion (2026-11-07; MN5 job 48711673, 2 min 5 s, 16 runs rc 0; `MOT-START-REF_*`, `MOT-DT-REF_dt0.00075/0.0014`)

**Mechanism.** The founder starts at x = -140 um, on a voxel face, next to the x_min wall (glucose Dirichlet 0; 0.65 mM at x_max). In the ATP-optimised solution (metabolic_mode 1) the ATP flux is 17.5 x the glucose uptake (10.37 = 17.5 x 0.59 mmol/gDW/h at 0.040 mM, Km 0.5, Vmax 8); the basal flux of 8.39 needs an uptake of 0.48 mmol/gDW/h, i.e. **at least about 0.032 mM of glucose in the voxel of the cell**; below that the solution is infeasible, the recorded ATP flux is 0, the motility flux is 0 and the cell does not move. Glucose in the founder voxel after 10 min, ATP flux (mmol/gDW/h) and outcome:

| start x (um) | dt 0.0005 | dt 0.001 | dt 0.002 |
|---|---|---|---|
| -142.5 | 0.0207 mM, 0 -> stays | 0.0298, 0 -> stays | 0.0462, 11.8 -> +142.6 |
| -140 (shipped; voxel face) | 0.0207, 0 -> stays | 0.040, 10.4 -> +126 | 0.056, 14.1 -> +140 |
| -137.5 | 0.0314, 0 -> stays | 0.0400, 10.4 -> +135.3 | 0.0559, 14.1 -> +140.8 (2 cells, 585 min) |
| -130 | 0.0526, 13.3 -> +132.0 (2 cells, 595 min) | 0.0607, 15.2 -> +136.8 (560) | 0.0754, 18.3 -> +141.0 (525) |
| -120 | 0.0739, 18.0 -> +130.4 (2 cells, 505) | 0.0813, 19.6 -> +136.7 (495) | 0.0948, 22.3 -> +138.9 (485) |

Every run is consistent with the threshold of 0.031-0.032 mM (0.0314 gives 0; 0.0400 gives 10.4). dt 0.0014 at -140: +135.5 um; dt 0.00075 at -140: -29.3 um (the glucose is at the threshold and the cell moves slowly). dx at -120: dx 2.5 +141.3 (2 cells, 490), dx 10 +130.0 (2 cells, 505).

**Why the glucose depends on dt.** Far from the wall the glucose is identical for all dt (0.330 mM at x = 2.5), near the Dirichlet-0 wall it is above the ideal steady profile (linear between the two boundary voxel centres) by about 15-19 x dt mM (first interior voxel at x = -142.5: +0.0097, +0.0187, +0.0354 mM at dt 0.0005, 0.001, 0.002; at -120: +0.008, +0.015, +0.029): the Dirichlet condition is enforced with an error that is first order in dt (the same one-voxel step as in the consortium, cross-feeding/RESULTS.md); the cause in the BioFVM scheme was not read from the code. As dt -> 0 the glucose at x = -140 tends to 0.011-0.022 mM, **below the threshold: in the limit the founder at x = -140 would not move**. The shipped motility-only demonstration therefore depends on the boundary error; MOT-GR is not affected (Km 0.02 mM: uptake 4 mmol/gDW/h at 0.02 mM).

**Convergence when the founder starts inside (x = -120 um).** The motile regime is present at every dt and dx tested; the first division is at 485-505 min (spread 20 min, criterion 30 min); the position at 10 h is +130.4 / +136.7 / +138.9 um at dt 0.0005 / 0.001 / 0.002 (dt 0.001 is 4.8% from the finest, within the 5% criterion; dt 0.002 is 6.5%) and +141.3 / +136.7 / +130.0 at dx 2.5 / 5 / 10 (dx 5 is 3.3% and dx 10 8.0% from the finest). At -130 the first division is 525-595 min (spread 70 min). The remaining dependence is the first-order boundary error (the glucose near the wall rises with dt).
Note: dt 0.00075 and 0.0014 do not divide the snapshot interval (5 min), so PhysiCell wrote only the initial and the final snapshot for these two runs; use dt that divide 5 min.

## MOT2-KM-REF_Km0.8_s50 (2026-11-07; MN5 job 48715446)
Motility-only variant with Km(glucose) 0.8 mM (the value of the manuscript), founder at x = -120 um, seed 50: ATP mode for the whole run, x = -114 / -95 / -58 / +6 / +105 / +135 um at 1 / 3 / 5 / 7 / 9 / 10 h, 1 cell, no division, mean speed 0.425 um/min; the mode switches to biomass at the last snapshot (600 min) on arrival. Reproduces migration without division inside 10 h, with the switch at the end of the window (Km 0.5: arrival about 6 h, division 495 min; Km 0.02: arrival about 4 h, division 400 min). Series in `results/mot2km_series.csv`.

## Block 6, MOT3-* (2026-11-07; MN5 job 48716078): regimes selected by the glucose supply (decision of the author: option 1)
Low supply (boundary 0.325 mM, below the switching threshold 0.35 mM), growth-recovery kinetics, founder at -120 um: the founder reaches the boundary at 235-240 min in every run (seeds 50-52, repeat, dt 0.0005 / 0.002, dx 2.5 / 10, pFBA), 1 cell, ATP mode 99%, no division. Km scan: arrival 225-255 min for Km 0.005-0.04 mM, 280 / 345 / 560 min at 0.08 / 0.16 / 0.32 mM, no movement at 0.64 and 1.28 mM (low supply); at the high supply 225-235 / 255 / 280 / 345 / 560 min for 0.005-0.04 / 0.08 / 0.16 / 0.32 / 0.64 mM and no movement at 1.28 mM, first division 395-495 min; the limit follows from the basal ATP flux (Km below about 15.7 x the start glucose). Details in `results/mot3_*.csv` and in the paper repo `Review_round_01/Motility_edits_plan.md` section 7.

## Block 7, MOT4-* (2026-11-08; MN5 job 48760763, 18 runs, rc 0, 3 min 44 s): high supply moved from 0.65 to 1.0 mM
Why: at 0.65 mM the glucose at the stopping position is 0.608 mM, 1.4% above the switching threshold (0.6 mM); the recorded mode of seed 50 stayed 1 (growth_rate 0) at the wall while seeds 51 and 52 switched. Same growth-recovery kinetics (Km 0.02, thresholds 0.35 / 0.6, cost 66.4), founder at -120 um, boundary 1.0 mM; design in `scripts/design_def.py` (block 7), `experiments.csv` (MOT4-*), analysis `results/mot4_{series,summary}.csv`.

| run | founder at 10 h (um) | cells at 10 h | first division (min) |
|---|---|---|---|
| MOT4-HIGH s50 (and repeat) / s51 / s52 | 37.6 / 37.9 / 36.7 | 8 / 8 / 8 | 290 / 285 / 295 |
| dt 0.0005 / 0.002 | 35.2 / 40.0 | 8 / 8 | 290 / 295 |
| dx 2.5 / 10 | 42.7 / 31.4 | 8 / 8 | 305 / 290 |
| pFBA | 37.6 (identical to s50) | 8 | 290 |

Km scan (time of the founder to x = 25 um, min): 150 (0.005-0.02 mM), 155, 165, 180, 210, 295 (0.04, 0.08, 0.16, 0.32, 0.64 mM), not reached at 1.28 and 2.56 mM; low supply (MOT3): 150, 155, 160, 175, 195, 250, 450 (0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32 mM), not reached at 0.64 and 1.28 mM. Cells at 10 h in the scan: 8 (Km <= 0.08), 16 (0.16), 10 (0.32), 4 (0.64), 1 (>= 1.28).
The cell stops at x = 30-31 um (the glucose reaches 0.6 mM at 30 um) at about 150-200 min, grows and divides at 285-295 min. Growth rate at the stop against Km (0.6 mM): 0.55 / 0.44 / 0.31 h-1 at 0.16 / 0.32 / 0.64 mM.

**Open: recorded mode at the threshold.** In MOT4-HIGH_s50 (and MOT2-GR_s50) the cell at 0.60-0.616 mM is recorded in mode 1 (growth_rate 0, speed 1, atp_flux 75) while its volume grows exactly like that of seed 51 (recorded mode 0, growth 0.683); after the first division the daughters are recorded in mode 1 with glucose above 0.6 mM, so population means of speed, growth and ATP flux show steps after the first division. Position, division time and cell number are not affected. Not traced to the code (custom.cpp `metabolic_bound_migration_rule`, inheritance of `metabolic_mode` and of `custom_optimization` by daughters). Consequence for the figure: panels A-F show seed 51 up to 280 min (before the first division).

**Fig. 7** (`results/fig7/`): `fig7_AF.png/svg` (`scripts/fig7_composite.py --scenario runs/MOT4-HIGH_s51::... --scenario runs/MOT3-LOW_s50::... --tmax-min 280`; `PYTHONPATH=<pctk>/src`), `panel_g.png/svg` (`scripts/fig7_panel_g.py`), `fig7_lineage.png` (lineage traces with forks, not used: daughters overlap the parent). The earlier version with 0.65 mM is superseded (MOT2-GR / MOT3-LOW with 0.65 mM kept as record).
