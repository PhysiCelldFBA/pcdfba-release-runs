# Metabolism-driven motility: results on the pinned code

Status 2026-11-06. Reference and seed runs done (8 runs); convergence (dt, dx), pFBA / FVA replay and the sensitivity block are not run yet.

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
